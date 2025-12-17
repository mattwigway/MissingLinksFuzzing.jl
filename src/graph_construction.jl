@kwdef mutable struct FuzzedGraphSettings
    width::Int64 = 6000
    height::Int64 = 3500

    # density per square meter
    intersection_density::Float64 = 5e-5

    # connectivity factor: what proportion of nodes we build SPTs from with reweighted graph
    connectivity_factor::Float64 = 0.2

    # detour factor: how much MST/already added edges get weighted down
    detour_factor::Float64 = 0.25

    # max_length: max edge length
    max_length::Float64 = 100

    # The angle used is the angle towards the point, plus a draw from a normal distribution
    # with this standard deviation
    angle_sd_degrees::Float64 = 15.0

    # After the first point, subsequent points choose an angle towards the destination
    # based on both the angle of the previous segment and the angle towards the destination
    # This determines bias towards destination. 1 means center the angle choice at the angle
    # towards the destination, 0 means continue in the chosen direction. Setting to zero will
    # likely result in improbable paths.
    destination_bias::Float64 = 0.5

    # Mean length of a segment of a curved edge
    segment_length_mean::Float64 = 10

    # Standard deviation of a segment of a curved edge
    segment_length_sd::Float64 = 15

    # Fraction of edges to randomly delete
    delete_fraction::Float64 = 0.3

    seed::UInt64 = 42
end

struct FuzzedGraph
    settings::FuzzedGraphSettings
    graph::Graph
    edges::DataFrame
    x::Vector{Float64}
    y::Vector{Float64}
end

function build_fuzzed_graph(settings::FuzzedGraphSettings)
    # everything is built based on this one StableRNG so it can be reproduced
    rng = StableRNG(settings.seed)

    # generate the nodes
    n_nodes = round(Int64, settings.width * settings.height * settings.intersection_density)
    x = rand(rng,  Uniform(0, settings.width), n_nodes)
    y = rand(rng,  Uniform(0, settings.height), n_nodes)

    geom_dmat = Matrix{Float64}(undef, n_nodes, n_nodes)
    for i in 1:n_nodes
        for j in i:n_nodes
            geom_dmat[i, j] = geom_dmat[j, i] = norm2((x[i] - x[j], y[i] - y[j]))
        end
    end

    # base graph: fully connected
    G = complete_graph(n_nodes)

    @debug "Complete graph has $(ne(G)) edges"

    # base network: MST
    mst = prim_mst(G, geom_dmat)

    # keep track of which edges go in the final network
    edges_in_final_network = BitMatrix(undef, (n_nodes, n_nodes))
    fill!(edges_in_final_network, false)

    # weight down distances from MST
    for edge in mst
        edges_in_final_network[edge.src, edge.dst] = edges_in_final_network[edge.dst, edge.src] = true
        new_dist = geom_dmat[edge.src, edge.dst] * settings.detour_factor
        geom_dmat[edge.src, edge.dst] = geom_dmat[edge.dst, edge.src] = new_dist
    end

    # ÷ 2 b/c undirected
    @debug "MST graph has $(sum(edges_in_final_network) ÷ 2) edges"

    # sample nodes for SPTs
    for node in sample(rng, 1:n_nodes, round(Int64, n_nodes * settings.connectivity_factor), replace=false)
        spt = dijkstra_shortest_paths(G, node, geom_dmat)

        for dest in eachindex(spt.parents)
            current = dest
            while current ≠ node
                prev = spt.parents[current]
                if !edges_in_final_network[prev, current]
                    # this edge is not yet in the final network but is in this SPT so should be added to final
                    # network and weighted down for already existing so future paths prefer to use it.
                    edges_in_final_network[prev, current] = edges_in_final_network[current, prev] = true
                    new_dist = geom_dmat[prev, current] * settings.detour_factor
                    geom_dmat[prev, current] = geom_dmat[current, prev] = new_dist
                end
                current = prev
            end
        end
    end

    @debug "SPT graph has $(sum(edges_in_final_network) ÷ 2) edges"

    # remove unused edges
    for i in 1:n_nodes
        for j in i:n_nodes
            if !edges_in_final_network[i, j]
                rem_edge!(G, i, j)
            end
        end
    end

    # randomly remove edges
    edges_to_remove = sample(rng, collect(edges(G)), round(Int64, ne(G) * settings.delete_fraction))

    for edge in edges_to_remove
        rem_edge!(G, edge)
    end

    @debug "Final graph has $(nv(G)) vertices, $(ne(G)) edges and $(length(connected_components(G))) components"

    # turn into GIS dataset
    return FuzzedGraph(
        settings,
        G,
        build_edges(settings, rng, G, x, y),
        x,
        y
    )
end

function build_edges(settings, rng, G, x, y)
    geoms = AG.IGeometry{AG.wkbLineString}[]
    sources = Int64[]
    targets = Int64[]

    for edge in edges(G)
        push!(geoms, build_edge(settings, rng, x[edge.src], y[edge.src], x[edge.dst], y[edge.dst]))
        push!(sources, edge.src)
        push!(targets, edge.dst)
    end

    gdf = DataFrame(fid=1:length(geoms), geometry=geoms, src=sources, dst=targets)
    metadata!(gdf, "geometrycolumns", (:geometry,))
    return gdf
end

function build_edge(settings::FuzzedGraphSettings, rng, src_x, src_y, dst_x, dst_y)
    coords = [[src_x, src_y]]

    seg_length_dist = truncated(Normal(settings.segment_length_mean, settings.segment_length_sd), settings.segment_length_mean * 0.05, settings.segment_length_mean * 3)
    angle_dist = Normal(0, settings.angle_sd_degrees)

    prev_dist = norm2([dst_x, dst_y] .- [src_x, src_y])
    prev_angle = atand(dst_x - src_x, dst_y - src_y)

    while true
        # calculate angle to destination
        angle_to_dest = atand(dst_y - coords[end][2], dst_x - coords[end][1])
        angle = midpoint_angle_degrees(prev_angle, angle_to_dest, settings.destination_bias)

        # fudge to make road curvey
        # TODO this could cause self intersections in extreme cases
        angle += rand(rng, angle_dist)

        # how far we go
        seg_len = rand(rng, seg_length_dist)

        # offset point
        candidate_point = coords[end] .+ [cosd(angle) * seg_len, sind(angle) * seg_len]

        new_dist = norm2([dst_x, dst_y] .- candidate_point)

        # if new_dist > orig_dist, this point is past the destination - just go straight there instead
        # and stop the algorithm
        if new_dist < settings.segment_length_mean * 2
            push!(coords, [dst_x, dst_y])
            break
        else
            push!(coords, candidate_point)
            prev_angle = angle
            prev_dist = new_dist
        end
    end

    geom = AG.createlinestring(coords)

    if AG.isvalid(geom)
        return geom
    else
        # try again, we may have a self-intersection etc
        return build_edge(settings, rng, src_x, src_y, dst_x, dst_y)
    end
end

function midpoint_angle_degrees(a1, a2, bias=0.5)
    # make sure a1 is larger
    if a2 > a1
        a1, a2 = a2, a1
        bias = 1 - bias
    end

    if a1 - a2 > 180
        # angle should go other way
        a2 += 360
    end

    (a2 * bias + a1 * (1 - bias)) % 360
end

# function Plots.plot(g::FuzzedGraph)
#     Plots.plot(g.edges.geometry, color="black", aspect_ratio=:equal)
#     Plots.scatter!(g.x, g.y, color="#4B9CD3", markersize=1.5, markerstrokewidth=0, legend=false)
# end