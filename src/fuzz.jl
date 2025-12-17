const MAX_GEOG_DIST = 100
const MIN_NET_DIST = 1000

"""
    fuzz(G, links)

Links should be the output of identify_missing_links, not deduplication;
since deduplication is a heuristic it is difficult to conclusively prove
correctness.
"""
function fuzz(G::FuzzedGraph, mlg, links)
    # make a weight matrix accounting for the now squiggly roads
    weights = fill(Inf64, length(G.x), length(G.x))
    for row in eachrow(G.edges)
        weights[row.src, row.dst] = weights[row.dst, row.src] = AG.geomlength(row.geometry)
    end

    dmat = fill(Inf64, (length(G.x), length(G.x)))
    Threads.@threads for src in 1:length(G.x)
        dmat[:, src] = dijkstra_shortest_paths(G.graph, src, weights).dists
    end
    
    found_links = BitVector(undef, length(links))
    fill!(found_links, false)

    # now naively check every pair of edges to see if they should have a link, and if so if they do
    for e1 in eachrow(G.edges)
        for e2 in eachrow(G.edges)
            geog_dist = AG.distance(e1.geometry, e2.geometry)
            if geog_dist < MAX_GEOG_DIST
                # these potentially should have a link
                # check net distance from closest points
                pt_e1, pt_e2 = LibGEOS.nearestPoints(e1.geometry, e2.geometry)

                # round everything - this accumulates a small amount of rounding error,
                # which we deem acceptable, and this way rounding error accumulates the same
                # way in fuzzing as in the main code.
                len_e1 = round(UInt16, AG.geomlength(e1.geometry))
                pos_e1 = round(UInt16, LibGEOS.project(e1.geometry, pt_e1))
                len_e2 = round(UInt16, AG.geomlength(e2.geometry))
                pos_e2 = round(UInt16, LibGEOS.project(e2.geometry, pt_e2))

                net_dist = min(
                    pos_e1 + dmat[e1.src, e2.src] + pos_e2,
                    pos_e1 + dmat[e1.src, e2.dst] + (len_e2 - pos_e2),
                    (len_e1 - pos_e1) + dmat[e1.dst, e2.src] + pos_e2,
                    (len_e1 - pos_e1) + dmat[e1.dst, e2.dst] + (len_e2 - pos_e2)
                )

                # if the net dist is actually say 1000.25m, we want to consider not greater than 1000m,
                # because distances all get rounded inside the missinglinks tool.
                if round(UInt16, min(net_dist, typemax(UInt16))) > MIN_NET_DIST
                    found_this_link = false

                    # These should have a link (unless they were part of an island - todo)
                    # TODO no shared IDs here what to do
                    for (i, link) in enumerate(links)
                        link_geom = map(pt -> round.(pt), get_xy(links_to_gis(mlg, [link]).geometry[1]))

                        if (
                           (abs(link.fr_dist_from_start - round(UInt16, pos_e1)) <= 2 && abs(link.fr_dist_to_end - (round(UInt16, len_e1) - round(UInt16, pos_e1))) <= 2 ) ||
                           # link could have been reversed
                           (abs(link.fr_dist_to_end - round(UInt16, pos_e1)) <= 2 && abs(link.fr_dist_from_start - (round(UInt16, len_e1) - round(UInt16, pos_e1))) <= 2 )
                        ) && (
                            (abs(link.to_dist_from_start - round(UInt16, pos_e2)) <= 2 && abs(link.to_dist_to_end - (round(UInt16, len_e2) - round(UInt16, pos_e2))) <= 2 ) ||
                            # link could have been reversed
                            (abs(link.to_dist_to_end - round(UInt16, pos_e2)) <= 2 && abs(link.to_dist_from_start - (round(UInt16, len_e2) - round(UInt16, pos_e2))) <= 2 )
                        ) &&
                        abs(link.geographic_length_m - round(UInt16, geog_dist) <= 2) &&

                        # checking the geometry
                        all(abs.(
                            Iterators.flatten(link_geom) .-
                            [LibGEOS.getGeomX(pt_e1), LibGEOS.getGeomY(pt_e1), LibGEOS.getGeomX(pt_e2), LibGEOS.getGeomY(pt_e2)]
                         ) .< 4) &&
                            !found_links[i]
                            found_links[i] = true
                            found_this_link = true
                            break
                        end
                    end

                    if !found_this_link
                        @warn "Link from $pt_e1 @ $pos_e1 / $len_e1 to $pt_e2 @ $pos_e2 / $len_e2 with length $geog_dist (network $net_dist) not found, seed $(G.settings.seed)"
                    end
                end
            end
        end
    end

    if !all(found_links)
        @warn "$(sum(.!found_links)) / $(length(found_links)) found but not expected, seed $(G.settings.seed)"
    end

    return(links)
end

function fuzzer(n; settings=FuzzedGraphSettings(), seed=nothing)
    # if seed is specified, use it to seed the RNG, but also use it as the seed for the first fuzz so things can be reproduced
    rng = StableRNG(isnothing(seed) ? rand(UInt64) : seed)
    seed = isnothing(seed) ? rand(rng, UInt64) : seed
    f = nothing
    G = nothing
    @showprogress for _ in 1:n
        settings.seed = seed 
        fuzzed = build_fuzzed_graph(settings)

        all_links = with_logger(NullLogger()) do
            G = graph_from_gdal(fuzzed.edges, max_edge_length=100_000)
            #G = remove_tiny_islands(G, 4) # TODO test this

            dmat = zeros(Float64, (nv(G), nv(G)))
            fill_distance_matrix!(G, dmat; maxdist=1000)

            identify_potential_missing_links(G, dmat, 100, 1000)
        end

        f = fuzz(fuzzed, G, all_links)

        # reset seed for next iteration
        seed = rand(rng, UInt64)
    end

    return (links=f, G=G)
end