const MAX_GEOG_DIST = 100
const MIN_NET_DIST = 1000

"""
    fuzz(G, links)

Links should be the output of identify_missing_links, not deduplication;
since deduplication is a heuristic it is difficult to conclusively prove
correctness.
"""
function fuzz(G::FuzzedGraph, mlg, all_links, links, scores, oweights, dweights, logfile)
    # make a weight matrix accounting for the now squiggly roads
    weights = fill(Inf64, length(G.x), length(G.x))
    for row in eachrow(G.edges)
        weights[row.src, row.dst] = weights[row.dst, row.src] = AG.geomlength(row.geometry)
    end

    dmat = fill(Inf64, (length(G.x), length(G.x)))
    Threads.@threads for src in 1:length(G.x)
        dmat[:, src] = dijkstra_shortest_paths(G.graph, src, weights).dists
    end
    
    found_links = BitVector(undef, length(all_links))
    fill!(found_links, false)

    # if two edges cross more than once, which location gets the link is undefined behavior.
    # this is unusual and even more so in the real world, so we just punt here, but in those rare
    # cases we expect to not find those edges (e.g. seed 7122159617084383647)
    expected_not_found = 0

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
                # because of rounding of offsets, the ends of the links may move a little before network distance
                # calculations, so allow a 1-meter buffer
                if round(UInt16, min(net_dist, typemax(UInt16))) > MIN_NET_DIST + 1
                    found_this_link = false
                    best_link = nothing
                    best_link_dist = Inf64

                    # if these two lines cross more than once, all bets are off - depending on whether any geometries got reversed,
                    # we might find a different one, and which one we find is undefined behavior. This is rare, so just punt in
                    # this case.
                    # seed 7122159617084383647 demonstrates this at (1523 456)
                    if geog_dist ≈ 0.0 && AG.ngeom(AG.intersection(e1.geometry, e2.geometry)) > 1
                        expected_not_found += 1 # not + 2, this link should exist in both directions but we will also find it in both directions here
                        continue
                    end

                    # These should have a link (unless they were part of an island - todo)
                    # TODO no shared IDs here what to do
                    for (i, link) in enumerate(all_links)
                        link_geom = map(pt -> round.(pt), get_xy(links_to_gis(mlg, [link]).geometry[1]))
                        geomdiff = abs.(Iterators.flatten(link_geom) .- [LibGEOS.getGeomX(pt_e1), LibGEOS.getGeomY(pt_e1), LibGEOS.getGeomX(pt_e2), LibGEOS.getGeomY(pt_e2)])

                        fromdiff = min(
                                max(
                                    abs(link.fr_dist_from_start - round(UInt16, pos_e1)) <= 2,
                                    abs(link.fr_dist_to_end - (round(UInt16, len_e1) - round(UInt16, pos_e1))) <= 2
                                ), max(
                                    # edge could have been reversed
                                    abs(link.fr_dist_to_end - round(UInt16, pos_e1)) <= 2,
                                    abs(link.fr_dist_from_start - (round(UInt16, len_e1) - round(UInt16, pos_e1))) <= 2
                                )
                        )

                        todiff = min(
                                max(
                                    abs(link.to_dist_from_start - round(UInt16, pos_e2)),
                                    abs(link.to_dist_to_end - (round(UInt16, len_e2) - round(UInt16, pos_e2)))
                                ), max(
                                    # edge could have been reversed
                                    abs(link.to_dist_to_end - round(UInt16, pos_e2)),
                                    abs(link.to_dist_from_start - (round(UInt16, len_e2) - round(UInt16, pos_e2)))
                                )
                        )

                        geogdiff = abs(link.geographic_length_m - round(UInt16, geog_dist))

                        this_link_diff = sum([geomdiff..., fromdiff, todiff, geogdiff]) 

                        if fromdiff <= 2 && todiff <= 2 && geogdiff <= 2 &&
                        # checking the geometry
                        all(abs.(geomdiff) .< 4) &&
                        !found_links[i] &&
                        this_link_diff < best_link_dist
                            # this link matches, but see if there is one that matches better before calling it a match,
                            # so we don't accidentally match things we shouldn't
                            best_link = i
                            best_link_dist = this_link_diff
                            found_this_link = true # we found a link that could be this link anyways, and if there isn't one closer it must be it
                        end
                    end

                    if found_this_link
                        found_links[best_link] = true
                    else
                        @warn "Link from $pt_e1 @ $pos_e1 / $len_e1 to $pt_e2 @ $pos_e2 / $len_e2 with length $geog_dist (network $net_dist) not found, seed $(G.settings.seed)"
                    end
                end
            end
        end
    end

    for (i, link) in enumerate(all_links)
        if link.network_length_m <= MIN_NET_DIST + 1
            # if very close to 1km we might think we should have it but not have it due to rounding,
            # so just say we found it regardless
            found_links[i] = true
        end
    end



    if sum(.!found_links) ≠ expected_not_found
        @warn "$(sum(.!found_links) - expected_not_found) / $(length(found_links)) found but not expected, seed $(G.settings.seed):" get_xy.(links_to_gis(mlg, all_links[.!found_links]).geometry)
    end

    # col vector x row vector -> square matrix for all combinations
    odweights = reshape(oweights, (:, 1)) * reshape(dweights, (1, :))
    base_score = naive_access(mlg, [], odweights)
    for (link, score) in zip(links, scores)
        new_score = naive_access(mlg, [link], odweights) - base_score
        geom = links_to_gis(mlg, [link]).geometry[1]

        # write this out for later statistical analysis
        # It is expected they will not all match perfectly; the algorithms differ in how they handle rounding.
        # In the MissingLinks algorithm, rouding happens to multiple components individually - distance to the link,
        # distance of the link, distance from the link to the destination. In the naive algorithm, nothing gets rounded
        # until the eend (other than the offsets, see below).

        # Note that it is expected that this will introduce a _very slight_ negative bias to the results. Julia
        # uses round halves to evens (bankers rounding), so in general rounding does not introduce any systematic
        # bias. However, in this case it results in the distances being (on average) marginally longer in the
        # realized graph than as calculated by the MissingLinks algorithm. I say on average because this small
        # bias is far smaller than the overall rounding error that could be positive or negative.

        # Several things are getting rounded in the algorithm - the distance matrix values, the link offsets
        # from the start and end of the link, and the length of the link itself. All of these use default rounding,
        # so none of them are biased. However, consider a path from point A to B via link L. The MissingLinks algorithm
        # calculates this as the distance from A to the end of the edge L connects to (it could be either end, it calculates
        # which is shortest), the distance from the end of L to where L connects, the length of L itself, the distance from
        # where L connect to the end of the other edge it is connected to, and the distance from the end of that edge to the
        # destination. These are all rounded independently. For the most part, all of that rounding should cancel out on average,
        # leading to unbiased results. There is one exception, however. When we round the offsets from the ends of the edges to L,
        # we may make the distance to get to L further from or closer to the origin/destination, with no bias (and the
        # algorithm ensures that the total of the offsets from each edge sum to the link length). However, rounding
        # these offsets _does_ bias the length of L itself very slightly. Since L is always at the closest point between
        # the edges, rounding the offsets will always make it longer, because moving in either direction moves away from the
        # optimal point, introducing up to 1m in additional length (if the true offsets were both x.5, the links were colinear,
        # and rounding moved the ends farther apart). The MissingLinks algorithm doesn't account for this; it treats the length
        # of the link as the geographic distance at the shortest point, rounded. But when we realize the graph,
        # the actual length from the rounded point is used. So the links in the realized graph are ever-so-slightly longer
        # leading to longer distances and lower access on average. The magnitude of this bias is so tiny it falls several orders
        # of magnitudes below many other sources of error in the algorithm, so really not worth worrying about - and in fact makes
        # the result probably better, as the shortest point between two edges is unlikely to be a place you could exactly build
        # a connection anyways. But it is statistically significant in large samples.
        if !isnothing(logfile)
            println(logfile, "$(G.settings.seed),$score,$new_score,$(new_score - score),\"$(AG.toWKT(geom))\"")
        elseif score ≠ new_score && (score > 0 && abs(new_score / score - 1) > 0.05)
            # if the score from the algo was 0, they should be the same. if it wasn't,
            # allow 5% variation for rounding - because the algorithm uses lengths rounded to meters, whereas
            # the length of the link when realized is a float. But if there's a logfile, don't spam the console;
            # any concerns can be identified post-hoc.
            link_gis = get_xy(geom)
            @warn "Link $(link_gis) has score $score from MissingLinks, expected $new_score; seed $(G.settings.seed)"
        end
    end

    # new approach to checking deduplication - realize graph with all links and with deduplicated links
    # create distance matrices for both. no trip should be longer due to deduplication by more than
    # 400 * number of links used in the original. Why 400? Each sphere of influence has a radius of
    # 100m around each end of the link that defined it, and then the shortest link is selected. If the
    # link in the shortest path was 100m in one direction from both ends of the defining link, and the
    # retained link 100m in the other direction, each end moves 200m. The link itself has to be shorter, but
    # maybe not by much (in the case of two parallel streets, could be by millimeters)

    # silence warnings about duplicate links, this is expected
    Gall = with_logger(() -> realize_graph(mlg, all_links), NullLogger())

    # _not_ nv(Gall) - realize_graph adds vertices at the end and we don't care what's going on with
    # added vertices
    dmatall = Matrix{Float64}(undef, nv(mlg), nv(mlg))
    nlinksused = zeros(Int64, nv(mlg), nv(mlg))
    for src in 1:nv(mlg)
        spt = dijkstra_shortest_paths(Gall, src)
        dmatall[:, src] = spt.dists[1:nv(mlg)]

        for dst in 1:nv(mlg)
            if isfinite(spt.dists[dst])
                node = dst
                while node != src
                    parent = spt.parents[node]
                    ltype = Gall[label_for(Gall, parent), label_for(Gall, node)].link_type
                    if !ismissing(ltype) && ltype == "candidate"
                        nlinksused[src, dst] += 1
                    end
                    node = parent
                end
            end
        end
    end

    Gdedupe = with_logger(() -> realize_graph(mlg, links), NullLogger())
    dmatdedupe = Matrix{Float64}(undef, nv(mlg), nv(mlg))
    for src in 1:nv(mlg)
        spt = dijkstra_shortest_paths(Gdedupe, src)
        dmatdedupe[:, src] = spt.dists[1:nv(mlg)]
    end

    # dedupe should make all trips the same or longer
    Δdmat = dmatdedupe .- dmatall
    # things that are unreachable should stay unreachable (since we don't have a distance limit)
    @assert all(isfinite.(dmatdedupe) .== isfinite.(dmatall))
    # get rid of NaNs (inf - inf == NaN)
    Δdmat[.!isfinite.(dmatdedupe)] .= 0.0
    all(Δdmat .> -1e-6) || @error "Dedupe dmat should always be same or larger than full dmat, but got values $(Δdmat[Δdmat .<= -1e-6]), seed $(G.settings.seed)"
    dedupe_problems = Δdmat .> nlinksused .* 400 .+ 1e-6

    dedupe_deets = DataFrame(full=reshape(dmatall[dedupe_problems], :), dedupe=reshape(dmatdedupe[dedupe_problems], :), n_links_used=reshape(nlinksused[dedupe_problems], :))

    !any(dedupe_problems) || @warn "deduplication caused $(sum(dedupe_problems)) paths to get too much longer, seed $(G.settings.seed)" dedupe_deets

    return(all_links)
end

function fuzzer(;settings=FuzzedGraphSettings(), seed=nothing, logfile=nothing, csvfile=nothing, once=false, headless=false)
    # if seed is specified, use it to seed the RNG, but also use it as the seed for the first fuzz so things can be reproduced
    rng = StableRNG(isnothing(seed) ? rand(UInt64) : seed)
    seed = isnothing(seed) ? rand(rng, UInt64) : seed

    if !isnothing(logfile) && logfile != ""
        logfile = replace(logfile, "%i" => "$seed")
        Base.global_logger(MinLevelLogger(TeeLogger(ConsoleLogger(), FileLogger(logfile)), Info))
    end

    if !isnothing(csvfile) && csvfile != ""
        csvfile = replace(csvfile, "%i" => "$seed")
        open(f -> _fuzzer(rng, settings, seed, f, once, headless), csvfile, "w")
    else
        _fuzzer(rng, settings, seed, nothing, once, headless)
    end
end

function _fuzzer(rng, settings, seed, csvfile, once, headless)
    if !isnothing(csvfile)
        println(csvfile, "seed,act_score,exp_score,difference,geom")
    end

    # in headless mode, don't use console progress bar
    progress = headless ? ProgressUnknownHeadless(10, "Fuzzed graphs tested:") : ProgressUnknown("Fuzzed graphs tested:")

    while true
        settings.seed = seed
        try
            fuzzed = build_fuzzed_graph(settings)

            local all_links, links, oweights_exp, dweights_exp, scores, G

            with_logger(NullLogger()) do
                G = graph_from_gdal(fuzzed.edges, max_edge_length=100_000)
                #G = remove_tiny_islands(G, 4) # TODO test this

                dmat = zeros(Float64, (nv(G), nv(G)))
                fill_distance_matrix!(G, dmat; maxdist=1000)

                all_links = identify_potential_missing_links(G, dmat, 100, 1000)
                links = deduplicate_links(G, all_links, dmat, 100)

                oweights_exp, oweightsdf = get_weights(rng, fuzzed, G, false)
                dweights_exp, dweightsdf = get_weights(rng, fuzzed, G, true)

                # We de-construct and re-construct the weights
                oweights_act = create_graph_weights(G, oweightsdf, [:weight], 1e-5)
                dweights_act = create_graph_weights(G, dweightsdf, [:weight], 1e-5)

                all(oweights_act .≈ oweights_exp) || @warn "Origin weights deviate from expected by up to $(max(abs(oweights_act .- oweights_exp))), seed $(seed)"
                all(dweights_act .≈ dweights_exp) || @warn "Origin weights deviate from expected by up to $(max(abs(dweights_act .- dweights_exp))), seed $(seed)"

                scores = score_links(x -> x < 1000, G, links, dmat, oweights_act, dweights_act, 1000) 
            end


            fuzz(fuzzed, G, all_links, links, scores, oweights_exp, dweights_exp, csvfile)

            if !isnothing(csvfile)
                flush(csvfile)
            end
        
        catch ex
            if ex isa InterruptException
                rethrow(e) # exit
            else
                @error "Error occurred with seed $seed" exception=(ex, catch_backtrace())
            end
        end

        if once
            break
        end

        # reset seed for next iteration
        seed = rand(rng, UInt64)
        next!(progress)
    end
end

# build a distance matrix
function build_uint16_dmat(G)
    result = Matrix{UInt16}(undef, nv(G), nv(G))

    for src in 1:nv(G)
        result[:, src] = round.(UInt16, min.(dijkstra_shortest_paths(G, src).dists, typemax(UInt16)))
    end

    return result
end

function naive_access(G, links, odweights)
    # splice the link into the graph and recalculate the distance matrix
    Gnew = isempty(links) ? G : realize_graph(G, links)
    # any newly-added vertices don't matter and will have weight 0.
    dmat = build_uint16_dmat(Gnew)[1:size(odweights, 1), 1:size(odweights, 2)]
    # filter based on distance
    sum(odweights .* (dmat .< 1000))
end

function get_link_xy(G, link)
    e1 = G[link.fr_edge_src, link.fr_edge_tgt]
    p1 = AG.pointalongline(e1.geom, link.fr_dist_from_start)

    e2 = G[link.to_edge_src, link.to_edge_tgt]
    p2 = AG.pointalongline(e2.geom, link.to_dist_from_start)

    p1, p2
end

"""
    get_weights(rng, fuzzed, G, poly)

Get random weights for a fuzzed graph. Return the expected weights and a geodataframe
that should produce them when run through the assignment procedure. If poly is true, return
polygons.
"""
function get_weights(rng, fuzzed, G, poly)
    # We create a point at each vertex of the graph; it should get snapped
    # to all edges connecting to that vertex
    # adding 1 so that weights can't be infinitesimal, we can easily see when something is wrong
    original = rand(rng, Float64, length(fuzzed.x)) .+ 1

    # create the geodataframe (collect converts tuple to vector for AG)
    df = DataFrame(weight=original, geometry=map(AG.createpoint ∘ collect, zip(fuzzed.x, fuzzed.y)))

    if poly
        # make them polygons
        df.geometry = AG.buffer.(df.geometry, 1e-2)
    end

    metadata!(df, "geometrycolumns", (:geometry,))

    # now, based on the graph structure, figure out what the computed weights should be
    # for each point in the dataframe, it will be divided among all edges connecting to the node
    # (assuming there aren't other edges extremely close, but that's unlikely). Then those edges
    # will be divided among the node and the neighbor, so half will go to this node, half evenly
    # divided among neighbors.

    # We do this first in fuzzed-graph space, then convert to MissingLinks-graph space

    total_expected_weight = zero(Float64)
    weights_fuzzed = zeros(Float64, length(fuzzed.x))
    for (i, weight) in enumerate(original)
        nbrs = neighbors(fuzzed.graph, i)

        # nodes with no neighbors will not end up in the missing links graph and as such will
        # have no weight
        if !isempty(nbrs)
            total_expected_weight += weight
            weights_fuzzed[i] += weight / 2
            for nbr in nbrs
                weights_fuzzed[nbr] += weight / 2 / length(nbrs)
            end
        end
    end

    @assert sum(weights_fuzzed) ≈ total_expected_weight

    # translate to missinglinks-graph space
    weights_mlg = zeros(Float64, nv(G))
    for i in 1:nv(G)
        pos = G[label_for(G, i)]
        found = false
        for (j, x, y) in zip(1:length(fuzzed.x), fuzzed.x, fuzzed.y)
            if norm2([x, y] .- pos) < 1e-5
                found = true
                @assert !isnan(weights_fuzzed[j])
                weights_mlg[i] = weights_fuzzed[j]
                # so we don't use it again
                weights_fuzzed[j] = NaN64
                break
            end
        end
        @assert found
    end

    @assert sum(weights_mlg) ≈ total_expected_weight

    return weights_mlg, df
end