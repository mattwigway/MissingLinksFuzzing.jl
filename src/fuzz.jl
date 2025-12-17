const MAX_GEOG_DIST = 100
const MIN_NET_DIST = 1000

"""
    fuzz(G, links)

Links should be the output of identify_missing_links, not deduplication;
since deduplication is a heuristic it is difficult to conclusively prove
correctness.
"""
function fuzz(G::FuzzedGraph, mlg, all_links, links, scores, oweights, dweights, logfile)
    if !isnothing(logfile)
        println(logfile, "seed,act_score,exp_score,difference,geom")
    end

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
                    for (i, link) in enumerate(all_links)
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
                        # this can happen when links have multiple points at the same distance (most often when they cross
                        # each other twice). If one of them got reversed in graph build the point found may be the same length but
                        # different location.
                        @warn "Link from $pt_e1 @ $pos_e1 / $len_e1 to $pt_e2 @ $pos_e2 / $len_e2 with length $geog_dist (network $net_dist) not found, seed $(G.settings.seed)"
                    end
                end
            end
        end
    end

    if !all(found_links)
        @warn "$(sum(.!found_links)) / $(length(found_links)) found but not expected, seed $(G.settings.seed)"
    end

    # col vector x row vector -> square matrix for all combinations
    odweights = reshape(oweights, (:, 1)) * reshape(dweights, (1, :))
    base_score = naive_access(mlg, [], odweights)
    for (link, score) in zip(links, scores)
        new_score = naive_access(mlg, [link], odweights) - base_score
        geom = links_to_gis(mlg, [link]).geometry[1]

        # write this out for later statistical analysis
        if !isnothing(logfile)
            println(logfile, "$(G.settings.seed),$score,$new_score,$(new_score - score),$(AG.toWKT(geom))")
        end

        # if the score from the algo was 0, they should be the same. if it wasn't,
        # allow 5% variation for rounding - because the algorithm uses lengths rounded to meters, whereas
        # the length of the link when realized is a float
        if score ≠ new_score && (score > 0 && abs(new_score / score - 1) > 0.05)
            link_gis = get_xy(geom)
            @warn "Link $(link_gis) has score $score from MissingLinks, expected $new_score; seed $(G.settings.seed)"
        end
    end

    # check the deduplication - every link found should have a link within 100 meters at both ends
    # this is a loose bound as it's using crow-flies not net distance
    link_xys = get_link_xy.(Ref(mlg), links)
    for link in all_links
        fr, to = get_link_xy(mlg, link)

        found = false

        for (fr2, to2) in link_xys
            if (AG.distance(fr, fr2) < 101 && AG.distance(to, to2) < 101) ||
                    # deduplication gets rid of opp-dir links as well
                    (AG.distance(to, fr2) < 101 && AG.distance(fr, to2) < 101)
                found = true
                break
            end
        end

        if !found
            @warn "Link from $fr to $to with length $(link.geographic_length_m) (network $(link.network_length_m)) removed in deduplication but no substitute found, seed $(G.settings.seed)"
        end

    end

    return(all_links)
end

function fuzzer(n; settings=FuzzedGraphSettings(), seed=nothing, logfile=nothing)
    # if seed is specified, use it to seed the RNG, but also use it as the seed for the first fuzz so things can be reproduced
    rng = StableRNG(isnothing(seed) ? rand(UInt64) : seed)
    seed = isnothing(seed) ? rand(rng, UInt64) : seed

    if !isnothing(logfile)
        open(f -> _fuzzer(rng, n, settings, seed, f), logfile, "w")
    else
        _fuzzer(rng, n, settings, seed, nothing)
    end
end

function _fuzzer(rng, n, settings, seed, logfile)
    @showprogress for _ in 1:n
        settings.seed = seed 
        fuzzed = build_fuzzed_graph(settings)

        local all_links, links, oweights, dweights, scores, f, G

        with_logger(NullLogger()) do
            G = graph_from_gdal(fuzzed.edges, max_edge_length=100_000)
            #G = remove_tiny_islands(G, 4) # TODO test this

            dmat = zeros(Float64, (nv(G), nv(G)))
            fill_distance_matrix!(G, dmat; maxdist=1000)

            all_links = identify_potential_missing_links(G, dmat, 100, 1000)
            links = deduplicate_links(G, all_links, dmat, 100)

            oweights = rand(rng, Float64, nv(G))
            dweights = rand(rng, Float64, nv(G))

            scores = score_links(x -> x < 1000, G, links, dmat, oweights, dweights, 1000) 
        end


        fuzz(fuzzed, G, all_links, links, scores, oweights, dweights, logfile)

        # reset seed for next iteration
        seed = rand(rng, UInt64)
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