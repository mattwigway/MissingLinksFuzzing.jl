struct VisualizerState
    fuzzed_graph::FuzzedGraph
    all_links::DataFrame
    links::DataFrame
end

function visualizer()
    state = Observable(next_state())

    fig = Figure()
    ax = Axis(fig[1, 1], aspect=@lift($state.fuzzed_graph.settings.width / $state.fuzzed_graph.settings.height))

    lines!(ax, @lift($state.fuzzed_graph.edges.geometry), color="black")
    scatter!(ax, @lift($state.fuzzed_graph.x), @lift($state.fuzzed_graph.y))

    # add all links in gray
    lines!(ax, @lift($state.all_links.geometry), color="gray")

    # and chosen in red
    lines!(ax, @lift($state.links.geometry), color="red")

    # set up the button
    fig[1, 2] = iface = GridLayout()
    colsize!(fig.layout, 1, Auto())
    colsize!(fig.layout, 2, Fixed(300))
    rowsize!(fig.layout, 1, Relative(1))

    sg = SliderGrid(
        iface[1, 1],
        (label="Density (int/km²)", range=10:1:100, startvalue=50)
    )

    seed = @lift("Seed: $($state.fuzzed_graph.settings.seed)")
    Label(iface[2, 1], seed)
    iface[3, 1] = nextbutton = Button(fig, label="Next graph")


    on(nextbutton.clicks) do _
        @debug("New state")
        state[] = next_state(FuzzedGraphSettings(intersection_density = sg.sliders[1].value[] / (1000^2)))
    end

    fig
end

function next_state(settings=FuzzedGraphSettings())
    settings.seed = rand(UInt64)
    settings.width = 6000
    settings.height = 3500
    fuzzed = build_fuzzed_graph(settings)

    G = graph_from_gdal(fuzzed.edges)

    dmat = zeros(Float64, (nv(G), nv(G)))
    fill_distance_matrix!(G, dmat; maxdist=1000)

    all_links = identify_potential_missing_links(G, dmat, 100, 1000)
    links = deduplicate_links(G, all_links, dmat, 100)

    VisualizerState(
        fuzzed,
        links_to_gis(G, all_links),
        links_to_gis(G, links)
    )
end