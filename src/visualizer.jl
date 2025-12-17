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

    # set up the ruler tool (h/t Google AI, but I wrote the code myself)
    ruler_p1 = Observable{Union{Point2f, Nothing}}(nothing)
    ruler_p2 = Observable{Union{Point2f, Nothing}}(nothing)
    cursorpos = Observable(Point2f(0, 0))

    on(events(fig).mousebutton) do event
        if event.button == Mouse.left && event.action == Mouse.press
            if !isnothing(ruler_p1[]) && isnothing(ruler_p2[])
                # finishing a measurement
                ruler_p2[] = mouseposition(ax.scene)
            else
                # starting a new measurement
                ruler_p1[] = mouseposition(ax.scene)
                ruler_p2[] = nothing
            end
        else
            cursorpos[] = mouseposition(ax.scene)
        end
    end

    sg = SliderGrid(
        iface[1, 1],
        (label="Density (int/km²)", range=10:1:100, startvalue=50)
    )

    distance = @lift(isnothing($ruler_p2) ? 0.0 : round(norm2($ruler_p1 - $ruler_p2), digits=1))
    distreadout = @lift("Distance: $($distance)m")
    Label(iface[2, 1], distreadout)

    ruler_line = @lift(if isnothing($ruler_p1)
        Point2f[]
    elseif isnothing($ruler_p2)
        [$ruler_p1]
    else
        [$ruler_p1, $ruler_p2]
    end)
    lines!(ax, ruler_line, color="pink")
    scatter!(ax, ruler_line, color="pink")

    seed = @lift("Seed: $($state.fuzzed_graph.settings.seed)")
    Label(iface[3, 1], seed)

    seedbox = Textbox(iface[4, 1], placeholder = "Set seed...", width=290, validator = UInt64)

    on(seedbox.stored_string) do s
        state[] = next_state(FuzzedGraphSettings(intersection_density = sg.sliders[1].value[] / (1000^2)), parse(UInt64, s))
    end

    iface[5, 1] = nextbutton = Button(fig, label="Random graph")

    on(nextbutton.clicks) do _
        @debug("New state with seed $(seedval[1])")
        state[] = next_state(FuzzedGraphSettings(intersection_density = sg.sliders[1].value[] / (1000^2)))
    end

    poslabel = @lift("Coordinate: $($cursorpos)")
    Label(iface[6, 1], poslabel)

    fig
end

function next_state(settings=FuzzedGraphSettings(), seed=nothing)
    settings.seed = isnothing(seed) ? rand(UInt64) : seed
    fuzzed = build_fuzzed_graph(settings)

    G = graph_from_gdal(fuzzed.edges)
    #G = remove_tiny_islands(G, 4)

    dmat = zeros(Float64, (nv(G), nv(G)))
    fill_distance_matrix!(G, dmat; maxdist=1000)

    all_links = identify_potential_missing_links(G, dmat, 100, 1000)
    links = deduplicate_links(G, all_links, dmat, 100)

    VisualizerState(
        fuzzed,
        isempty(all_links) ? DataFrame(geometry=[]) : links_to_gis(G, all_links),
        isempty(all_links) ? DataFrame(geometry=[]) : links_to_gis(G, links)
    )
end