struct VisualizerState
    fuzzed_graph::FuzzedGraph
end

function visualizer()
    state = Observable(next_state())

    fig = Figure()
    ax = Axis(fig[1, 1], aspect=1)

    lines!(ax, @lift($state.fuzzed_graph.edges.geometry), color="black")
    scatter!(ax, @lift($state.fuzzed_graph.x), @lift($state.fuzzed_graph.y))

    # set up the button
    fig[2, 1] = iface = GridLayout(tellwidth=false)
    iface[1, 1] = nextbutton = Button(fig, label="Next graph")
    seed = @lift("Seed: $($state.fuzzed_graph.settings.seed)")
    Label(iface[2, 1], seed)

    on(nextbutton.clicks) do _
        @debug("New state")
        state[] = next_state()
    end

    fig
end

function next_state()
    VisualizerState(get_graph())
end
get_graph() = build_fuzzed_graph(FuzzedGraphSettings(seed=rand(UInt64)))