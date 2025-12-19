struct VisualizerState
    fuzzed_graph::FuzzedGraph
    all_links::DataFrame
    links::DataFrame
end

function visualizer()
    # we do this here to avoid loading the OpenGL ecosystem and opening a Makie window
    # when running the fuzzer
    @eval import GLMakie

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

    linkinfo = Observable("")

    on(events(fig).mousebutton) do event
        if event.button == Mouse.left && event.action == Mouse.press
            io = IOBuffer()

            if (Keyboard.left_shift == events(fig).keyboardbutton[].key || Keyboard.right_shift == events(fig).keyboardbutton[].key) &&
                    events(fig).keyboardbutton[].action == Keyboard.press
                # link identification
                clickpos = AG.createpoint(collect(mouseposition(ax.scene).data)) # in data coordinates

                linkinfos = String[]
                linkdists = []
                
                for link in eachrow(state[].all_links)
                    candidate_dist = AG.distance(link.geometry, clickpos)
                    if candidate_dist < 10
                        push!(linkinfos, """
                            Link @ $(round(candidate_dist, digits=2)):
                                Length:
                                    Geographic: $(link.geographic_length_m)m
                                    Network: $(link.network_length_m)m
                                From edge:
                                    From start: $(link.fr_dist_from_start)m
                                    To end: $(link.fr_dist_to_end)m
                                To edge:
                                    From start: $(link.to_dist_from_start)m
                                    To end: $(link.to_dist_to_end)m
                                """)
                        push!(linkdists, candidate_dist)
                    end
                end

                linkinfo[] = join(linkinfos[sortperm(linkdists)][1:min(2, length(linkinfos))], "\n")
            else
                # distance measurement
                if !isnothing(ruler_p1[]) && isnothing(ruler_p2[])
                    # finishing a measurement
                    ruler_p2[] = mouseposition(ax.scene)
                else
                    # starting a new measurement
                    ruler_p1[] = mouseposition(ax.scene)
                    ruler_p2[] = nothing
                end
            end
        end
    end

    on(events(fig).mouseposition) do _
        if is_mouseinside(ax.scene)
            cursorpos[] = mouseposition(ax.scene)
        end
    end

    sg = SliderGrid(
        iface[1, 1],
        (label="Density (int/km²)", range=10:1:100, startvalue=50)
    )

    iface[2, 1] = toggles = GridLayout()

    cbox = Checkbox(toggles[1, 1], checked=false)
    Label(toggles[1, 2], "Show edge lengths")

    # find middle of line segment
    lengths = @lift(AG.geomlength.($state.fuzzed_graph.edges.geometry))
    labelpts = @lift(AG.pointalongline.($state.fuzzed_graph.edges.geometry, $lengths * 0.5))
    labelx = @lift(AG.getx.($labelpts, 0))
    labely = @lift(AG.gety.($labelpts, 0))
    lenround = @lift($(cbox.checked) ? repr.(round.($lengths, digits=2)) : fill("", length($labelx)))
    text!(ax, labelx, labely, text=lenround)

    distance = @lift(isnothing($ruler_p2) ? 0.0 : round(norm2($ruler_p1 - $ruler_p2), digits=1))
    distreadout = @lift("Distance: $($distance)m")
    Label(iface[3, 1], distreadout)

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
    Label(iface[4, 1], seed)

    seedbox = Textbox(iface[5, 1], placeholder = "Set seed...", width=290, validator = UInt64)

    on(seedbox.stored_string) do s
        linkinfo[] = ""
        state[] = next_state(FuzzedGraphSettings(intersection_density = sg.sliders[1].value[] / (1000^2)), parse(UInt64, s))
    end

    iface[6, 1] = nextbutton = Button(fig, label="Random graph")

    on(nextbutton.clicks) do _
        @debug("New state with seed $(seedval[1])")
        linkinfo[] = ""
        state[] = next_state(FuzzedGraphSettings(intersection_density = sg.sliders[1].value[] / (1000^2)))
    end

    poslabel = @lift("Coordinate: $(round($cursorpos[1], digits=2)), $(round($cursorpos[2], digits=2))")
    Label(iface[7, 1], poslabel)

    Label(iface[8, 1], linkinfo, justification=:left)

    Label(iface[10, 1], """
        Usage
        Click two points to measure
        Shift-click to see link info for links within 10m
        """, justification=:left)

    # because display() uses GLMakie, and we loaded GLMakie at runtime with eval() so that it is only
    # loaded if the interactive visualizer is requested, we need to use invokelatest so that it is
    # available to render the display.
    Base.invokelatest(display, fig)

    while true
        sleep(0.1)
        yield()
    end
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
        isempty(all_links) ? DataFrame(geometry=[]) : links_to_gis(G, all_links,
            :fr_dist_from_start => map(l -> l.fr_dist_from_start, all_links),
            :fr_dist_to_end => map(l -> l.fr_dist_to_end, all_links),
            :to_dist_from_start => map(l -> l.to_dist_from_start, all_links),
            :to_dist_to_end => map(l -> l.to_dist_to_end, all_links),
            :geographic_length_m => map(l -> l.geographic_length_m, all_links)
        ),
        isempty(all_links) ? DataFrame(geometry=[]) : links_to_gis(G, links)
    )
end