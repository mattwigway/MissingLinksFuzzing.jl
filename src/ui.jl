function ui(raw_args=ARGS)
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--interactive", "-i"
            help = "Start the interactive visualizer"
            action = "store_true"
    end

    args = parse_args(raw_args, s)

    if args["interactive"]
        visualizer()
    end
end