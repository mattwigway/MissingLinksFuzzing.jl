function ui(raw_args=ARGS)
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--interactive", "-i"
            help = "Start the interactive visualizer"
            action = "store_true"
        "--once", "-o"
            help = "Run once (useful to re-check a problematic seed)"
            action = "store_true"
        "--seed", "-s"
            help = "Seed the algorithm for reproducible runs"
            default = rand(UInt64)
            arg_type=UInt64
        "--log", "-l"
            help = "CSV file to log score differences to, for later statistical analysis. If filename contains %i, will be replaced with initial seed."
            default=nothing
            arg_type = String
        "--headless"
            help = "Format output appropriately for a headless environment, using progress log messages instead of a progress bar"
            action = "store_true"
    end

    args = parse_args(raw_args, s)

    if args["interactive"]
        display(visualizer())
        while true
            sleep(0.1)
            yield()
        end
    else
        fuzzer(seed=args["seed"], logfile=args["log"], once=args["once"], headless=args["headless"])
    end
end