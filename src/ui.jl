function ui(raw_args=ARGS)
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--interactive", "-i"
            help = "Start the interactive visualizer"
            action = "store_true"
        "--count", "-c"
            help = "Number of times to fuzz the algorithms"
            default = 100
            arg_type = Int64
        "--seed", "-s"
            help = "Seed the algorithm for reproducible runs"
            default = rand(UInt64)
            arg_type=UInt64
        "--log", "-l"
            help = "CSV file to log score differences to, for later statistical analysis"
            default=nothing
            arg_type = String
    end

    args = parse_args(raw_args, s)

    if args["interactive"]
        display(visualizer())
        while true
            sleep(0.1)
            yield()
        end
    else
        fuzzer(args["count"]; seed=args["seed"], logfile=args["log"])
    end
end