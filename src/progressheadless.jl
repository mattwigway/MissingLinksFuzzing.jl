# headless version of ProgressUnknown, for environments without a console (e.g. slurm)
mutable struct ProgressUnknownHeadless
    i::Int64
    mod::Int64
    message::String
end

ProgressUnknownHeadless(mod::Int64, msg::String) = ProgressUnknownHeadless(0, mod, msg)

function ProgressMeter.next!(p::ProgressUnknownHeadless)
    p.i += 1

    if p.i % p.mod == 0
        @info "$(p.message) $(p.i)"
    end
end