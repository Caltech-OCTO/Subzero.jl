using Dates
using Logging
using Crayons

import Logging: shouldlog, min_enabled_level, catch_exceptions, handle_message

const RED    = Crayon(foreground=:red)
const YELLOW = Crayon(foreground=:light_yellow)
const CYAN   = Crayon(foreground=:cyan)
const BLUE   = Crayon(foreground=:blue)

const BOLD      = Crayon(bold=true)
const UNDERLINE = Crayon(underline=true)

struct SubzeroLogger <: Logging.AbstractLogger
    stream :: IO
    min_level :: Logging.LogLevel
    message_limits :: Dict{Any,Int}
    message_per_timestep::Int
    timestep:: Int
    show_info_source :: Bool
end

SubzeroLogger(sim, show_info_source=false) =
    SubzeroLogger(
        "./log/$(sim.name).log",
        Logging.Info,
        Dict{Any,Int}(),
        1,  # number of messages per timestep (per message)
        0,  # starting timestep,
        show_info_source,
    )

shouldlog(logger::OceananigansLogger, level, _module, group, id) =
    get(logger.message_limits, id, 1) > 0

min_enabled_level(logger::SubzeroLogger) = Logging.Info

catch_exceptions(logger::SubzeroLogger) = false

function handle_message(
    logger::SubzeroLogger,
    args...;
    tstep = nothing,
    kwargs...,
)
    # I only want it to log max_log times *per timestep*
    if !isnothing(tstep) # gotta add something here
        remaining = get!(logger.message_limits, id, maxlog)
        logger.message_limits[id] = remaining - 1
        remaining > 0 || return nothing
    end

    buf = IOBuffer()
    iob = IOContext(buf, logger.stream)

    level_name = level_to_string(level)
    crayon = level_to_crayon(level)

    file_name   = something(filepath, "nothing")
    line_number = something(line, "nothing")
    msg_timestamp = Dates.format(Dates.now(), "[yyyy/mm/dd HH:MM:SS.sss]")

    formatted_message = "$(crayon(msg_timestamp)) $(BOLD(crayon(level_name))) $message"

    if logger.show_info_source || level != Logging.Info
        formatted_message *= " $(BOLD(crayon("-@->"))) $(UNDERLINE("$file_name:$line_number"))"
    end
    # If timestep exists, should add timestep to message
    println(iob, formatted_message)
    write(logger.stream, take!(buf))

    return
end