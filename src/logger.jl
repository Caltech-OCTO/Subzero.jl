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
    file_logger::Base.SimpleLogger = SimpleLogger(open(log_file, "w"))
    console_logger::Base.SimpleLogger = SimpleLogger(stdout)
    seen_message::Dict{Any,Bool}
end

SubzeroLogger(log_file::String) =
    SubzeroLogger(
        SimpleLogger(open(log_file, "w")),
        SimpleLogger(stdout),
        Dict{Any,Bool}(),
    )

# what should I make the default? Should everything go to terminal? Not go to file? Make default file?

"""
Based on: https://juliacheat.codes/how-to/get-started-with-logging/
"""
shouldlog(logger::SubzeroLogger, level, _module, group, id) = true

min_enabled_level(logger::SubzeroLogger) = Logging.Info

catch_exceptions(logger::SubzeroLogger) = false

function handle_message(logger::SubzeroLogger, args...; kwargs...)
    # Write to logger message to file
    Logging.handle_message(logger.file_logger, args...; kwargs...)
    flush(logger.file_logger.stream)
    # If message hasn't been seen before, write to terminal as well
    seen_message = get!(logger.seen_message, id, false)
    # Do we need a lock for this if statement??
    if !seen_message  # if message hasn't been seen before
        println(iob, formatted_message)
        logger.seen_message[id] = true
        Logging.handle_message(logger.console_logger, args...; kwargs...)
        flush(logger.console_logger.stream)
    end
    return
end