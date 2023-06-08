using Dates
using Logging
using Crayons

import Logging: shouldlog, min_enabled_level, catch_exceptions, handle_message

struct SubzeroLogger <: Logging.AbstractLogger
    stream::IO
    min_level::Logging.LogLevel
    message_limits::Dict{Any,Int}
    messages_per_tstep::Int
    show_info_source::Bool

    function SubzeroLogger(
        stream,
        min_level,
        message_limits,
        messages_per_tstep,
        show_info_source
    )
        if !haskey(message_limits, "tstep")
            message_limits["tstep"] = 0
        end
        new(
            stream,
            min_level,
            message_limits,
            messages_per_tstep,
            show_info_source,
        )
    end
end

SubzeroLogger(sim; show_info_source=false, messages_per_tstep = 1) =
    SubzeroLogger(
        open("./log/$(sim.name).log", "w+"),
        Logging.Info,
        Dict{Any,Int}(),
        messages_per_tstep,  # number of messages per timestep (per message)
        show_info_source,
    )

shouldlog(logger::SubzeroLogger, level, _module, group, id) = true

min_enabled_level(logger::SubzeroLogger) = Logging.Info

catch_exceptions(logger::SubzeroLogger) = false

function level_to_string(level)
    level == Logging.Error && return "ERROR"
    level == Logging.Warn  && return "WARN "
    level == Logging.Info  && return "INFO "
    level == Logging.Debug && return "DEBUG"
    return string(level)
end

function handle_message(
    logger::SubzeroLogger,
    level,
    message,
    _module,
    group,
    id,
    filepath,
    line;
    tstep = nothing,
    kwargs...,
)
    # Each message should only write `messages_per_tstep` times per timestep
    if !isnothing(tstep)
        if tstep > logger.message_limits["tstep"]
            for id in keys(logger.message_limits)
                logger.message_limits[id] = logger.messages_per_tstep
            end
            logger.message_limits["tstep"] = tstep
        end
        remaining = get!(
            logger.message_limits,
            id,
            logger.messages_per_tstep,
        )
        remaining > 0 || return nothing
        logger.message_limits[id] = remaining - 1
    end

    buf = IOBuffer()
    iob = IOContext(buf, logger.stream)

    # Log type and message
    level_name = level_to_string(level)
    formatted_message = "$level_name $message"
    # Add simulation timestep to message
    if !isnothing(tstep)
        formatted_message *= " --> timestep $tstep"
    end
    # Add wall clock time to message
    msg_timestamp = Dates.format(Dates.now(), "[yyyy/mm/dd HH:MM:SS.sss]")
    formatted_message *= " $msg_timestamp"
    # Add log location in code to message
    if logger.show_info_source || level != Logging.Info
        file_name   = something(filepath, "nothing")
        line_number = something(line, "nothing")
        formatted_message *= " --> $file_name:$line_number"
    end

    # Write message
    println(iob, formatted_message)
    write(logger.stream, take!(buf))
    return
end