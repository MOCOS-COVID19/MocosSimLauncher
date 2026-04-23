
struct DetectionCallback
    detection_times::Vector{OptTimePoint}
    detection_types::Vector{UInt8}

    tracing_times::Vector{OptTimePoint}
    tracing_sources::Vector{UInt32}
    tracing_types::Vector{UInt8}

    transmission_times::Vector{OptTimePoint}
    transmission_sources::Vector{UInt32}
    transmission_types::Vector{UInt8}

    max_num_infected::UInt32
    time_limit::MocosSim.TimePoint
end

DetectionCallback(sz::Integer, max_num_infected::Integer=10^8, time_limit::MocosSim.TimePoint=typemax(MocosSim.TimePoint)) = DetectionCallback(
    Vector{OptTimePoint}(missing, sz),
    fill(UInt8(0), sz),
    Vector{OptTimePoint}(missing, sz),
    fill(UInt32(0), sz),
    fill(UInt8(0), sz),
    Vector{OptTimePoint}(missing, sz),
    fill(UInt32(0), sz),
    fill(UInt8(0), sz),
    max_num_infected,
    time_limit
)

mutable struct TrajectoryContext
  trajectory_id::Int
  checkpoint_path::Union{Nothing,String}
  checkpoint_time::OptTimePoint
  checkpoint_written::Bool
  window_start::MocosSim.TimePoint
end

TrajectoryContext() = TrajectoryContext(0, nothing, missing, false, MocosSim.TimePoint(0))

function reset!(context::TrajectoryContext)
  context.trajectory_id = 0
  context.checkpoint_path = nothing
  context.checkpoint_time = missing
  context.checkpoint_written = false
  context.window_start = MocosSim.TimePoint(0)
end

function reset!(cb::DetectionCallback)
  fill!(cb.detection_times, missing)
  fill!(cb.detection_types, 0)
  fill!(cb.tracing_sources, 0)
  fill!(cb.tracing_types, 0)
  fill!(cb.transmission_times, missing)
  fill!(cb.transmission_sources, 0)
  fill!(cb.transmission_types, 0)
end

function save_checkpoint(path::AbstractString, state::MocosSim.SimState, callback::DetectionCallback, context::TrajectoryContext)
  jldopen(path, "w"; compress=true) do f
    f["state"] = state
    f["callback"] = callback
    f["trajectory_id"] = context.trajectory_id
    f["checkpoint_time"] = Float64(MocosSim.time(state))
    f["window_start"] = Float64(context.window_start)
  end
  nothing
end

function maybe_write_checkpoint!(cb::DetectionCallback, state::MocosSim.SimState, context::TrajectoryContext)
  if context.checkpoint_written || isnothing(context.checkpoint_path) || ismissing(context.checkpoint_time)
    return
  end
  if MocosSim.time(state) >= context.checkpoint_time
    save_checkpoint(context.checkpoint_path, state, cb, context)
    context.checkpoint_written = true
  end
  nothing
end

function (cb::DetectionCallback)(event::MocosSim.Event, state::MocosSim.SimState, params::MocosSim.SimParams, context::TrajectoryContext=TrajectoryContext())
  eventkind = MocosSim.kind(event)
  contactkind = MocosSim.contactkind(event)
  subject = MocosSim.subject(event)
  if MocosSim.isdetection(eventkind)
    cb.detection_times[subject] = MocosSim.time(event)
    cb.detection_types[subject] = MocosSim.detectionkind(event) |> UInt8
  elseif MocosSim.istracing(eventkind)
    cb.tracing_times[subject] = MocosSim.time(event)
    cb.tracing_sources[subject] = MocosSim.source(event)
    cb.tracing_types[subject] = MocosSim.tracingkind(event) |> UInt8
  elseif MocosSim.istransmission(eventkind)
    cb.transmission_times[subject] = MocosSim.time(event)
    cb.transmission_sources[subject] = MocosSim.source(event)
    cb.transmission_types[subject] = MocosSim.contactkind(event) |> UInt8
  end
  maybe_write_checkpoint!(cb, state, context)
  return MocosSim.numinfected(state.stats) < cb.max_num_infected && MocosSim.time(event) < cb.time_limit
end

function saveparams(dict, cb::DetectionCallback, prefix::AbstractString="")
  dict[prefix*"detection_times"] = optreal2float32.(cb.detection_times)
  dict[prefix*"detection_types"] = cb.detection_types

  dict[prefix*"tracing_times"] = optreal2float32.(cb.tracing_times)
  dict[prefix*"tracing_sources"] = cb.tracing_sources
  dict[prefix*"tracing_types"] = cb.tracing_types

  dict[prefix*"transmission_times"] = optreal2float32.(cb.transmission_times)
  dict[prefix*"transmission_sources"] = cb.transmission_sources
  dict[prefix*"transmission_types"] = cb.transmission_types
end

function save_infections_and_detections(path::AbstractString, writelock::Base.AbstractLock, simstate::MocosSim.SimState, callback::DetectionCallback)
  try lock(writelock)
    f = jldopen(path, "w"; iotype=IOStream, compress=true)
    try
      # MocosSim.saveparams(f, simstate)
      saveparams(f, callback)
    finally
      close(f)
    end
  finally
    unlock(writelock)
  end
  nothing
end