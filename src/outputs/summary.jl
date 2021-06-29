struct Summary <: Output
  filename::String
  last_infections::Vector{Float32}
  num_infections::Vector{UInt32}
  peak_daily_infections::Vector{UInt32}
  peak_daily_detections::Vector{UInt32}
  Summary(filename::AbstractString, num_trajectories::Integer) = new(
    filename,
    fill(NaN32, num_trajectories),
    Vector{UInt32}(undef, num_trajectories),
    Vector{UInt32}(undef, num_trajectories),
    Vector{UInt32}(undef, num_trajectories))
end

time_or_nan(e) = MocosSim.istransmission(e) ? Float32(MocosSim.time(e)) : NaN32

function pushtrajectory!(s::Summary, trajectory_id::Integer, ::Base.AbstractLock, state::MocosSim.SimState, ::MocosSim.SimParams, cb::DetectionCallback)
  max_time = maximum(MocosSim.time, state.forest.inedges)
  s.last_infections[trajectory_id] = max_time
  s.num_infections[trajectory_id] = count(MocosSim.istransmission, state.forest.inedges)
  counts = Vector{UInt32}(undef, ceil(Int, max_time))
  s.peak_daily_infections[trajectory_id] = maximum(daily!(time_or_nan, counts, state.forest.inedges), init=0)
  s.peak_daily_detections[trajectory_id] = maximum(daily!(counts, skipmissing(cb.detection_times)), init=0)
  nothing
end

function aftertrajectories(s::Summary, ::MocosSim.SimParams)
  file = jldopen(s.filename, "w")
  try
    file["num_infections"] = s.num_infections
    file["last_infections"] = s.last_infections
    file["peak_daily_detections"] = s.peak_daily_detections
    file["peak_daily_infections"] = s.peak_daily_infections
  finally
    close(file)
  end
  nothing
end