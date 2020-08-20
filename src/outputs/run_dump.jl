struct RunDump <: Output
  path_prefix::String
  RunDump(path_prefix::AbstractString) = new(path_prefix)
end
RunDump(path_prefix::AbstractString, ::Integer) = RunDump(path_prefix)


function pushtrajectory!(d::RunDump, trajectory_id::Integer, writelock::Base.AbstractLock, state::MocosSim.SimState, ::MocosSim.SimParams, callback::DetectionCallback)
  try lock(writelock)
    f = jldopen(d.path_prefix*"_$trajectory_id.jld2", "w", compress=true)
    try
      MocosSim.saveparams(f, state)
      saveparams(f, callback)
    finally
      close(f)
    end
  finally
    unlock(writelock)
  end
  nothing
end