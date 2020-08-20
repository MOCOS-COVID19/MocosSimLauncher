struct Summary <: Output
  filename::String
  last_infections::Vector{Float32}
  num_infections::Vector{UInt32}
  Summary(filename::AbstractString, num_trajectories::Integer) = new(filename, fill(NaN32, num_trajectories), Vector{UInt32}(undef, num_trajectories))
end

function pushtrajectory!(s::Summary, trajectory_id::Integer, writelock::Base.AbstractLock, state::MocosSim.SimState, ::MocosSim.SimParams, ::DetectionCallback)
  try 
    lock(writelock)
    s.last_infections[trajectory_id] = maximum(MocosSim.time, state.forest.inedges)
    s.num_infections[trajectory_id] = count(MocosSim.istransmission, state.forest.inedges)
  finally
    unlock(writelock)
  end
  nothing
end

function aftertrajectories(s::Summary, ::MocosSim.SimParams)
  file = jldopen(s.filename, "w")
  try 
    file["num_infections"] = s.num_infections
    file["last_infections"] = s.last_infections
  finally
    close(file)
  end
  nothing
end