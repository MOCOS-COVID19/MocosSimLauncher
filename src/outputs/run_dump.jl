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
      dict = JLD2.Group(d.path_prefix*"_$trajectory_id.jld2", string(trajectory_id))
      num_individuals = MocosSim.numindividuals(state)
      infection_times = Vector{OptTimePoint}(missing, num_individuals)
      death_times = Vector{OptTimePoint}(missing, num_individuals)
      hosp_times = Vector{OptTimePoint}(missing, num_individuals)
      hospitalization_progressions = getproperty.(state.progressions, :severe_symptoms_time)
      death_progressions = getproperty.(state.progressions, :death_time)
      for i in 1:num_individuals
        event = MocosSim.backwardinfection(state, i)
        kind = contactkind(event)
        infection_times[i] = ifelse(kind == MocosSim.NoContact, missing, time(event))
        death_times[i] = infection_times[i] + death_progressions[i]
        hosp_times[i] = infection_times[i] + hospitalization_progressions[i]
      end
      dict["infections"] = infection_times
      dict["deaths"] = death_times
      dict["hospitalizations"] = hosp_times
    finally
      close(f)
    end
  finally
    unlock(writelock)
  end
  nothing
end
