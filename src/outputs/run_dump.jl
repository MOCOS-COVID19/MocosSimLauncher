struct RunDump <: Output
  path_prefix::String
  RunDump(path_prefix::AbstractString) = new(path_prefix)
end
RunDump(path_prefix::AbstractString, ::Integer) = RunDump(path_prefix)


function pushtrajectory!(d::RunDump, trajectory_id::Integer, writelock::Base.AbstractLock, state::MocosSim.SimState, params::MocosSim.SimParams, callback::DetectionCallback)
  try lock(writelock)
    f = jldopen(d.path_prefix*"_$trajectory_id.jld2", "w", compress=true)
    try
      dict = JLD2.Group(f, string(trajectory_id))
      num_individuals = MocosSim.numindividuals(state)
      infection_times = Vector{OptTimePoint}(missing, num_individuals)
      detection_times = Vector{OptTimePoint}(missing, num_individuals)
      death_times = Vector{OptTimePoint}(missing, num_individuals)
      hosp_times = Vector{OptTimePoint}(missing, num_individuals)
      recovery_times = Vector{OptTimePoint}(missing, num_individuals)
      incubation_times = Vector{OptTimePoint}(missing, num_individuals)

      # Prepare an array to store the contact kind as strings
      contact_kinds = Vector{String}(undef, num_individuals)

      source_ids = Vector{UInt32}(undef, num_individuals)

      hospitalization_progressions = getproperty.(state.progressions, :severe_symptoms_time)
      death_progressions = getproperty.(state.progressions, :death_time)
      recovery_progressions = getproperty.(state.progressions, :recovery_time)
      incubation_progressions = getproperty.(state.progressions, :incubation_time)
      for i in 1:num_individuals
        event = MocosSim.backwardinfection(state, i)
        kind = contactkind(event)
        infection_times[i] = ifelse(kind == MocosSim.NoContact, missing, time(event))
        detection_times[i] = callback.detection_times[i]
        death_times[i] = infection_times[i] + death_progressions[i]
        hosp_times[i] = infection_times[i] + hospitalization_progressions[i]
        recovery_times[i] = infection_times[i] + recovery_progressions[i]
        incubation_times[i] = infection_times[i] + incubation_progressions[i]
        # Convert enum to string
        contact_kinds[i] = string(kind)

        source_ids[i] = ifelse(kind == MocosSim.NoContact, undef, source(event))
      end
      dict["detections"] = detection_times
      dict["infections"] = infection_times
      dict["deaths"] = death_times
      dict["hospitalizations"] = hosp_times
      dict["recovery_times"] = recovery_times
      dict["incubation_times"] = incubation_times
      dict["contact_kinds"] = contact_kinds
      dict["source_ids"] = source_ids
    finally
      close(f)
    end
  finally
    unlock(writelock)
  end
  nothing
end
