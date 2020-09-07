using StatsBase
using JLD2

daily(times::AbstractVector{T} where T, max_days::Integer) = fit(Histogram, times, 0:(max_days) ).weights
daily(times::AbstractVector{T} where T) = daily(times, maximum(times))

struct DailyTrajectories <: Output
  file::JLD2.JLDFile
end

DailyTrajectories(fname::AbstractString) = DailyTrajectories(JLD2.jldopen(fname, "w+", compress=true))
DailyTrajectories(fname::AbstractString, ::Integer) = DailyTrajectories(fname)

function pushtrajectory!(d::DailyTrajectories, trajectory_id::Integer, writelock::Base.AbstractLock, state::MocosSim.SimState, params::MocosSim.SimParams, cb::DetectionCallback)
  try 
    lock(writelock)
    trajectory_group = JLD2.Group(d.file, string(trajectory_id))
    save_daily_trajectories(trajectory_group, state, params, cb)
  finally
    unlock(writelock)
  end
  nothing
end

aftertrajectories(d::DailyTrajectories) = close(d.file)

function save_daily_trajectories(dict, state::MocosSim.SimState, params::MocosSim.SimParams, cb::DetectionCallback)
  max_days = MocosSim.time(state) |> floor |> Int
  num_individuals = MocosSim.numindividuals(state)
  
  infection_times = Vector{OptTimePoint}(missing, num_individuals)
  contact_kinds = Vector{MocosSim.ContactKind}(undef, num_individuals)

  for i in 1:num_individuals
    event = MocosSim.backwardinfection(state, i)
    kind = contactkind(event)
    contact_kinds[i] = kind
    infection_times[i] = ifelse(kind == MocosSim.NoContact, missing, time(event))
  end

  hospitalization_progressions = getproperty.(params.progressions, :severe_symptoms_time)
  death_progressions = getproperty.(params.progressions, :death_time)

  dict["daily_infections"] = daily(filter(!ismissing, infection_times), max_days)
  dict["daily_detections"] = daily(filter(!ismissing, cb.detection_times), max_days)
  dict["daily_deaths"] = daily(filter(!ismissing, infection_times.+death_progressions), max_days)
  dict["daily_hospitalizations"] = daily(filter(!ismissing, infection_times.+hospitalization_progressions), max_days)
  for kind in instances(ContactKind)
    if kind != NoContact
      dict["daily_" * lowercase(string(kind))] = daily(infection_times[contact_kinds.==Int(kind)], max_days)
    end
  end
end