using StatsBase
using JLD2

function daily!(f, counts::AbstractVector{T} where T<:Integer, arr)
  fill!(counts, 0)
  max_day = length(counts)
  for a in arr
    t = f(a)
    if ismissing(t) || !isfinite(t)
      continue
    end
    day = floor(Int, t) + 1
    if day <= max_day
      counts[day] += 1
    end
  end
  counts
end

daily!(counts, arr) = daily!(identity, counts, arr)
daily(f::Function, times) = daily!(f, fill(0, ceil(Int, maximum(f, times))), times)
daily(times, max_time=ceil(Int, maximum(times))) = daily!(fill(0, max_time), times)

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
  non_asymptomatic = Vector{OptTimePoint}(missing, num_individuals)

  for i in 1:num_individuals
    event = MocosSim.backwardinfection(state, i)
    kind = contactkind(event)
    contact_kinds[i] = kind
    infection_times[i] = ifelse(kind == MocosSim.NoContact, missing, time(event))
    severity = MocosSim.severityof(params, i)
    non_asymptomatic[i] = ifelse(severity == MocosSim.Asymptomatic, missing, 1.0)
  end

  hospitalization_progressions = getproperty.(params.progressions, :severe_symptoms_time)
  recovery_progressions = getproperty.(params.progressions, :recovery_time)
  death_progressions = getproperty.(params.progressions, :death_time)
  release_progressions = coalesce.(recovery_progressions, death_progressions)
  hospital_release_progressions = (!ismissing).(hospitalization_progressions) .* release_progressions

  dict["daily_infections"] = daily(filter(!ismissing, infection_times), max_days)
  dict["daily_detections"] = daily(filter(!ismissing, cb.detection_times), max_days)
  dict["daily_infections_non_asymptomatic"] = daily(filter(!ismissing, infection_times .* non_asymptomatic), max_days)
  dict["daily_detections_non_asymptomatic"] = daily(filter(!ismissing, cb.detection_times .* non_asymptomatic), max_days)
  dict["daily_deaths"] = daily(filter(!ismissing, infection_times.+death_progressions), max_days)
  dict["daily_hospitalizations"] = daily(filter(!ismissing, infection_times.+hospitalization_progressions), max_days)
  dict["daily_hospital_releases"] = daily(filter(!ismissing, infection_times.+hospital_release_progressions), max_days)
  dict["daily_hospitalizations_non_asymptomatic"] = daily(filter(!ismissing, (infection_times .+ hospitalization_progressions) .* non_asymptomatic), max_days)
  dict["daily_hospital_releases_non_asymptomatic"] = daily(filter(!ismissing, (infection_times .+ hospital_release_progressions) .* non_asymptomatic), max_days)
  for kind in instances(ContactKind)
    if kind != NoContact
      dict["daily_" * lowercase(string(kind))] = daily(infection_times[contact_kinds.==Int(kind)], max_days)
    end
  end
end