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

  infections_immunity_kind = zeros(Int, 6, max_days + 1)
  death_immunity_kind = zeros(Int, 6, max_days + 1)
  hospitalization_immunity_kind = zeros(Int, 6, max_days + 1)
  hospitalization_release_immunity_kind = zeros(Int, 6, max_days + 1)
  for i in 1:num_individuals
    event = MocosSim.backwardinfection(state, i)
    kind = contactkind(event)
    contact_kinds[i] = kind
    infection_times[i] = ifelse(kind == MocosSim.NoContact, missing, time(event))
    severity = state.progressions[i].severity
    non_asymptomatic[i] = ifelse(severity == MocosSim.Asymptomatic, missing, 1.0)
  end
  hospitalization_progressions = getproperty.(state.progressions, :severe_symptoms_time)
  recovery_progressions = getproperty.(state.progressions, :recovery_time)
  death_progressions = getproperty.(state.progressions, :death_time)
  release_progressions = coalesce.(recovery_progressions, death_progressions)
  hospital_release_progressions = (hospitalization_progressions .- hospitalization_progressions) .+ release_progressions
  for i in 1:num_individuals
    if non_asymptomatic[i] !== missing
      if infection_times[i] !== missing && infection_times[i] <= max_days
        immunity_int = state.individuals[i].immunity |> UInt8
        time_int = infection_times[i] + 1 |> floor |> Int
        infections_immunity_kind[immunity_int,time_int] += 1
      end
      if death_progressions[i] !== missing && infection_times[i] + death_progressions[i] <= max_days
        immunity_int = state.individuals[i].immunity |> UInt8
        time_int = infection_times[i] + death_progressions[i] + 1 |> floor |> Int
        death_immunity_kind[immunity_int,time_int] += 1
      end
      if hospitalization_progressions[i] !== missing && infection_times[i] +  hospitalization_progressions[i] <= max_days
        immunity_int = state.individuals[i].immunity |> UInt8
        time_int = infection_times[i] + hospitalization_progressions[i] + 1 |> floor |> Int
        hospitalization_immunity_kind[immunity_int,time_int] += 1
      end
      if hospital_release_progressions[i] !== missing && infection_times[i] +  hospital_release_progressions[i] <= max_days
        immunity_int = state.individuals[i].immunity |> UInt8
        time_int = infection_times[i] + hospital_release_progressions[i] + 1 |> floor |> Int
        hospitalization_release_immunity_kind[immunity_int,time_int] += 1
      end
    end
  end
  dict["daily_infections"] = daily(filter(!ismissing, infection_times .* non_asymptomatic), max_days)
  dict["daily_detections"] = daily(filter(!ismissing, cb.detection_times .* non_asymptomatic), max_days)
  dict["daily_deaths"] = daily(filter(!ismissing, infection_times.+death_progressions), max_days)
  dict["daily_hospitalizations"] = daily(filter(!ismissing, (infection_times.+hospitalization_progressions) .* non_asymptomatic), max_days)
  dict["daily_hospital_releases"] = daily(filter(!ismissing, (infection_times.+hospital_release_progressions) .* non_asymptomatic), max_days)
  for kind in instances(MocosSim.ContactKind)
    if kind != NoContact
      dict["daily_" * lowercase(string(kind))] = daily(infection_times[contact_kinds.==Int(kind)], max_days)
    end
  end
  for immunity in instances(MocosSim.ImmunityState)
    if immunity !== MocosSim.NullImmunity
      immunity_int = immunity |> UInt8
      dict["daily_infections_" * lowercase(string(immunity))] = infections_immunity_kind[immunity_int,:]
      dict["daily_death_" * lowercase(string(immunity))] = death_immunity_kind[immunity_int,:]
      dict["daily_hospitalizations_" * lowercase(string(immunity))] = hospitalization_immunity_kind[immunity_int,:]
      dict["daily_hospital_releases_" * lowercase(string(immunity))] = hospitalization_release_immunity_kind[immunity_int,:]
    end
  end
end