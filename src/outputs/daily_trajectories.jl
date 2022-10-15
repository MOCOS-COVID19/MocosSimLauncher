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
  #max_ages = maximum(params.ages)
  thresholds = [0, 12, 18, 60]
  num_agegroup = length(thresholds)
  infection_times = Vector{OptTimePoint}(missing, num_individuals)
  contact_kinds = Vector{MocosSim.ContactKind}(undef, num_individuals)
  non_asymptomatic = Vector{OptTimePoint}(missing, num_individuals)
  # infections_immunity_kind = zeros(Int, 6, max_days + 1)
  infections_ages = zeros(Int, num_agegroup, max_days + 1)
  # infections_strain_kind = zeros(Int, MocosSim.NUM_STRAINS, max_days + 1)
  # detections_immunity_kind = zeros(Int, 6, max_days + 1)
  detections_ages = zeros(Int, num_agegroup, max_days + 1)
  # detections_strain_kind = zeros(Int, MocosSim.NUM_STRAINS, max_days + 1)
  # death_immunity_kind = zeros(Int, 6, max_days + 1)
  # hospitalization_immunity_kind = zeros(Int, 6, max_days + 1)
  # hospitalization_release_immunity_kind = zeros(Int, 6, max_days + 1)
  death_ages = zeros(Int, num_agegroup, max_days + 1)
  hospitalization_admissions_ages = zeros(Int, num_agegroup, max_days + 1)
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
      group_ids = MocosSim.agegroup(thresholds,params.ages[i]) |> Int
      if infection_times[i] !== missing && infection_times[i] <= max_days
        # immunity_int = state.individuals[i].immunity |> UInt8
        time_int = infection_times[i] + 1 |> floor |> Int
        # infections_immunity_kind[immunity_int,time_int] += 1
        infections_ages[group_ids,time_int] += 1
        # strain_int = state.individuals[i].strain |> UInt8
        # if strain_int != 0#czasami pojawiają się nullstrain
        #   infections_strain_kind[strain_int,time_int] += 1
        # end
        if cb.detection_times[i] !== missing && cb.detection_times[i] <=max_days
          time_int = cb.detection_times[i] + 1 |> floor |> Int
          # detections_immunity_kind[immunity_int,time_int] += 1
          detections_ages[group_ids,time_int] += 1
          # detections_strain_kind[strain_int,time_int] += 1
        end
      end
      if death_progressions[i] !== missing && infection_times[i] + death_progressions[i] <= max_days
        # immunity_int = state.individuals[i].immunity |> UInt8
        time_int = infection_times[i] + death_progressions[i] + 1 |> floor |> Int
        # death_immunity_kind[immunity_int,time_int] += 1
        death_ages[group_ids,time_int] += 1
      end
      if hospitalization_progressions[i] !== missing && infection_times[i] +  hospitalization_progressions[i] <= max_days
        # immunity_int = state.individuals[i].immunity |> UInt8
        time_int = infection_times[i] + hospitalization_progressions[i] + 1 |> floor |> Int
        # hospitalization_immunity_kind[immunity_int,time_int] += 1
        hospitalization_admissions_ages[group_ids,time_int] += 1
      end
      # if hospital_release_progressions[i] !== missing && infection_times[i] +  hospital_release_progressions[i] <= max_days
      #   immunity_int = state.individuals[i].immunity |> UInt8
      #   time_int = infection_times[i] + hospital_release_progressions[i] + 1 |> floor |> Int
      #   hospitalization_release_immunity_kind[immunity_int,time_int] += 1
      # end
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
  # for immunity in instances(MocosSim.ImmunityState)
  #   if immunity !== MocosSim.NullImmunity
  #     immunity_int = immunity |> UInt8
  #     dict["daily_infections_" * lowercase(string(immunity))] = infections_immunity_kind[immunity_int,:]
  #     dict["daily_detections_" * lowercase(string(immunity))] = detections_immunity_kind[immunity_int,:]
  #     dict["daily_death_" * lowercase(string(immunity))] = death_immunity_kind[immunity_int,:]
  #     dict["daily_hospitalizations_" * lowercase(string(immunity))] = hospitalization_immunity_kind[immunity_int,:]
  #     dict["daily_hospital_releases_" * lowercase(string(immunity))] = hospitalization_release_immunity_kind[immunity_int,:]
  #   end
  # end
  #  for strain in instances(MocosSim.StrainKind)
  #   if strain !== MocosSim.NullStrain
  #     strain_int = strain |> UInt8
  #     dict["daily_infections_" * lowercase(string(strain))] = infections_strain_kind[strain_int,:]
  #     dict["daily_detections_" * lowercase(string(strain))] = detections_strain_kind[strain_int,:]
  #   end
  # end
  for group_ids in 1:num_agegroup
    dict["daily_infections_" * string(thresholds[group_ids])] = infections_ages[group_ids,:]
    dict["daily_detections_" * string(thresholds[group_ids])] = detections_ages[group_ids,:]
    dict["daily_death_" * string(thresholds[group_ids])] = death_ages[group_ids,:]
    dict["daily_hospitalization_admissions_" * string(thresholds[group_ids])] = hospitalization_admissions_ages[group_ids,:]
  end
end