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

aftertrajectories(d::DailyTrajectories, ::MocosSim.SimParams) = close(d.file)

function save_daily_trajectories(dict, state::MocosSim.SimState, params::MocosSim.SimParams, cb::DetectionCallback)
  max_days = MocosSim.time(state) |> floor |> Int
  num_individuals = MocosSim.numindividuals(state)
  #max_ages = maximum(params.ages)
  # Match the RKI age bands: 00-04, 05-14, 15-34, 35-59, 60-79, 80+.
  thresholds = [0, 5, 15, 35, 60, 80]
  age_labels = ["00_04", "05_14", "15_34", "35_59", "60_79", "80_plus"]
  num_agegroup = length(thresholds)
  infection_times = Vector{OptTimePoint}(missing, num_individuals)
  contact_kinds = Vector{MocosSim.ContactKind}(undef, num_individuals)
  non_asymptomatic = Vector{OptTimePoint}(missing, num_individuals)
  attending_schools = Vector{OptTimePoint}(missing, num_individuals)
  strain_instances = [strain for strain in instances(MocosSim.StrainKind) if strain !== MocosSim.NullStrain]
  num_strains = length(strain_instances)
  strain_index = Dict{MocosSim.StrainKind,Int}(strain => idx for (idx, strain) in enumerate(strain_instances))
  strain_per_individual = zeros(Int, num_individuals)
  # infections_immunity_kind = zeros(Int, 6, max_days + 1)
  infections_ages = zeros(Int, num_agegroup, max_days)
  infections_strain_kind = zeros(Int, num_strains, max_days + 1)
  # detections_immunity_kind = zeros(Int, 6, max_days + 1)
  detections_ages = zeros(Int, num_agegroup, max_days)
  detections_strain_kind = zeros(Int, num_strains, max_days + 1)
  # death_immunity_kind = zeros(Int, 6, max_days + 1)
  # hospitalization_immunity_kind = zeros(Int, 6, max_days + 1)
  # hospitalization_release_immunity_kind = zeros(Int, 6, max_days + 1)
  death_ages = zeros(Int, num_agegroup, max_days)
  hospitalization_admissions_ages = zeros(Int, num_agegroup, max_days)
  hospitalization_releases_ages = zeros(Int, num_agegroup, max_days)
  for i in 1:num_individuals
    event = MocosSim.backwardinfection(state, i)
    kind = contactkind(event)
    contact_kinds[i] = kind
    infection_times[i] = ifelse(kind == MocosSim.NoContact, missing, time(event))
    attending_schools[i] = ifelse(params.attending_schools[i], 1.0, missing)
    severity = state.progressions[i].severity
    non_asymptomatic[i] = ifelse(severity == MocosSim.Asymptomatic, missing, 1.0)
    strain_kind = MocosSim.strainkind(event)
    strain_per_individual[i] = get(strain_index, strain_kind, 0)
  end
  hospitalization_progressions = getproperty.(state.progressions, :severe_symptoms_time)
  recovery_progressions = getproperty.(state.progressions, :recovery_time)
  death_progressions = getproperty.(state.progressions, :death_time)
  release_progressions = coalesce.(recovery_progressions, death_progressions)
  hospital_release_progressions = (hospitalization_progressions .- hospitalization_progressions) .+ release_progressions
  for i in 1:num_individuals
    group_ids = MocosSim.agegroup(thresholds, params.ages[i]) |> Int
    if non_asymptomatic[i] !== missing
      if infection_times[i] !== missing && infection_times[i] < max_days
        time_int = infection_times[i] + 1 |> floor |> Int
        infections_ages[group_ids,time_int] += 1
        strain_idx = strain_per_individual[i]
        if strain_idx > 0
          infections_strain_kind[strain_idx, time_int] += 1
        end
        if cb.detection_times[i] !== missing && cb.detection_times[i] < max_days
          time_int = cb.detection_times[i] + 1 |> floor |> Int
          detections_ages[group_ids,time_int] += 1
          if strain_idx > 0
            detections_strain_kind[strain_idx, time_int] += 1
          end
        end
      end
    end
    if infection_times[i] !== missing && death_progressions[i] !== missing &&
       infection_times[i] + death_progressions[i] < max_days
      time_int = infection_times[i] + death_progressions[i] + 1 |> floor |> Int
      death_ages[group_ids,time_int] += 1
    end
    if non_asymptomatic[i] !== missing &&
       infection_times[i] !== missing &&
       hospitalization_progressions[i] !== missing &&
       infection_times[i] + hospitalization_progressions[i] < max_days
      time_int = infection_times[i] + hospitalization_progressions[i] + 1 |> floor |> Int
      hospitalization_admissions_ages[group_ids,time_int] += 1
    end
    if non_asymptomatic[i] !== missing &&
       infection_times[i] !== missing &&
       hospital_release_progressions[i] !== missing &&
       infection_times[i] + hospital_release_progressions[i] < max_days
      time_int = infection_times[i] + hospital_release_progressions[i] + 1 |> floor |> Int
      hospitalization_releases_ages[group_ids,time_int] += 1
    end
  end
  dict["daily_infections"] = daily(skipmissing(infection_times), max_days)
  dict["daily_symptomatic_infections"] = daily(filter(!ismissing, infection_times .* non_asymptomatic), max_days)
  dict["daily_detections"] = daily(filter(!ismissing, cb.detection_times .* non_asymptomatic), max_days)
  dict["daily_deaths"] = daily(filter(!ismissing, infection_times.+death_progressions), max_days)
  dict["daily_hospitalizations"] = daily(filter(!ismissing, (infection_times.+hospitalization_progressions) .* non_asymptomatic), max_days)
  dict["daily_hospital_releases"] = daily(filter(!ismissing, (infection_times.+hospital_release_progressions) .* non_asymptomatic), max_days)
  dict["daily_student_detections"] = daily(filter(!ismissing, cb.detection_times .* attending_schools), max_days)
  dict["daily_age_total_detections"] = dict["daily_detections"]
  dict["daily_age_total_deaths"] = dict["daily_deaths"]
  for kind in instances(MocosSim.ContactKind)
    if kind != NoContact
      dict["daily_" * lowercase(string(kind))] = daily(infection_times[contact_kinds.==kind], max_days)
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
  for group_id in 1:num_agegroup
    label = age_labels[group_id]
    dict["daily_age_$(label)_detections"] = detections_ages[group_id,:]
    dict["daily_age_$(label)_deaths"] = death_ages[group_id,:]
    dict["daily_age_$(label)_infections"] = infections_ages[group_id,:]
    dict["daily_age_$(label)_hospitalization_admissions"] = hospitalization_admissions_ages[group_id,:]
    dict["daily_age_$(label)_hospitalization_releases"] = hospitalization_releases_ages[group_id,:]

    # Preserve the historical threshold-based keys for existing consumers.
    dict["daily_infections_" * string(thresholds[group_id])] = infections_ages[group_id,:]
    dict["daily_detections_" * string(thresholds[group_id])] = detections_ages[group_id,:]
    dict["daily_death_" * string(thresholds[group_id])] = death_ages[group_id,:]
    dict["daily_hospitalization_admissions_" * string(thresholds[group_id])] = hospitalization_admissions_ages[group_id,:]
    dict["daily_hospitalization_releases_" * string(thresholds[group_id])] = hospitalization_releases_ages[group_id,:]
  end
end
