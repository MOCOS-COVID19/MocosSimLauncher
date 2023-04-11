using Random

function create_modulation(modulation_dict)
  if isnothing(modulation_dict)
    return nothing
  end
  params_dict = get(modulation_dict, "params", Dict{String,Any}())
  modulation_name = modulation_dict["function"]
  if modulation_name == "IntervalsModulations"
    params_dict["interval_values"] = params_dict["interval_values"] |> Vector{Float64}
    params_dict["interval_times"] = params_dict["interval_times"] |> Vector{MocosSim.TimePoint}
  end
  modulation_params =  NamedTuple{Tuple(Symbol.(keys(params_dict)))}(values(params_dict))

  MocosSim.make_infection_modulation( modulation_name; modulation_params...)
end

function read_params(config, rng::AbstractRNG)
  population_path = config["population_path"] # <= JSON
  population_path::AbstractString # checks if it was indeed a string
  individuals_df = load(population_path)["individuals_df"]

  effectiveness_table = Float64[0.0 0.0 0.0 0.0]

  hospitalization_ratio = get(config["initial_conditions"], "hospitalization_ratio", 1.0) |> Float64
  hospitalization_multiplier = get(config["initial_conditions"], "hospitalization_multiplier", 1.0) |> Float64
  death_multiplier = get(config["initial_conditions"], "death_multiplier", 1.0) |> Float64
  incubation_ratio = get(config["initial_conditions"], "incubation_ratio", [1.0, 1.0, 1.0, 1.0]) |> Vector{Float64}
  progression_params = MocosSim.make_progression_params(hospitalization_ratio, hospitalization_multiplier, death_multiplier, incubation_ratio)

  infection_modulation = get(config, "infection_modulation", nothing) |> create_modulation
  mild_detection_modulation = get(config, "mild_detection_modulation", nothing) |> create_modulation
  tracing_modulation = get(config, "tracing_modulation", nothing) |> create_modulation

  constant_kernel_param = config["transmission_probabilities"]["constant"]  |> float
  household_kernel_param = config["transmission_probabilities"]["household"] |> float
  class_kernel_param = get(config["transmission_probabilities"], "class", 0.0) |> float
  school_kernel_param = get(config["transmission_probabilities"], "school", 0.0) |> float
  hospital_detections = get(config, "hospital_detections", true) |> Bool
  mild_detection_prob = config["mild_detection_prob"]  |> float
  mild_detection_delay = get(config, "mild_detection_delay", 2.0) |> float

  tracing_prob = config["contact_tracing"]["probability"]  |> float
  tracing_delay = config["contact_tracing"]["detection_delay"]  |> float
  testing_time = config["contact_tracing"]["testing_time"]  |> float
  quarantine_length = get(config, "quarantine_length", 14.0) |> float

  age_coupling_kernel_param = get(config["transmission_probabilities"], "age_coupling_param", nothing)
  age_coupling_data_path = get(config["transmission_probabilities"], "age_coupling_data_path", nothing)
  @assert isnothing(age_coupling_kernel_param) == isnothing(age_coupling_data_path)
  age_coupling_thresholds, age_coupling_weights, age_coupling_use_genders =
    isnothing(age_coupling_data_path) ? (nothing, nothing, false) : load(age_coupling_data_path, "age_thresholds", "contact_mat", "uses_genders")

  screening_params = if !haskey(config, "screening")
      nothing
    else
      screen = config["screening"]
      NamedTuple{Tuple(Symbol.(keys(screen)))}(values(screen))
      MocosSim.ScreeningParams(;NamedTuple{Tuple(Symbol.(keys(screen)))}(values(screen))...)
    end

  spreading = get(config, "spreading", nothing)
  spreading_alpha = isnothing(spreading) ? nothing : spreading["alpha"]
  spreading_x0 = isnothing(spreading) ? 1 : get(spreading, "x0", 1)
  spreading_truncation = isnothing(spreading) ? Inf : get(spreading, "truncation", Inf)

  household_params = if !haskey(config, "household_params")
      nothing
    else
      hparams = config["household_params"]
      NamedTuple{Tuple(Symbol.(keys(hparams)))}(values(hparams))
      MocosSim.HouseholdParams(;NamedTuple{Tuple(Symbol.(keys(hparams)))}(values(hparams))...)
    end

  phone_tracing = get(config, "phone_tracing", nothing)
  phone_tracing_usage = isnothing(phone_tracing) ? 0.0 : phone_tracing["usage"] |> float
  phone_tracing_testing_delay = isnothing(phone_tracing) ? 1.0 : phone_tracing["detection_delay"] |> float
  phone_tracing_usage_by_household = isnothing(phone_tracing) ? false : phone_tracing["usage_by_household"] |> Bool

  british_strain_multiplier = get(config["transmission_probabilities"], "british_strain_multiplier", 1.5) |> float
  delta_strain_multiplier = get(config["transmission_probabilities"], "delta_strain_multiplier", 1.5 * 1.5) |> float
  omicron_strain_multiplier = get(config["transmission_probabilities"], "omicron_strain_multiplier", 5.1) |> float
  @info omicron_strain_multiplier
  hospital_kernel_param = get(config["transmission_probabilities"], "hospital", 0.0) |> float
  healthcare_detection_prob, healthcare_detection_delay =  if !haskey(config, "healthare_detections")
    0.8, 1.0
  else
    float(config["healthare_detections"]["probability"]), float(config["healthare_detections"]["delay"])
  end

  MocosSim.make_params(
    rng;
    individuals_df = individuals_df,
    progression_params = progression_params,

    infection_modulation = infection_modulation,
    mild_detection_modulation = mild_detection_modulation,
    forward_tracing_modulation = tracing_modulation,
    backward_tracing_modulation = tracing_modulation,

    constant_kernel_param = constant_kernel_param,
    household_kernel_param = household_kernel_param,
    class_kernel_param = class_kernel_param,
    school_kernel_param = school_kernel_param,

    hospital_detections = hospital_detections,
    mild_detection_prob = mild_detection_prob,
    mild_detection_delay = mild_detection_delay,

    backward_tracing_prob = tracing_prob,
    backward_detection_delay = tracing_delay,

    forward_tracing_prob = tracing_prob,
    forward_detection_delay = tracing_delay,

    testing_time = testing_time,
    quarantine_length = quarantine_length,

    age_coupling_param = age_coupling_kernel_param,
    age_coupling_thresholds = age_coupling_thresholds,
    age_coupling_weights = age_coupling_weights,
    age_coupling_use_genders = age_coupling_use_genders,

    screening_params = screening_params,
    household_params = household_params,

    spreading_alpha=spreading_alpha,
    spreading_x0=spreading_x0,
    spreading_truncation=spreading_truncation,

    phone_tracing_usage = phone_tracing_usage,
    phone_detection_delay = phone_tracing_testing_delay,
    phone_tracing_usage_by_household = phone_tracing_usage_by_household,

    british_strain_multiplier = british_strain_multiplier,
    delta_strain_multiplier = delta_strain_multiplier,
    omicron_strain_multiplier = omicron_strain_multiplier,

    hospital_kernel_param = hospital_kernel_param,
    healthcare_detection_prob = healthcare_detection_prob,
    healthcare_detection_delay = healthcare_detection_delay
  )
end
