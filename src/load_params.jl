using Random
#using Distributions

function read_params(json, rng::AbstractRNG)
  constant_kernel_param = json["transmission_probabilities"]["constant"]  |> float
  household_kernel_param = json["transmission_probabilities"]["household"] |> float
  hospital_kernel_param = get(json["transmission_probabilities"], "hospital", 0.0) |> float
  # friendship_kernel_param = get(json["transmission_probabilities"], "friendship", 0.0) |> float
  british_strain_multiplier = get(json["transmission_probabilities"], "british_strain_multiplier", 1.7) |> float
  delta_strain_multiplier = get(json["transmission_probabilities"], "delta_strain_multiplier", 1.7* 1.5) |> float

  age_coupling_data_path = get(json["transmission_probabilities"], "age_coupling_data_path", nothing)
  age_coupling_tresholds, age_coupling_weights, age_coupling_use_genders =
    isnothing(age_coupling_data_path) ? (nothing, nothing, false) : load(age_coupling_data_path, "age_thresholds", "age_coupling_weights", "age_coupling_use_genders")
  mild_detection_prob = json["detection_mild_proba"]  |> float

  tracing_prob = json["contact_tracking"]["probability"]  |> float
  tracing_backward_delay = json["contact_tracking"]["backward_detection_delay"]  |> float
  tracing_forward_delay = json["contact_tracking"]["forward_detection_delay"]  |> float
  testing_time = json["contact_tracking"]["testing_time"]  |> float

  phone_tracing = get(json, "phone_tracking", nothing)
  phone_tracing_usage = isnothing(phone_tracing) ? 0.0 : phone_tracing["usage"] |> float
  phone_tracing_testing_delay = isnothing(phone_tracing) ? 1.0 : phone_tracing["detection_delay"] |> float
  phone_tracing_usage_by_household = isnothing(phone_tracing) ? false : phone_tracing["usage_by_household"] |> Bool

  population_path = json["population_path"] # <= JSON
  population_path::AbstractString # checks if it was indeed a string

  individuals_df = load(population_path)["individuals_df"]

  infection_modulation_name, infection_modulation_params = if !haskey(json, "modulation")
    nothing, NamedTuple{}()
  else
    modulation = json["modulation"]
    params = get(modulation, "params", Dict{String,Any}())
    modulation["function"], NamedTuple{Tuple(Symbol.(keys(params)))}(values(params))
  end
  travels_frequency::MocosSim.TimePoint = 0.0
  infection_travels_name, infection_travels_params = if !haskey(json, "travels")
    nothing, NamedTuple{}()
  else
    travels = json["travels"]
    travels_frequency = get(travels, "frequency", 0.05) |> MocosSim.TimePoint
    params2 = get(travels, "params", Dict{String,Any}())
    travels["function"], NamedTuple{Tuple(Symbol.(keys(params2)))}(values(params2))
  end

  screening_params = if !haskey(json, "screening")
    nothing
  else
    screen = json["screening"]
    NamedTuple{Tuple(Symbol.(keys(screen)))}(values(screen))
    MocosSim.ScreeningParams(;NamedTuple{Tuple(Symbol.(keys(screen)))}(values(screen))...)
  end

  spreading = get(json, "spreading", nothing)
  spreading_alpha = isnothing(spreading) ? nothing : spreading["alpha"]
  spreading_x0 = isnothing(spreading) ? 1 : get(spreading, "x0", 1)
  spreading_truncation = isnothing(spreading) ? Inf : get(spreading, "truncation", Inf)

  MocosSim.load_params(
    rng,
    population = individuals_df,

    mild_detection_prob = mild_detection_prob,

    constant_kernel_param = constant_kernel_param,
    household_kernel_param = household_kernel_param,
    hospital_kernel_param = hospital_kernel_param,
    age_coupling_tresholds=age_coupling_tresholds,
    age_coupling_weights=age_coupling_weights,
    age_coupling_use_genders=age_coupling_use_genders,

    backward_tracing_prob = tracing_prob,
    backward_detection_delay = tracing_backward_delay,

    forward_tracing_prob = tracing_prob,
    forward_detection_delay = tracing_forward_delay,

    testing_time = testing_time,

    phone_tracing_usage = phone_tracing_usage,
    phone_detection_delay = phone_tracing_testing_delay,
    phone_tracing_usage_by_household = phone_tracing_usage_by_household,

    infection_modulation_name=infection_modulation_name,
    infection_modulation_params=infection_modulation_params,

    screening_params = screening_params,

    spreading_alpha=spreading_alpha,
    spreading_x0=spreading_x0,
    spreading_truncation=spreading_truncation,
    british_strain_multiplier=british_strain_multiplier,
    delta_strain_multiplier=delta_strain_multiplier
  )
end
