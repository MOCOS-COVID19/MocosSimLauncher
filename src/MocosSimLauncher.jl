#push!(LOAD_PATH, "../MocosSim")

module MocosSimLauncher

using ArgParse
using Base.Threads
using CodecZlib
using DataFrames
using FileIO
using JLD2
using JSON
using ProgressMeter
using Random
using Setfield
using TOML

import MocosSim
import MocosSim: ContactKind, NoContact, contactkind, time

const OptTimePoint = Union{Missing, MocosSim.TimePoint}
optreal2float32(optreal::Union{Missing,T} where T<:Real) = ismissing(optreal) ? NaN32 : Float32(optreal)

include("cmd_parsing.jl")
include("callback.jl")
include("load_params.jl")
include("outputs.jl")

export launch

function string2enum(enum_group::Type{<:Enum{T}}, str::AbstractString) where {T<:Integer}
    sym = Symbol(str)
    for val in instances(enum_group)
        if sym == Symbol(val)
            return val
        end
    end
    error("not found $str in $enum_group")
end

dict2kwargs(dict::Dict{S, Any} where S<:AbstractString) = NamedTuple{Tuple(Symbol.(keys(dict)))}(values(dict))

function launch(args::AbstractVector{T} where T<:AbstractString)
  @info "Stated" nthreads()
  if nthreads() == 1
    @warn "using single thread, set more threads by passing --threads agrument to julia or setting JULIA_NUM_THREADS environment variable"
  end

  cmd_args = parse_commandline(args)
  @info "Parsed args" cmd_args
  config_path = cmd_args["CONFIG"]
  config = endswith(config_path, ".json") ? JSON.parsefile(config_path) : TOML.parsefile(config_path)

  max_num_infected = config["stop_simulation_threshold"] |> Int
  time_limit = get(config, "stop_simulation_time", typemax(MocosSim.TimePoint)) |> MocosSim.TimePoint
  num_trajectories = config["num_trajectories"] |> Int
  params_seed = get(config, "params_seed", 0)

  @info "loading population and setting up parameters" params_seed
  rng = MersenneTwister(params_seed)
  GC.gc()
  params = read_params(config, rng)
  GC.gc()
  num_individuals =  MocosSim.numindividuals(params)

  immunization = nothing
  immune = nothing
  immune_ages = nothing
  if haskey(config, "initial_conditions")
    initial_conditions = config["initial_conditions"]
    if haskey(initial_conditions, "immunization")
      immunization_cfg = initial_conditions["immunization"]

      if haskey(immunization_cfg, "age_groups")
        immunization_thresholds = get(immunization_cfg["age_groups"], "immunization_thresholds", [0, 12, 18, 60]) |> Vector{Int32}
        immunization_tables = get(immunization_cfg["age_groups"], "immunization_tables", [0.0, 0.0, 39.0, 3.85, 64.0, 28.8, 80.0, 52.13]) |> Vector{Float32}
        immunization_tables = reshape(immunization_tables, 2,length(immunization_thresholds))' |> Matrix{Float32}
        immunization_previously_infected = get(immunization_cfg["age_groups"], "immunization_previously_infected", [0.24, 0.24, 0.24, 0.24]) |> Vector{Float32}
        immune_ages = [immunization_thresholds, immunization_tables, immunization_previously_infected]
      end

      if haskey(immunization_cfg, "order_file")
        immunization = load(immunization_cfg["order_file"], "immunization")::MocosSim.Immunization
        enqueue_immunizations = get(immunization_cfg, "enqueue", true) |> Bool
      end

      # keeping legacy immunization for a while
      if haskey(immunization_cfg, "level")
        @warn "legacy immunization in use"
        immunization_ordering::AbstractVector{T} where T<: Integer = load(immunization_cfg["order_data"], "ordering")
        immunization_level::Real = immunization_cfg["level"]

        num_immune = round(UInt, immunization_level * num_individuals)
        immune_ids = @view immunization_ordering[begin : num_immune]
        immune = falses(num_individuals)
        immune[immune_ids] .= true
      end
    end
  end
  @info "allocating simulation states"
  states = [MocosSim.SimState(num_individuals) for _ in 1:nthreads()]
  callbacks = [DetectionCallback(num_individuals, max_num_infected, time_limit) for _ in 1:nthreads()]
  outputs = make_outputs(cmd_args, num_trajectories)

  for o in outputs
    beforetrajectories(o, params)
  end

  @info "starting simulation" num_trajectories
  writelock = ReentrantLock()
  progress = ProgressMeter.Progress(num_trajectories)
  GC.gc()

  outside_case_imports = MocosSim.AbstractOutsideCases[]
  if haskey(config, "imported_cases")
    for outside_import in config["imported_cases"]
      name = outside_import["function"]
      params_dict = get(outside_import, "params", Dict{String,Any}())
      if haskey(params_dict, "strain")
        strain_str = params_dict["strain"] * "Strain"
        params_dict["strain"] = string2enum(MocosSim.StrainKind, strain_str)
      end
      outside_cases = MocosSim.make_imported_cases(name; dict2kwargs(params_dict)...)
      push!(outside_case_imports, outside_cases)
    end
  else
    error("the import function was not used!")
  end

  @threads for trajectory_id in 1:num_trajectories
    state = states[threadid()]
    MocosSim.reset!(state, trajectory_id)

    for outside_fun in outside_case_imports
      outside_fun(state, params)
    end
    if params.screening_params !== nothing
      MocosSim.add_screening!(state, params)
    end

    if immune_ages !== nothing
      MocosSim.immunize!(state, params, immune_ages[1], immune_ages[2], immune_ages[3])
    end

    if immune !== nothing
      immune::AbstractVector{Bool}
      for i in 1:num_individuals
        if !immune[i]
          continue
        end
        individual = state.individuals[i]
        state.individuals[i] = @set individual.health = MocosSim.Recovered
      end
    end

    if immunization !== nothing
      immunization::MocosSim::Immunization
      MocosSim.immunize!(state, immunization, enqueue_immunizations)
    end

    callback = callbacks[threadid()]
    reset!(callback)
    try
      MocosSim.simulate!(state, params, callback)
      for o in outputs
        pushtrajectory!(o, trajectory_id, writelock, state, params, callback)
      end
    catch err
      @warn "Failed on thread " threadid() trajectory_id err
      foreach(x -> println(stderr, x), stacktrace(catch_backtrace()))
    end

    ProgressMeter.next!(progress) # is thread-safe
  end

  for o in outputs
    aftertrajectories(o, params)
  end
end

function julia_main()::Cint
  try
    launch(ARGS)
  catch err
    Base.invokelatest(Base.display_error, Base.catch_stack())
    return 1
  end
  return 0
end

precompile(MocosSim.simulate!, (MocosSim.SimState, MocosSim.SimParams, DetectionCallback))
precompile(launch, (Vector{String},))
precompile(julia_main, ())

end

