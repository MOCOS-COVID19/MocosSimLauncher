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

import MocosSim
import MocosSim: ContactKind, NoContact, contactkind, time

const OptTimePoint = Union{Missing, MocosSim.TimePoint}
optreal2float32(optreal::Union{Missing,T} where T<:Real) = ismissing(optreal) ? NaN32 : Float32(optreal)

include("cmd_parsing.jl")
include("callback.jl")
include("load_params.jl")
include("outputs.jl")

export launch

function launch(args::AbstractVector{T} where T<:AbstractString)
  @info "Stated" nthreads()
  if nthreads() == 1
    @warn "using single thread, set more threads by passing --threads agrument to julia or setting JULIA_NUM_THREADS environment variable"
  end

  cmd_args = parse_commandline(args)
  @info "Parsed args" cmd_args
  json = JSON.parsefile(cmd_args["JSON"])

  max_num_infected = json["stop_simulation_threshold"] |> Int
  time_limit = get(json, "stop_simulation_time", typemax(MocosSim.TimePoint)) |> MocosSim.TimePoint
  num_trajectories = json["num_trajectories"] |> Int
  params_seed = get(json, "params_seed", 0)

  @info "loading population and setting up parameters" params_seed
  rng = MersenneTwister(params_seed)
  GC.gc()
  params = read_params(json, rng)
  GC.gc()
  num_individuals =  MocosSim.numindividuals(params)

  immune = nothing
  if haskey(json, "initial_conditions")
    initial_conditions = json["initial_conditions"]
    if haskey(initial_conditions, "immunization")
      immunization = initial_conditions["immunization"]
      immunization_ordering::AbstractVector{T} where T<: Integer = load(immunization["order_data"], "ordering")
      immunization_level::Real = immunization["level"]

      num_immune = round(UInt, immunization_level * num_individuals)
      immune_ids = @view immunization_ordering[begin : num_immune]
      immune = falses(num_individuals)
      immune[immune_ids] .= true
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

  import_funcs = MocosSim.AbstractOutsideCases[]
  if haskey(json,"imported_cases")
    for fun in json["imported_cases"]
      import_name = fun["function"]
      import_preparams = get(fun, "params", Dict{String,Any}())
      if haskey(import_preparams,"strain")
        strain = import_preparams["strain"] *"Strain"
        import_preparams["strain"] = getfield(MocosSim, Symbol(strain))::MocosSim.StrainKind
      end
      import_params = NamedTuple{Tuple(Symbol.(keys(import_preparams)))}(values(import_preparams))
      import_fun = MocosSim.make_imported_cases(import_name; import_params...)
      push!(import_funcs,import_fun)
    end
  else
    error("the import function was not used!")
  end

  @threads for trajectory_id in 1:num_trajectories
    state = states[threadid()]
    MocosSim.reset!(state, trajectory_id)

    for fun in import_funcs
      fun(state, params)
    end
    if params.screening_params !== nothing
      MocosSim.add_screening!(state, params)
    end

    if immune !== nothing
      immune::AbstractVector{Bool}
      @info "immunizing" count(immune)
      for i in 1:num_individuals
        if !immune[i]
          continue
        end
        individual = state.individuals[i]
        state.individuals[i] = @set individual.health = MocosSim.Recovered
      end
      @info "immunized" count(getproperty.( state.individuals, :health) .== MocosSim.Recovered)
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

