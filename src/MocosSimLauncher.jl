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
include("scoring.jl")

export launch, evaluate_config_window

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

function checkpoint_path_for(base::AbstractString, trajectory_id::Integer)
  trajectory_id == 1 && return base
  root, ext = splitext(base)
  isempty(ext) && (ext = ".jld2")
  string(root, "_", trajectory_id, ext)
end

function load_resume_state(path::AbstractString)
  isfile(path) || throw(ArgumentError("checkpoint file does not exist: $(path)"))
  lock(_JLD2_LOAD_LOCK) do
    f = jldopen(path, "r")
    try
      state = read(f, "state")
      callback = read(f, "callback")
      trajectory_id = haskey(f, "trajectory_id") ? Int(read(f, "trajectory_id")) : 0
      checkpoint_time = haskey(f, "checkpoint_time") ? MocosSim.TimePoint(read(f, "checkpoint_time")) : MocosSim.TimePoint(0)
      window_start = haskey(f, "window_start") ? MocosSim.TimePoint(read(f, "window_start")) : checkpoint_time
      return state, callback, trajectory_id, checkpoint_time, window_start
    finally
      close(f)
    end
  end
end

const _LAUNCH_SETUP_LOCK = ReentrantLock()

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
  checkpoint_time = haskey(config, "checkpoint_time") ? config["checkpoint_time"] |> MocosSim.TimePoint : missing
  checkpoint_output = get(cmd_args, "output-checkpoint", nothing)
  resume_checkpoint = get(cmd_args, "resume-from-checkpoint", nothing)
  window_start = get(config, "window_start", 0) |> MocosSim.TimePoint
  if resume_checkpoint !== nothing
    num_trajectories == 1 || throw(ArgumentError("resuming from one checkpoint requires num_trajectories = 1"))
    isfile(resume_checkpoint) || throw(ArgumentError("checkpoint file does not exist: $(resume_checkpoint)"))
  end

  # Serialise all file I/O so concurrent launch() calls don't race on JLD2/HDF5.
  # Lock covers the setup phase only; released before the @threads simulation loop.
  local params, num_individuals, immune, immune_events,
        states, callbacks, contexts, outputs,
        writelock, writelock2, progress, outside_case_imports
  lock(_LAUNCH_SETUP_LOCK) do

  @info "loading population and setting up parameters" params_seed
  rng = MersenneTwister(params_seed)
  GC.gc()
  params = read_params(config, rng)
  GC.gc()
  num_individuals =  MocosSim.numindividuals(params)

  immune = nothing
  immune_events = nothing
  if haskey(config, "initial_conditions")
    initial_conditions = config["initial_conditions"]
    if haskey(initial_conditions, "immunization")
      immunization_cfg = initial_conditions["immunization"]

      if haskey(immunization_cfg, "immunity_events")
        immune_events = load(immunization_cfg["immunity_events"], "events")::MocosSim.ImmunizationEvents
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
  contexts = [TrajectoryContext() for _ in 1:nthreads()]
  outputs = make_outputs(cmd_args, num_trajectories)

  for o in outputs
    beforetrajectories(o, params)
  end

  @info "starting simulation" num_trajectories
  writelock = ReentrantLock()
  writelock2 = ReentrantLock()
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

  end  # lock(_LAUNCH_SETUP_LOCK) — setup complete, simulation runs concurrently from here

  trajectory_errors = Tuple{Int,Any,Any}[]
  error_lock = ReentrantLock()
  @threads for trajectory_id in 1:num_trajectories
    state = states[threadid()]
    callback = callbacks[threadid()]
    context = contexts[threadid()]
    reset!(context)
    context.trajectory_id = trajectory_id
    context.window_start = window_start

    if isnothing(resume_checkpoint)
      MocosSim.reset!(state, trajectory_id)
      reset!(callback)
      if checkpoint_output !== nothing
        context.checkpoint_path = checkpoint_path_for(checkpoint_output, trajectory_id)
        context.checkpoint_time = checkpoint_time
      end

      for outside_fun in outside_case_imports
        outside_fun(state, params)
      end
      if params.screening_params !== nothing
        MocosSim.add_screening!(state, params)
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
      if immune_events !== nothing
        MocosSim.immunize!(state, immune_events)
      end
    else
      resumed_state, resumed_callback, _, resumed_time, resumed_window_start = load_resume_state(resume_checkpoint)
      states[threadid()] = resumed_state
      callbacks[threadid()] = resumed_callback
      state = resumed_state
      callback = resumed_callback
      context.window_start = resumed_window_start
      if !ismissing(checkpoint_time) && resumed_time >= checkpoint_time
        context.checkpoint_written = true
      end
      if checkpoint_output !== nothing
        context.checkpoint_path = checkpoint_path_for(checkpoint_output, trajectory_id)
        context.checkpoint_time = checkpoint_time
      end
    end

    try
      local callback_wrapper = (event, simstate, simparams) -> callback(event, simstate, simparams, context)
      MocosSim.simulate!(state, params, callback_wrapper)
      for o in outputs
        pushtrajectory!(o, trajectory_id, writelock, state, params, callback)
      end
      if cmd_args["output-run-dump-prefix"] !== nothing
        path = cmd_args["output-run-dump-prefix"] * "_t$(trajectory_id).jld2"
        save_infections_and_detections(path, writelock2, state, callback)
      end
    catch err
      bt = catch_backtrace()
      lock(error_lock) do
        push!(trajectory_errors, (trajectory_id, err, bt))
      end
      @error "Trajectory failed" thread_id=threadid() trajectory_id exception=(err, bt)
    end
    ProgressMeter.next!(progress) # is thread-safe
  end

  for o in outputs
    aftertrajectories(o, params)
  end
  if !isempty(trajectory_errors)
    failed_ids = sort!(first.(trajectory_errors))
    error("$(length(trajectory_errors)) trajectory/trajectories failed: $(join(failed_ids, ", "))")
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
