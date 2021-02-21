push!(LOAD_PATH, "../MocosSim")

module MocosSimLauncher

using ArgParse
using Base.Threads
using DataFrames
using FileIO
using JLD2
using JSON
using ProgressMeter
using Random

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
  num_trajectories = json["num_trajectories"] |> Int
  num_initial_infected = json["initial_conditions"]["cardinalities"]["infectious"] |> Int
  params_seed = get(json, "params_seed", 0)

  @info "loading population and setting up parameters" params_seed
  rng = MersenneTwister(params_seed)
  params = read_params(json, rng)
  num_individuals =  MocosSim.numindividuals(params)

  states = [MocosSim.SimState(num_individuals) for _ in 1:nthreads()]
  callbacks = [DetectionCallback(num_individuals, max_num_infected) for _ in 1:nthreads()]

  outputs = make_outputs(cmd_args, num_trajectories)
  foreach(o->beforetrajectories(o, params), outputs)

  @info "starting simulation" num_trajectories
  writelock = ReentrantLock()
  progress = ProgressMeter.Progress(num_trajectories)
  GC.gc()
  @threads for trajectory_id in 1:num_trajectories
    state = states[threadid()]
    MocosSim.reset!(state, trajectory_id)
    MocosSim.initialfeed!(state, num_initial_infected)

    callback = callbacks[threadid()]
    reset!(callback)
    try
      MocosSim.simulate!(state, params, callback)
      for o in outputs
        o->pushtrajectory!(o, trajectory_id, writelock, state, params, callback)
      end
    catch err
      @warn "Failed on thread " threadid() trajectory_id err
      foreach(x -> println(stderr, x), stacktrace(catch_backtrace()))
    end

    ProgressMeter.next!(progress) # is thread-safe
  end
  foreach(o->aftertrajectories(o, params), outputs)
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

end
