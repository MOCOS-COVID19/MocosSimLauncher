abstract type Output end

beforetrajectories(::Output, ::MocosSim.SimParams) = nothing
pushtrajectory!(::Output, ::Integer, ::Base.AbstractLock, ::MocosSim.SimState, ::MocosSim.SimParams, ::DetectionCallback) = nothing
aftertrajectories(::Output, ::MocosSim.SimParams) = nothing

include("outputs/daily_trajectories.jl")
include("outputs/params_dump.jl")
include("outputs/run_dump.jl")
include("outputs/summary.jl")

const cmd_to_output = Dict{String,Type}(
  "output-summary" => Summary,
  "output-daily" => DailyTrajectories,
  "output-params-dump" => ParamsDump,
  "output-run-dump-prefix" => RunDump
)

function make_outputs(cmd_args::Dict{String,Any}, num_trajectories::Integer)
  outputs = Output[]
  for (opt, output_type) in cmd_to_output
    if cmd_args[opt] |> !isnothing
      push!(outputs, output_type(cmd_args[opt], num_trajectories))
    end
  end
  outputs
end