using DataFrames
using JLD2
using MocosSimLauncher
using Random
using Test
using TOML

import MocosSim

const Launcher = MocosSimLauncher

function write_test_population(path::AbstractString)
  individuals = DataFrame(
    household_index=[1, 1, 2],
    attending_school=[1, 0, 0],
    school_index=[1, 0, 0],
    class_index=[1, 0, 0],
    age=[20, 40, 60],
    gender=[0, 1, 0],
    ishealthcare=[false, false, false],
  )
  jldopen(path, "w") do file
    file["individuals_df"] = individuals
  end
  path
end

function minimal_config(population_path::AbstractString)
  Dict{String,Any}(
    "population_path" => population_path,
    "num_trajectories" => 1,
    "mild_detection_prob" => 0.0,
    "stop_simulation_threshold" => 100,
    "stop_simulation_time" => 2,
    "initial_conditions" => Dict{String,Any}(),
    "screening" => Dict{String,Any}(
      "adherence_pmf" => [1.0],
    ),
    "contact_tracing" => Dict{String,Any}(
      "probability" => 0.0,
      "detection_delay" => 1.0,
      "testing_time" => 0.25,
    ),
    "transmission_probabilities" => Dict{String,Any}(
      "constant" => 0.0,
      "household" => 0.0,
    ),
    "imported_cases" => [Dict{String,Any}(
      "function" => "InstantOutsideCases",
      "params" => Dict{String,Any}(
        "import_time" => 0.0,
        "num_infections" => 1,
        "strain" => "Chinese",
      ),
    )],
  )
end

@testset "callback reset" begin
  callback = Launcher.DetectionCallback(2)
  callback.detection_times[1] = 1.0
  callback.detection_types[1] = 2
  callback.tracing_times[1] = 3.0
  callback.tracing_sources[1] = 2
  callback.tracing_types[1] = 1
  callback.transmission_times[1] = 4.0
  callback.transmission_sources[1] = 3
  callback.transmission_types[1] = 2

  Launcher.reset!(callback)

  @test all(ismissing, callback.detection_times)
  @test all(iszero, callback.detection_types)
  @test all(ismissing, callback.tracing_times)
  @test all(iszero, callback.tracing_sources)
  @test all(iszero, callback.tracing_types)
  @test all(ismissing, callback.transmission_times)
  @test all(iszero, callback.transmission_sources)
  @test all(iszero, callback.transmission_types)
end

@testset "trajectory context reset" begin
  context = Launcher.TrajectoryContext(7, "checkpoint.jld2", 12.0, true, 4.0)

  Launcher.reset!(context)

  @test context.trajectory_id == 0
  @test context.checkpoint_path === nothing
  @test ismissing(context.checkpoint_time)
  @test !context.checkpoint_written
  @test context.window_start == 0
end

@testset "summary initialization" begin
  summary = Launcher.Summary("unused.jld2", 3)
  @test summary.num_infections == zeros(UInt32, 3)
  @test summary.peak_daily_infections == zeros(UInt32, 3)
  @test summary.peak_daily_detections == zeros(UInt32, 3)
end

@testset "checkpoint paths" begin
  @test Launcher.checkpoint_path_for("checkpoint.jld2", 1) == "checkpoint.jld2"
  @test Launcher.checkpoint_path_for("checkpoint.jld2", 2) == "checkpoint_2.jld2"
  @test Launcher.checkpoint_path_for("checkpoint", 2) == "checkpoint_2.jld2"
end

@testset "checkpoint round trip" begin
  mktempdir() do directory
    path = joinpath(directory, "checkpoint.jld2")
    state = MocosSim.SimState(2)
    callback = Launcher.DetectionCallback(2)
    callback.detection_times[2] = 3.5
    context = Launcher.TrajectoryContext(4, path, 0.0, false, 2.0)

    Launcher.save_checkpoint(path, state, callback, context)
    loaded_state, loaded_callback, trajectory_id, checkpoint_time, window_start =
      Launcher.load_resume_state(path)

    @test MocosSim.time(loaded_state) == MocosSim.time(state)
    @test isequal(loaded_callback.detection_times, callback.detection_times)
    @test trajectory_id == 4
    @test checkpoint_time == MocosSim.time(state)
    @test window_start == 2
  end

  @test_throws ArgumentError Launcher.load_resume_state("not-a-checkpoint.jld2")
end

@testset "conditional checkpoint writing" begin
  mktempdir() do directory
    path = joinpath(directory, "checkpoint.jld2")
    state = MocosSim.SimState(1)
    callback = Launcher.DetectionCallback(1)
    context = Launcher.TrajectoryContext(1, path, 1.0, false, 0.0)

    Launcher.maybe_write_checkpoint!(callback, state, context)
    @test !isfile(path)
    @test !context.checkpoint_written

    state.time = MocosSim.TimePoint(1)
    Launcher.maybe_write_checkpoint!(callback, state, context)
    @test isfile(path)
    @test context.checkpoint_written

    callback.detection_times[1] = 2.0
    Launcher.maybe_write_checkpoint!(callback, state, context)
    _, saved_callback, _, _, _ = Launcher.load_resume_state(path)
    @test ismissing(saved_callback.detection_times[1])
  end

  state = MocosSim.SimState(1)
  callback = Launcher.DetectionCallback(1)
  Launcher.maybe_write_checkpoint!(callback, state, Launcher.TrajectoryContext())
end

@testset "legacy checkpoint defaults" begin
  mktempdir() do directory
    path = joinpath(directory, "legacy.jld2")
    state = MocosSim.SimState(1)
    callback = Launcher.DetectionCallback(1)
    jldopen(path, "w") do file
      file["state"] = state
      file["callback"] = callback
    end

    _, _, trajectory_id, checkpoint_time, window_start = Launcher.load_resume_state(path)
    @test trajectory_id == 0
    @test checkpoint_time == 0
    @test window_start == checkpoint_time
  end
end

@testset "callback parameter serialization" begin
  callback = Launcher.DetectionCallback(2)
  callback.detection_times[1] = 1.5
  callback.detection_types[1] = 2
  callback.tracing_sources[2] = 8
  callback.transmission_times[2] = 4.0
  output = Dict{String,Any}()

  Launcher.saveparams(output, callback, "callback_")

  @test output["callback_detection_times"][1] == Float32(1.5)
  @test isnan(output["callback_detection_times"][2])
  @test output["callback_detection_types"] == UInt8[2, 0]
  @test output["callback_tracing_sources"] == UInt32[0, 8]
  @test isequal(output["callback_transmission_times"], Float32[NaN, 4.0])
end

@testset "callback records supported event kinds" begin
  mktempdir() do directory
    population_path = write_test_population(joinpath(directory, "population.jld2"))
    params = Launcher.read_params(minimal_config(population_path), MersenneTwister(1))
    state = MocosSim.SimState(MocosSim.numindividuals(params))
    callback = Launcher.DetectionCallback(3)

    detection = MocosSim.Event(
      Val(MocosSim.DetectionEvent), 1.0, 1, MocosSim.FromTracingDetection)
    tracing = MocosSim.Event(
      Val(MocosSim.TracedEvent), 1.5, 2, 1, MocosSim.PhoneTraced)
    transmission = MocosSim.Event(
      Val(MocosSim.TransmissionEvent), 2.0, 3, 2,
      MocosSim.HouseholdContact, MocosSim.ChineseStrain)

    @test callback(detection, state, params)
    @test callback(tracing, state, params)
    @test callback(transmission, state, params)
    @test callback.detection_times[1] == 1.0
    @test callback.detection_types[1] == UInt8(MocosSim.FromTracingDetection)
    @test callback.tracing_times[2] == 1.5
    @test callback.tracing_sources[2] == 1
    @test callback.tracing_types[2] == UInt8(MocosSim.PhoneTraced)
    @test callback.transmission_times[3] == 2.0
    @test callback.transmission_sources[3] == 2
    @test callback.transmission_types[3] == UInt8(MocosSim.HouseholdContact)
  end
end

@testset "minimal parameter loading" begin
  mktempdir() do directory
    population_path = write_test_population(joinpath(directory, "population.jld2"))
    params = Launcher.read_params(minimal_config(population_path), MersenneTwister(7))

    @test MocosSim.numindividuals(params) == 3
    @test params.ages == [20, 40, 60]
    @test params.constant_kernel_param == 0.0
    @test params.household_kernel_param == 0.0
    @test params.screening_params !== nothing
    @test params.phone_tracing_params === nothing
  end
end

@testset "seasonality interval validation" begin
  modulation = Dict{String,Any}(
    "function" => "IntervalsModulations",
    "params" => Dict{String,Any}(
      "interval_times" => [7, 14],
      "interval_values" => [1.0, 1.0],
    ),
  )
  seasonality = Dict{String,Any}(
    "factor" => 0.5,
    "start_day_offset" => 0,
    "end_day_offset" => 7,
  )
  result = Launcher.apply_seasonality_to_intervals!(deepcopy(modulation), seasonality)
  @test result["params"]["interval_values"] == [0.5, 1.0]

  open_ended = deepcopy(modulation)
  push!(open_ended["params"]["interval_values"], 1.0)
  @test_throws ArgumentError Launcher.apply_seasonality_to_intervals!(open_ended, seasonality)

  unsorted = deepcopy(modulation)
  unsorted["params"]["interval_times"] = [14, 7]
  @test_throws ArgumentError Launcher.apply_seasonality_to_intervals!(unsorted, seasonality)

  @test Launcher.apply_seasonality_to_intervals!(deepcopy(modulation), nothing) == modulation
  ignored = merge(seasonality, Dict("target" => "tracing_modulation"))
  @test Launcher.apply_seasonality_to_intervals!(deepcopy(modulation), ignored) == modulation
end

@testset "seasonality handles partial intervals and years" begin
  modulation = Dict{String,Any}(
    "function" => "IntervalsModulations",
    "params" => Dict{String,Any}(
      "interval_times" => [10, 20, 370, 380],
      "interval_values" => [1.0, 1.0, 1.0, 1.0],
    ),
  )
  seasonality = Dict{String,Any}(
    "factor" => 0.5,
    "start_day_offset" => 5,
    "end_day_offset" => 15,
    "years" => [0, 1],
  )

  result = Launcher.apply_seasonality_to_intervals!(modulation, seasonality)

  @test result["params"]["interval_values"] == [0.75, 0.75, 1.0, 0.5]
end

@testset "metric aggregation preserves the longest trajectory" begin
  mktempdir() do directory
    path = joinpath(directory, "daily.jld2")
    jldopen(path, "w") do file
      short = JLD2.Group(file, "1")
      short["daily_detections"] = [2, 2]
      long = JLD2.Group(file, "2")
      long["daily_detections"] = [4, 4, 4, 4]
    end

    @test Launcher.aggregate_metric(path, "daily_detections") == [3.0, 3.0, 2.0, 2.0]
    @test_throws ErrorException Launcher.aggregate_metric(path, "daily_deaths")
  end
end

@testset "daily counts handle boundaries and invalid times" begin
  times = Union{Missing,Float64}[0.0, 0.99, 1.0, 2.0, missing, NaN, Inf]
  @test Launcher.daily!(zeros(Int, 2), times) == [2, 1]

  events = [(time=0.5,), (time=1.5,), (time=3.0,)]
  @test Launcher.daily!(event -> event.time, zeros(Int, 2), events) == [1, 1]
end

@testset "scoring helpers" begin
  @test Launcher.rolling_average([1, 2, 3, 4], 2) == [1.0, 1.5, 2.5, 3.5]
  @test Launcher.rolling_average(Int[], 3) == Float64[]
  @test Launcher.rmse_window([1.0], [1.0, 2.0], 1, 2) == sqrt(2.0)
  @test Launcher.rmse_window([9.0], [1.0], 3, 4) |> isnan
  @test Launcher.mean_value(Float64[]) |> isnan

  mktempdir() do directory
    path = joinpath(directory, "observed.csv")
    write(path, "date,value,notes\n2021-01-01,1.5,ok\n2021-01-02,invalid,skip\n2021-01-03,2.5,ok\n")
    @test Launcher.load_ground_truth_series(path) == [1.5, 2.5]
    @test Launcher.load_ground_truth_series(path; column=1, start_row=2) == Float64[]
  end
end

@testset "configuration scoring" begin
  mktempdir() do directory
    daily_path = joinpath(directory, "daily.jld2")
    jldopen(daily_path, "w") do file
      trajectory = JLD2.Group(file, "1")
      trajectory["daily_detections"] = [1.0, 2.0, 3.0]
      trajectory["daily_deaths"] = [2.0, 2.0, 2.0]
    end
    detections_path = joinpath(directory, "detections.csv")
    deaths_path = joinpath(directory, "deaths.csv")
    write(detections_path, "day,value\n1,1.0\n2,1.5\n3,2.0\n")
    write(deaths_path, "day,value\n1,1.0\n2,1.0\n3,1.0\n")

    result = Launcher.evaluate_config_window(
      daily_path,
      Dict{String,Any}(
        "daily_detections" => detections_path,
        "daily_deaths" => deaths_path,
      );
      weights=Dict("daily_detections" => 2.0, "daily_deaths" => 1.0),
      global_penalty=0.5,
      temporal_penalty=0.25,
    )

    @test result["daily_detections"]["rmse"] == 0.0
    @test result["daily_deaths"]["rmse"] == 1.0
    @test result["combined_rmse"] ≈ 1 / 3
    @test result["combined_score"] ≈ 1 / 3 + 0.75
    @test result["window"] == Dict("start_day" => 1, "end_day" => typemax(Int))
  end
end

@testset "command-line parsing" begin
  args = Launcher.parse_commandline([
    "config.toml",
    "--output-summary", "summary.jld2",
    "--output-checkpoint", "checkpoint.jld2",
  ])

  @test args["CONFIG"] == "config.toml"
  @test args["output-summary"] == "summary.jld2"
  @test args["output-checkpoint"] == "checkpoint.jld2"
  @test args["output-daily"] === nothing
end

@testset "output selection" begin
  mktempdir() do directory
    args = Dict{String,Any}(
      "output-summary" => joinpath(directory, "summary.jld2"),
      "output-daily" => joinpath(directory, "daily.jld2"),
      "output-params-dump" => joinpath(directory, "params.jld2"),
      "output-run-dump-prefix" => joinpath(directory, "run"),
    )

    outputs = Launcher.make_outputs(args, 1)

    @test Set(typeof.(outputs)) == Set([
      Launcher.Summary,
      Launcher.DailyTrajectories,
      Launcher.ParamsDump,
      Launcher.RunDump,
    ])
    @test only(output for output in outputs if output isa Launcher.Summary).filename ==
      args["output-summary"]
    @test only(output for output in outputs if output isa Launcher.ParamsDump).path ==
      args["output-params-dump"]
    @test only(output for output in outputs if output isa Launcher.RunDump).path_prefix ==
      args["output-run-dump-prefix"]

    close(only(output for output in outputs if output isa Launcher.DailyTrajectories).file)
  end
end

@testset "single-trajectory launch and output round trips" begin
  mktempdir() do directory
    population_path = write_test_population(joinpath(directory, "population.jld2"))
    config_path = joinpath(directory, "config.toml")
    summary_path = joinpath(directory, "summary.jld2")
    daily_path = joinpath(directory, "daily.jld2")
    open(config_path, "w") do io
      TOML.print(io, minimal_config(population_path))
    end

    Launcher.launch([
      config_path,
      "--output-summary", summary_path,
      "--output-daily", daily_path,
    ])

    jldopen(summary_path, "r") do file
      @test length(file["num_infections"]) == 1
      @test length(file["last_infections"]) == 1
      @test length(file["peak_daily_infections"]) == 1
      @test length(file["peak_daily_detections"]) == 1
      @test file["num_infections"][1] >= 1
    end
    jldopen(daily_path, "r") do file
      @test collect(keys(file)) == ["1"]
      trajectory = file["1"]
      @test haskey(trajectory, "daily_infections")
      @test haskey(trajectory, "daily_detections")
      @test haskey(trajectory, "daily_deaths")
      @test haskey(trajectory, "daily_hospitalizations")
      @test sum(trajectory["daily_infections"]) >= 1
    end
  end
end
