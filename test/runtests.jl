using JLD2
using MocosSimLauncher
using Test

const Launcher = MocosSimLauncher

@testset "callback reset" begin
  callback = Launcher.DetectionCallback(2)
  callback.tracing_times[1] = 3.0
  callback.tracing_sources[1] = 2
  callback.tracing_types[1] = 1

  Launcher.reset!(callback)

  @test all(ismissing, callback.tracing_times)
  @test all(iszero, callback.tracing_sources)
  @test all(iszero, callback.tracing_types)
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
  end
end

@testset "daily counts include all finite infection times" begin
  @test Launcher.daily(skipmissing(Union{Missing,Float64}[0.0, missing, 1.0]), 2) == [1, 1]
end
