using DelimitedFiles

function load_ground_truth_series(path::AbstractString; column::Integer=2, start_row::Integer=2)
  mat = readdlm(path, ',', Any, '\n')  # 2D matrix: rows × cols
  values = Float64[]
  for i in start_row:size(mat, 1)
    size(mat, 2) < column && continue
    value = tryparse(Float64, string(mat[i, column]))
    value === nothing && continue
    push!(values, value)
  end
  values
end

function rolling_average(values::AbstractVector{<:Real}, window::Integer=7)
  result = zeros(Float64, length(values))
  acc = 0.0
  for idx in eachindex(values)
    acc += values[idx]
    if idx > window
      acc -= values[idx - window]
    end
    result[idx] = acc / min(idx, window)
  end
  result
end

function aggregate_metric(path::AbstractString, metric::AbstractString)
  agg = Float64[]
  count = 0
  jldopen(path, "r") do f
    for key in keys(f)
      group = f[key]
      haskey(group, metric) || continue
      data = Float64.(group[metric])
      rolled = rolling_average(data)
      if isempty(agg)
        resize!(agg, length(rolled))
        fill!(agg, 0.0)
      elseif length(rolled) > length(agg)
        old_length = length(agg)
        resize!(agg, length(rolled))
        fill!(@view(agg[old_length + 1:end]), 0.0)
      end
      limit = min(length(agg), length(rolled))
      @inbounds for idx in 1:limit
        agg[idx] += rolled[idx]
      end
      count += 1
    end
  end
  count == 0 && error("Metric $(metric) not found in $(path)")
  agg ./ count
end

mean_value(values::AbstractVector{<:Real}) = isempty(values) ? NaN : sum(values) / length(values)

function rmse_window(predicted::AbstractVector{<:Real}, observed::AbstractVector{<:Real}, start_day::Integer, end_day::Integer)
  start_idx = max(1, start_day)
  end_idx   = min(end_day, length(observed))
  end_idx < start_idx && return NaN
  # pad predicted with zeros if simulation ended early (e.g. epidemic died out)
  pred = if length(predicted) < end_idx
    vcat(predicted, zeros(Float64, end_idx - length(predicted)))
  else
    predicted
  end
  diff = @view(pred[start_idx:end_idx]) .- @view(observed[start_idx:end_idx])
  sqrt(mean_value(diff .^ 2))
end

function evaluate_config_window(daily_path::AbstractString, ground_truth::Dict{String,Any}; start_day::Integer=1, end_day::Integer=typemax(Int), weights::Dict{String,Float64}=Dict{String,Float64}(), global_penalty::Float64=0.0, temporal_penalty::Float64=0.0)
  metrics = Dict(
    "daily_detections" => get(ground_truth, "daily_detections", nothing),
    "daily_deaths" => get(ground_truth, "daily_deaths", nothing),
    "daily_hospitalizations" => get(ground_truth, "daily_hospitalizations", nothing),
  )
  results = Dict{String,Any}()
  weighted_sum = 0.0
  total_weight = 0.0
  for (metric, gt_path) in metrics
    isnothing(gt_path) && continue
    observed = load_ground_truth_series(string(gt_path))
    predicted = aggregate_metric(daily_path, metric)
    score = rmse_window(predicted, observed, start_day, end_day)
    weight = get(weights, metric, 1.0)
    results[metric] = Dict("rmse" => score, "weight" => weight)
    if !isnan(score)
      weighted_sum += score * weight
      total_weight += weight
    end
  end
  combined = total_weight == 0 ? NaN : weighted_sum / total_weight
  results["combined_rmse"] = combined
  results["combined_score"] = isnan(combined) ? NaN : combined + global_penalty + temporal_penalty
  results["global_penalty"] = global_penalty
  results["temporal_penalty"] = temporal_penalty
  results["window"] = Dict("start_day" => start_day, "end_day" => end_day)
  results
end
