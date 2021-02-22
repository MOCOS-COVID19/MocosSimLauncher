push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using MocosSimLauncher

project_dir = joinpath(@__DIR__, "..")
cd(project_dir) do
  launch(["--help"])

  outdir = mktempdir()
  args = String[
    joinpath("example", "example.json"),
    "--output-summary", joinpath(outdir, "summary.jld2"),
    "--output-run-dump-prefix", joinpath(outdir, "run_dump"),
    "--output-daily", joinpath(outdir, "daily.jld2"),
    "--output-params-dump", joinpath(outdir, "params_dump.jld2"),
  ]
  launch(args)
end
