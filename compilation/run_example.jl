push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using MocosSimLauncher

project_dir = joinpath(@__DIR__, "..")
cd(project_dir) do
  launch(["--help"])
  json_path = joinpath("example", "example.json")
  args = String[json_path]
  launch(args)
end
