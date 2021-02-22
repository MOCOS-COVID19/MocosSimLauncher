using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()
using PackageCompiler
create_app(
  ".",
  joinpath(@__DIR__, "..", "build"),
  #precompile_execution_file=joinpath(@__DIR__, "run_example.jl"),
  filter_stdlibs=true,
  audit=true,
  force=true,
  )