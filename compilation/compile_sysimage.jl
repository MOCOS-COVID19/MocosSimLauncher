using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()
using PackageCompiler
create_sysimage(
  :MocosSimLauncher,
  sysimage_path=joinpath(@__DIR__, "..", "sysimage.img"),
  precompile_execution_file=joinpath(@__DIR__, "run_example.jl"),
  #precompile_statements_file=joinpath(@__DIR__, "precompile_statements.jl"),
  project=joinpath(@__DIR__, ".."),
  filter_stdlibs=true,
  incremental=false,
  )