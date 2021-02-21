#!/bin/bash
#julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --startup-file=no --trace-compile precompile_trace.jl run_example.jl