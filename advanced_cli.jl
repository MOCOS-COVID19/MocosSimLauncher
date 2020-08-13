#push!(LOAD_PATH, joinpath(@__DIR__,"Launcher"))
#push!(LOAD_PATH, joinpath(@__DIR__,"Simulation"))
#println(LOAD_PATH)
println("Launching")
using MocosSimLauncher

launch()
