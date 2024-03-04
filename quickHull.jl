# install if needed
# using Pkg; Pkg.add("QHull");
# using Pkg; Pkg.add("TimerOutputs"); Pkg.add("BenchmarkTools");

using TimerOutputs
using QHull
using BenchmarkTools

include("utils.jl")

# Create a TimerOutput, this is the main type that keeps track of everything.
const TimeOutputQH = TimerOutput()

function callQHull(points, exportCH =false, exportFile="")
#    println("Calling QHull")
    @timeit TimeOutputQH "calling QHull" begin
        convexHull = chull(points)
    end

    if (exportCH)
        exportResult(points, convexHull.vertices, exportFile)
    end

end
