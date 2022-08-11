# install if needed
# import Pkg; Pkg.add("Distributions")

using Random
using Distributions
using JLD2

include("octagonal.jl")

function main()
    sizes = [100000, 1000000, 10000000]
    # sizes = [1000, 3000, 10000, 30000, 100000, 300000, 1000000, 3000000, 10000000, 30000000, 100000000]

    setNumbers = 3

    resultDirectory = "result/"

    # set benchmarking = true if want to benchmark
    benchmarking = true
    # only export results if exportResult = true and benchmarking = false
    exportResult = true

    instanceNames = Vector{String}(undef, length(sizes)*setNumbers)
    runningTime = Matrix{Float32}(undef, length(instanceNames), 7)

    # create data, use the same seed number in mainOctagonal.jl and mainQH.jl
    Random.seed!(42)
    for i in 1:length(sizes)
        for j in 1:setNumbers
            k = (i-1)*setNumbers+j
            instanceNames[k] = string(sizes[i], "_", j)

            println()
            println("Consider instance ", instanceNames[k])
            points = rand(Uniform(-1500, 1500), sizes[i], 2)
            points[:,1] *= 5

            if benchmarking
                # https://discourse.julialang.org/t/output-of-benchmark-to-string-or-table/27977/4
                # call atmost 30 times until execution time of 10000 is reached
                bm = run(@benchmarkable octagonalCH($points) samples=30 seconds=10000)
                # create report ignoring the worst 10 runs
                runningTime[k,:] = report(bm, 10)
            else
                exportFile = string(resultDirectory, instanceNames[k], "_Octagonal")
                octagonalCH(points, exportResult, exportFile)
            end
        end
    end

    if benchmarking
        exportReport(instanceNames, runningTime, string(resultDirectory,"octagonal_running_time.csv"))
    end
end

main()

# Print the timings in the default way
show(TimeOutput)
