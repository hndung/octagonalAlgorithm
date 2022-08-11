# using Pkg; Pkg.add("JLD2"); Pkg.add("Printf"); Pkg.add("CSV");
# using Pkg; Pkg.add("DataFrames"); Pkg.add("StatsBase");

using JLD2
using Printf
using CSV
using DataFrames
using StatsBase

function loadData(filename::String)
    println("Loading jld2 file")
    points = Matrix{Float64}[]
    @timeit TimeOutput "load jld2 file" begin
        points = load(filename, "points")
    end
    return points
end

function report(bm, ignore::Int)
    # times in nano seconds hence multiplied by 1.e-9
    sorted_times = sort(bm.times)
    s = max(1, length(sorted_times)-ignore)
    min_time = minimum(sorted_times[1:s])*1.e-9
    max_time = maximum(sorted_times[1:s])*1.e-9
    mean_time = mean(sorted_times[1:s])*1.e-9
    geomean_time = geomean(sorted_times[1:s])*1.e-9
    median_time = median(sorted_times[1:s])*1.e-9
    runs = length(bm)
    @printf("running time (seconds):  min=%.5f", min_time)
    @printf(" max=%.5f", max_time)
    @printf(" mean=%.5f", mean_time)
    @printf(" geomean=%.5f", geomean_time)
    @printf(" median=%.5f", median_time)
    print(" runs=", runs)
    println(" ignore=", ignore)
    return [min_time, max_time, mean_time, geomean_time, median_time, runs, ignore]
end

function exportReport(names, runningTime, file_name)
    df = DataFrame(instance = names,
                    min = map(i -> runningTime[i,1], collect(1:length(names))),
                    max = map(i -> runningTime[i,2], collect(1:length(names))),
                    mean = map(i -> runningTime[i,3], collect(1:length(names))),
                    geomean = map(i -> runningTime[i,4], collect(1:length(names))),
                    median = map(i -> runningTime[i,5], collect(1:length(names))),
                    runs = map(i -> runningTime[i,6], collect(1:length(names))),
                    ignore = map(i -> runningTime[i,7], collect(1:length(names))))
    CSV.write(file_name, df)
end

function exportResult(points, vertices, exportFile)
    println("Writing jld2 and csv files of vertices")
    jldsave(string(exportFile, ".jld2"); vertices)
    df = DataFrame(index = vertices,
                    x = map(v -> points[v,1], vertices),
                    y = map(v -> points[v,2], vertices))
    CSV.write(string(exportFile, ".csv"), df)
end
