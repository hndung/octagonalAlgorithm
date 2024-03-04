# using Pkg; Pkg.add("JLD2"); Pkg.add("Printf"); Pkg.add("CSV");
# using Pkg; Pkg.add("DataFrames"); Pkg.add("StatsBase");

using JLD2
using Printf
using CSV
using DataFrames
using StatsBase
using Base.Threads

function loadData(filename::String)
    println("Loading jld2 file")
    points = Matrix{Float64}[]
    points = load(filename, "points")
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

function mySplit(startIndex, endIndex, numberOfSet)
    step = Int(ceil((endIndex-startIndex+1)/numberOfSet))
    grid = collect(startIndex:step:endIndex)
    push!(grid, endIndex+1)
end

function outsidePoints(points::Matrix{Float64},
                        u::Int,
                        v::Int,
                        startIndex::Int,
                        endIndex::Int)
    if (points[u]==points[v])
        return [u]
    end
    a = points[v,2] - points[u,2]
    b = points[u,1] - points[v,1]
    c = a*points[u,1] + b*points[u,2]
    outside_points = Vector{Int}(undef, 0)
    for i in startIndex:endIndex
        if (a*points[i,1]+b*points[i,2]>=c)
            push!(outside_points, i)
        end
    end
    # Include u and v in case of numerical errors that exclude those points
    if (u>=startIndex && u<=endIndex && a*points[u,1]+b*points[u,2]<c)
        push!(outside_points, u)
    end
    if (v>=startIndex && v<=endIndex && a*points[v,1]+b*points[v,2]<c)
        push!(outside_points, v)
    end
    return outside_points
end

function outsidePoints(points::Matrix{Float64}, u::Int, v::Int, splitFactor::Int)
    if (points[u]==points[v])
        return [u]
    end
    n = size(points, 1)
    threadNr = nthreads()
    if (threadNr == 1)
        return outsidePoints(points, u, v, 1, n)
    else
        grid = mySplit(1, n, threadNr*splitFactor)
        result = Vector{Vector{Int}}(undef, length(grid)-1)
        @threads for i in 1:length(grid)-1
            result[i] = outsidePoints(points, u, v, grid[i], grid[i+1]-1)
        end
        outside_points = Vector{Int}(undef, 0)
        for i in 1:length(grid)-1
            append!(outside_points, result[i])
        end
        return outside_points
    end
end

function createRectangleData(n::Int64, x_min::Number, x_max::Number,
                                        y_min::Number, y_max::Number)
    @assert x_min<x_max
    @assert y_min<y_max

    x_mean = (x_max+x_min)/2
    y_mean = (y_max+y_min)/2
    ratio = (x_max-x_min)/(y_max-y_min)

    points = rand(Uniform(y_min-y_mean, y_max-y_mean), n, 2)

    for i in 1:n
        points[i,1] *= ratio
        points[i,1] += x_mean
        points[i,2] += y_mean
    end

    return points
end

# points in the hollow rectangle
#          [x_min,x_max]x[y_min,y_max] \ [z_min,z_max]x[t_min,t_max]
function createHollowRectangleData(n::Int64, x_min::Number, x_max::Number,
                                            y_min::Number, y_max::Number,
                                            z_min::Number, z_max::Number,
                                            t_min::Number, t_max::Number)

    @assert x_min<z_min
    @assert z_min<z_max
    @assert z_max<x_max
    @assert y_min<t_min
    @assert t_min<t_max
    @assert t_max<y_max

    bigArea = (x_max-x_min)*(y_max-y_min)
    smallArea = (z_max-z_min)*(t_max-t_min)

    ratio = bigArea/(bigArea-smallArea)

    if (n>100000)
        ratio *= 1.1
    elseif (n>=10000)
        ratio *= 1.2
    else
        ratio *= 1.5
    end

    while(true)
        m = Int64(ceil(n*ratio))
        points_init = createRectangleData(m, x_min, x_max, y_min, y_max)
        points = Matrix{Float64}(undef, n, 2)
        pos = 1
        for i in 1:m
            if (pos>n)
                break
            end
            if (points_init[i,1]<= z_min || points_init[i,1] >= z_max || points_init[i,2]<= t_min || points_init[i,2] >= t_max)
                points[pos,:] = points_init[i,:]
                pos += 1
            end
        end

        if (pos>n)
            return points
        end
        ratio *= 1.1
    end
end

function createRectangleIntData(n::Int64, x_min::Int64, x_max::Int64, y_min::Int64, y_max::Int64)
    @assert x_min<x_max
    @assert y_min<y_max

    yMaxIndex = y_max-y_min+1
    maxIndex = (x_max-x_min+1)*yMaxIndex
    randomIndex = rand(1:maxIndex, n)

    p = collect(product(x_min:x_max, y_min:y_max))
    points = Matrix{Float64}(undef, n, 2)

    for i in 1:length(randomIndex)
        row = div(randomIndex[i]-1, yMaxIndex)+1
        column = mod(randomIndex[i]-1, yMaxIndex)+1
        points[i,1] = p[row,column][1]
        points[i,2] = p[row,column][2]
    end
    return points
end
