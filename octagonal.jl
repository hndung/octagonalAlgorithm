# install if needed
# using Pkg; Pkg.add("TimerOutputs");
# using Pkg; Pkg.add("QHull"); Pkg.add("BenchmarkTools");

using TimerOutputs
using QHull
using BenchmarkTools
using Base.Threads

include("utils.jl")

const splitFactor = 4
const TimeOutput = TimerOutput()

struct MinMaxType
    min::Float64
    minIndex::Int
    max::Float64
    maxIndex::Int
end

function mySplit(startIndex, endIndex, numberOfSet)
    step = Int(ceil((endIndex-startIndex+1)/numberOfSet))
    grid = collect(startIndex:step:endIndex)
    push!(grid, endIndex+1)
end

function minMaxCollect(result::Vector{MinMaxType})
    mi = result[1].min
    ma = result[1].max
    minIndex = result[1].minIndex
    maxIndex = result[1].maxIndex
    for i in 2:length(result)
        if (mi>result[i].min)
            mi = result[i].min
            minIndex = result[i].minIndex
        end
        if (ma<result[i].max)
            ma = result[i].max
            maxIndex = result[i].maxIndex
        end
    end
    return MinMaxType(mi, minIndex, ma, maxIndex)
end

function findMinMax(v::Matrix{Float64},
                    coord::Int,
                    startIndex::Int,
                    endIndex::Int)
    mi = v[startIndex, coord]
    ma = v[startIndex, coord]
    minIndex = startIndex
    maxIndex = startIndex
    for i in startIndex+1:endIndex
        if (mi>v[i, coord])
            mi = v[i, coord]
            minIndex = i
        elseif (ma<v[i, coord])
            ma = v[i, coord]
            maxIndex = i
        end
    end
    return MinMaxType(mi,minIndex,ma,maxIndex)
end

function findMinMax(v::Matrix{Float64}, coord::Int)
    n = size(v, 1)
    threadNr = nthreads()
    if (threadNr == 1)
        return findMinMax(v, coord, 1, n)
    else
        grid = mySplit(1, n, threadNr*splitFactor)
        result = Vector{MinMaxType}(undef, length(grid)-1)
        @threads for i in 1:length(grid)-1
            result[i] = findMinMax(v, coord, grid[i], grid[i+1]-1)
        end
        return minMaxCollect(result)
    end
end

function findDirectionalMinMax(v::Matrix{Float64},
                                a::Vector{Float64},
                                startIndex::Int,
                                endIndex::Int)
    mi = a[1]*v[startIndex, 1] + a[2]*v[startIndex, 2]
    ma = mi
    minIndex = startIndex
    maxIndex = startIndex
    for i in startIndex+1:endIndex
        value = a[1]*v[i, 1] + a[2]*v[i, 2]
        if (mi>value)
            mi = value
            minIndex = i
        elseif (ma<value)
            ma = value
            maxIndex = i
        end
    end
    return MinMaxType(mi, minIndex, ma, maxIndex)
end

function findDirectionalMinMax(v::Matrix{Float64}, a::Vector{Float64})
    n = size(v, 1)
    threadNr = nthreads()
    if (threadNr == 1)
        return findDirectionalMinMax(v, a, 1, n)
    else
        grid = mySplit(1, n, threadNr*splitFactor)
        result = Vector{MinMaxType}(undef, length(grid)-1)
        @threads for i in 1:length(grid)-1
            result[i] = findDirectionalMinMax(v, a, grid[i], grid[i+1]-1)
        end
        return minMaxCollect(result)
    end
end

function findCrossDirectionalMinMax(v::Matrix{Float64},
                                    a::Vector{Float64},
                                    startIndex::Int,
                                    endIndex::Int)
    min1 = a[1]*v[startIndex, 1] + a[2]*v[startIndex, 2]
    min2 = a[2]*v[startIndex, 2] - a[1]*v[startIndex, 1]
    max1 = min1
    max2 = min2

    minIndex1 = startIndex; minIndex2 = startIndex
    maxIndex1 = startIndex; maxIndex2 = startIndex
    for i in startIndex+1:endIndex
        c = a[1]*v[i, 1]
        d = a[2]*v[i, 2]
        value1 = c+d
        value2 = d-c
        if (min1>value1)
            min1 = value1
            minIndex1 = i
        elseif (max1<value1)
            max1 = value1
            maxIndex1 = i
        end
        if (min2>value2)
            min2 = value2
            minIndex2 = i
        elseif (max2<value2)
            max2 = value2
            maxIndex2 = i
        end
    end
    return [MinMaxType(min1, minIndex1, max1, maxIndex1); MinMaxType(min2, minIndex2, max2, maxIndex2)]
end

function findCrossDirectionalMinMax(v::Matrix{Float64}, a::Vector{Float64})
    n = size(v, 1)
    threadNr = nthreads()
    if (threadNr == 1)
        return findCrossDirectionalMinMax(v, a, 1, n)
    else
        grid = mySplit(1, n, threadNr*splitFactor)
        result = Matrix{MinMaxType}(undef, length(grid)-1, 2)
        @threads for i in 1:length(grid)-1
            result[i,:] = findCrossDirectionalMinMax(v, a, grid[i], grid[i+1]-1)
        end
        return [minMaxCollect(result[:,1]); minMaxCollect(result[:,2])]
    end
end

function find8ExtremePoints(points::Matrix{Float64})
    @timeit TimeOutput "find 8 EP" begin
        @timeit TimeOutput "first 4 EP" begin
            horizontalMinMax = findMinMax(points, 1)
            minX = horizontalMinMax.min
            maxX = horizontalMinMax.max
            v1 = horizontalMinMax.maxIndex
            v5 = horizontalMinMax.minIndex

            verticalMinMax = findMinMax(points, 2)
            minY = verticalMinMax.min
            maxY = verticalMinMax.max
            v3 = verticalMinMax.maxIndex
            v7 = verticalMinMax.minIndex
        end

        # @timeit TimeOutput "second 4 EP" begin
        #     diagonalMinMax = findDirectionalMinMax(points, [maxY - minY, maxX - minX])
        #     v2 = diagonalMinMax.maxIndex
        #     v6 = diagonalMinMax.minIndex
        #     diagonalMinMax = findDirectionalMinMax(points, [minY - maxY, maxX - minX])
        #     v4 = diagonalMinMax.maxIndex
        #     v8 = diagonalMinMax.minIndex
        # end

        # calling findCrossDirectionalMinMax is a little bit faster than
        # call two times findDirectionalMinMax (see above)
        @timeit TimeOutput "second 4 EP" begin
            diagonalMinMax = findCrossDirectionalMinMax(points, [maxY - minY, maxX - minX])
            v2 = diagonalMinMax[1].maxIndex
            v6 = diagonalMinMax[1].minIndex
            v4 = diagonalMinMax[2].maxIndex
            v8 = diagonalMinMax[2].minIndex
        end
    end
    return [v1,v2,v3,v4,v5,v6,v7,v8,v1]
end

function outsidePoints(points::Matrix{Float64},
                        u::Int,
                        v::Int,
                        startIndex::Int,
                        endIndex::Int)
    if (u==v)
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

function outsidePoints(points::Matrix{Float64}, u::Int, v::Int)
    if (u==v)
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

function filterRemainingPoints(points::Matrix{Float64},
                                extremalVertices::Vector{Int})
#    println("Get remaining points")
    @timeit TimeOutput "filter remaining points" begin
        remaining_points = Vector{Int}(undef, 0)
        for i in 1:8
            u = extremalVertices[i]
            v = extremalVertices[i+1]
            append!(remaining_points, outsidePoints(points, u, v))
        end
    end
    return remaining_points
end

function octagonalCH(points, exportCH = false, exportFile = "")
    v = find8ExtremePoints(points)
    remaining_points_index = filterRemainingPoints(points, v)

    remaining_points = Matrix{Float64}(undef, length(remaining_points_index), 2)
    @timeit TimeOutput "get remaining points" begin
        for i in 1:length(remaining_points_index)
            remaining_points[i,:] = points[remaining_points_index[i], :]
        end
    end

#    println("Calling QHull")
    @timeit TimeOutput "calling QHull" begin
        convex_hull = chull(remaining_points)
    end

    if (exportCH)
        vertices = map(v -> remaining_points_index[v], convex_hull.vertices)
        exportResult(points, vertices, exportFile)
    end
end
