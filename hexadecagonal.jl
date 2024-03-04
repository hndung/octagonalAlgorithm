# install if needed
# using Pkg; Pkg.add("TimerOutputs");
# using Pkg; Pkg.add("QHull"); Pkg.add("BenchmarkTools");

using TimerOutputs
using QHull
using BenchmarkTools
using Base.Threads

include("utils.jl")

const splitFactorHex = 16
# Create a TimerOutput, this is the main type that keeps track of everything.
const TimeOutputHex = TimerOutput()

struct MinMaxType2
    min::Float64
    minIndex::Set{Int}
    max::Float64
    maxIndex::Set{Int}
end

function minMaxCollect(result::Vector{MinMaxType2})
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
    return MinMaxType2(mi, minIndex, ma, maxIndex)
end

function findMinMax2(v::Matrix{Float64},
                    coord::Int,
                    startIndex::Int,
                    endIndex::Int)
    mi = v[startIndex, coord]
    ma = v[startIndex, coord]
    minIndex = Set(startIndex)
    maxIndex = Set(startIndex)
    for i in startIndex+1:endIndex
        if (mi>v[i, coord])
            mi = v[i, coord]
            minIndex = Set(i)
        elseif (ma<v[i, coord])
            ma = v[i, coord]
            maxIndex = Set(i)
        else
            if (mi==v[i, coord])
                push!(minIndex, i)
            end
            if (ma==v[i, coord])
                push!(maxIndex, i)
            end
        end
    end
    return MinMaxType2(mi,minIndex,ma,maxIndex)
end

function findMinMax2(v::Matrix{Float64}, coord::Int)
    n = size(v, 1)
    threadNr = nthreads()
    if (threadNr == 1)
        return findMinMax2(v, coord, 1, n)
    else
        grid = mySplit(1, n, threadNr*splitFactorHex)
        result = Vector{MinMaxType2}(undef, length(grid)-1)
        @threads for i in 1:length(grid)-1
            result[i] = findMinMax2(v, coord, grid[i], grid[i+1]-1)
        end
        return minMaxCollect(result)
    end
end

function findDirectionalMinMax2(v::Matrix{Float64},
                                a::Vector{Float64},
                                startIndex::Int,
                                endIndex::Int)
    mi = a[1]*v[startIndex, 1] + a[2]*v[startIndex, 2]
    ma = mi
    minIndex = Set(startIndex)
    maxIndex = Set(startIndex)
    for i in startIndex+1:endIndex
        value = a[1]*v[i, 1] + a[2]*v[i, 2]
        if (mi>value)
            mi = value
            minIndex = Set(i)
        elseif (ma<value)
            ma = value
            maxIndex = Set(i)
        else
            if (mi==value)
                push!(minIndex, i)
            end
            if (ma==value)
                push!(maxIndex, i)
            end
        end
    end
    return MinMaxType2(mi, minIndex, ma, maxIndex)
end

function findDirectionalMinMax2(v::Matrix{Float64}, a::Vector{Float64})
    n = size(v, 1)
    threadNr = nthreads()
    if (threadNr == 1)
        return findDirectionalMinMax2(v, a, 1, n)
    else
        grid = mySplit(1, n, threadNr*splitFactorHex)
        result = Vector{MinMaxType2}(undef, length(grid)-1)
        @threads for i in 1:length(grid)-1
            result[i] = findDirectionalMinMax2(v, a, grid[i], grid[i+1]-1)
        end
        return minMaxCollect(result)
    end
end

function findCrossDirectionalMinMax2(v::Matrix{Float64},
                                    a::Vector{Float64},
                                    startIndex::Int,
                                    endIndex::Int)
    min1 = a[1]*v[startIndex, 1] + a[2]*v[startIndex, 2]
    min2 = a[2]*v[startIndex, 2] - a[1]*v[startIndex, 1]
    max1 = min1
    max2 = min2

    minIndex1 = Set(startIndex); minIndex2 = Set(startIndex)
    maxIndex1 = Set(startIndex); maxIndex2 = Set(startIndex)
    for i in startIndex+1:endIndex
        c = a[1]*v[i, 1]
        d = a[2]*v[i, 2]
        value1 = c+d
        value2 = d-c
        if (min1>value1)
            min1 = value1
            minIndex1 = Set(i)
        elseif (max1<value1)
            max1 = value1
            maxIndex1 = Set(i)
        else
            if (min1==value1)
                push!(minIndex1, i)
            end
            if (max1==value1)
                push!(maxIndex1, i)
            end
        end
        if (min2>value2)
            min2 = value2
            minIndex2 = Set(i)
        elseif (max2<value2)
            max2 = value2
            maxIndex2 = Set(i)
        else
            if (min2==value2)
                push!(minIndex2, i)
            end
            if (max2==value2)
                push!(maxIndex2, i)
            end
        end
    end
    return [MinMaxType2(min1, minIndex1, max1, maxIndex1); MinMaxType2(min2, minIndex2, max2, maxIndex2)]
end

function findCrossDirectionalMinMax2(v::Matrix{Float64}, a::Vector{Float64})
    n = size(v, 1)
    threadNr = nthreads()
    if (threadNr == 1)
        return findCrossDirectionalMinMax2(v, a, 1, n)
    else
        grid = mySplit(1, n, threadNr*splitFactorHex)
        result = Matrix{MinMaxType2}(undef, length(grid)-1, 2)
        @threads for i in 1:length(grid)-1
            result[i,:] = findCrossDirectionalMinMax2(v, a, grid[i], grid[i+1]-1)
        end
        return [minMaxCollect(result[:,1]); minMaxCollect(result[:,2])]
    end
end

function getLineVertices(points::Matrix{Float64}, S::Set{Int}, coord::Int)
    @assert length(S)>=1
    i = pop!(S)
    if (length(S)==0)
        return [i,i]
    end
    minVal = points[i,coord]
    minIndex = i
    maxVal = minVal
    maxIndex = i
    for j in S
        if (points[j,coord] < minVal)
            minVal = points[j,coord]
            minIndex = j
        elseif (points[j,coord] > maxVal)
            maxVal = points[j,coord]
            maxIndex = j
        end
    end
    return [minIndex,maxIndex]
end

function find16ExtremePoints(points::Matrix{Float64})
    @timeit TimeOutputHex "find 16 EP" begin
        @timeit TimeOutputHex "first 8 EP" begin
            horizontalMinMax = findMinMax2(points, 1)
            minX = horizontalMinMax.min
            maxX = horizontalMinMax.max
            v1 = getLineVertices(points, horizontalMinMax.maxIndex, 2)
            v1first = v1[1]
            v1last = v1[2]
            v5 = getLineVertices(points, horizontalMinMax.minIndex, 2)
            v5first = v5[2]
            v5last = v5[1]

            verticalMinMax = findMinMax2(points, 2)
            minY = verticalMinMax.min
            maxY = verticalMinMax.max
            v3 = getLineVertices(points, verticalMinMax.maxIndex, 1)
            v3first = v3[2]
            v3last = v3[1]
            v7 = getLineVertices(points, verticalMinMax.minIndex, 1)
            v7first = v7[1]
            v7last = v7[2]
        end

        # calling findCrossDirectionalMinMax2 is a little bit faster than
        # call two times findDirectionalMinMax2 (see above)
        @timeit TimeOutputHex "second 8 EP" begin
            diagonalMinMax = findCrossDirectionalMinMax2(points, [maxY - minY, maxX - minX])
            v2 = getLineVertices(points, diagonalMinMax[1].maxIndex, 1)
            v2first = v2[2]
            v2last = v2[1]
            v6 = getLineVertices(points, diagonalMinMax[1].minIndex, 1)
            v6first = v6[1]
            v6last = v6[2]
            v4 = getLineVertices(points, diagonalMinMax[2].maxIndex, 1)
            v4first = v4[2]
            v4last = v4[1]
            v8 = getLineVertices(points, diagonalMinMax[2].minIndex, 1)
            v8first = v8[1]
            v8last = v8[2]
        end
    end
    return [v1last, v2first, v2last, v3first, v3last, v4first, v4last, v5first, v5last, v6first, v6last, v7first, v7last, v8first, v8last, v1first]
end

function filterRemainingPoints2(points::Matrix{Float64},
                                extremalVertices::Vector{Int})
#    println("Get remaining points")
    @assert length(extremalVertices)==16
    @timeit TimeOutputHex "filter remaining points" begin
        remaining_points = Vector{Int}(undef, 0)
        for i in 1:8
            u = extremalVertices[2*i-1]
            v = extremalVertices[2*i]
            append!(remaining_points, outsidePoints(points, u, v, splitFactorHex))
        end
    end
    return remaining_points
end

function hexadecagonalCH(points, exportCH = false, exportFile = "")
    v = find16ExtremePoints(points)
    remaining_points_index = filterRemainingPoints2(points, v)
    # println("remaining points ", length(remaining_points_index))
    remaining_points = Matrix{Float64}(undef, length(remaining_points_index), 2)
    @timeit TimeOutputHex "get remaining points" begin
        for i in 1:length(remaining_points_index)
            remaining_points[i,:] = points[remaining_points_index[i], :]
        end
    end

#    println("Calling QHull")
    @timeit TimeOutputHex "calling QHull" begin
        convex_hull = chull(remaining_points)
        # vertices_index = map(v -> remaining_points_index[v], convex_hull.vertices)
        # vertices = map(v -> points[v,:], vertices_index)
        # println("vertices Octagonal ", length(vertices))
    end

    if (exportCH)
        vertices = map(v -> remaining_points_index[v], convex_hull.vertices)
        exportResult(points, vertices, exportFile)
    end
end
