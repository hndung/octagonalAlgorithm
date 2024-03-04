# import Pkg; Pkg.add("CGAL")
# import Pkg; Pkg.add("IterTools")

using Random
using Distributions
using JLD2
using CGAL
using IterTools

include("utils.jl")
include("quickHull.jl")
include("octagonal.jl")
include("hexadecagonal.jl")

function main()
    # sizes = [100000, 1000000]
    sizes = [1000, 3000, 10000, 30000, 100000, 300000, 1000000, 3000000, 10000000, 30000000, 100000000]

    setNumbers = 10

    resultDirectory = "result/"

    # set benchmarking = true if want to benchmark
    benchmarking = true
    # only export results if exportResult = true and benchmarking = false
    exportResult = true

    # dataType: 1 for rectangle type, 2 for hollow rectangle type
    dataType = 1
    dataTypeName = ["rectangleInt", "hollowRectangleInt"]
    # input for rectangle data type
    # points in [x_min,x_max]x[y_min,y_max]
    x_min = -10000
    x_max = 10000
    y_min = -10000
    y_max = 10000

    x_min = -7500
    x_max = 7500
    y_min = -1500
    y_max = 1500
    # input for hollow rectangle data type
    # points in [x_min,x_max]x[y_min,y_max] \ [z_min,z_max]x[t_min,t_max]
    # [-7000, 7000]*[-1400,1400] is about 0.8711 the area of
    # [-7500, 7500]*[-1500,1500]
    z_min = -8000
    z_max = 8000
    t_min = -8000
    t_max = 8000

    instanceNames = Vector{String}(undef, length(sizes)*setNumbers)

    # running times
    runningTimeQH = Matrix{Float32}(undef, length(instanceNames), 7)
    runningTimeOct = Matrix{Float32}(undef, length(instanceNames), 7)
    runningTimeHex = Matrix{Float32}(undef, length(instanceNames), 7)
    runningTimeAklToussaint = Matrix{Float32}(undef, length(instanceNames), 7)
    runningTimeBykat = Matrix{Float32}(undef, length(instanceNames), 7)
    runningTimeEddy = Matrix{Float32}(undef, length(instanceNames), 7)
    runningTimeGrahamAndrew = Matrix{Float32}(undef, length(instanceNames), 7)

    Random.seed!(42)
    for i in 1:length(sizes)
        for j in 1:setNumbers
            k = (i-1)*setNumbers+j
            instanceNames[k] = string(dataTypeName[dataType], "_", sizes[i], "_", j)

            println()
            println("Consider instance ", instanceNames[k])

            # create random data
            points = Matrix{Float64}(undef, sizes[i], 2)
            if (dataType == 1)
                points = createRectangleIntData(sizes[i], x_min,x_max,y_min,y_max)
            #else
            #    points = createHollowRectangleData(sizes[i], x_min, x_max, y_min, y_max, z_min, z_max, t_min, t_max)
            end

            # create input for CGAL's algorithm from points
            points_CGAL = Vector{Point2}(undef, size(points, 1))
            for i in 1:length(points_CGAL)
                points_CGAL[i] = Point2(points[i,1],points[i,2])
            end

            if benchmarking
                # QuickHull
                # https://discourse.julialang.org/t/output-of-benchmark-to-string-or-table/27977/4
                # call atmost 30 times until execution time of 10000 is reached
                bmQH = run(@benchmarkable callQHull($points) samples=30 seconds=10000)
                # create report ignoring the worst 10 runs
                runningTimeQH[k,:] = report(bmQH, 10)

                # Octagonal algorithm
                bmOct = run(@benchmarkable octagonalCH($points) samples=30 seconds=10000)
                runningTimeOct[k,:] = report(bmOct, 10)

                # hexadecagonal algorithm
                bmHex = run(@benchmarkable hexadecagonalCH($points) samples=30 seconds=10000)
                runningTimeHex[k,:] = report(bmHex, 10)

                # Akl Toussaint
                bmAklToussaint = run(@benchmarkable ch_akl_toussaint($points_CGAL) samples=30 seconds=10000)
                runningTimeAklToussaint[k,:] = report(bmAklToussaint, 10)

                # Bykat
                bmBykat = run(@benchmarkable ch_bykat($points_CGAL) samples=30 seconds=10000)
                runningTimeBykat[k,:] = report(bmBykat, 10)

                # Eddy
                bmEddy = run(@benchmarkable ch_eddy($points_CGAL) samples=30 seconds=10000)
                runningTimeEddy[k,:] = report(bmEddy, 10)

                # Graham Andrew
                bmGrahamAndrew = run(@benchmarkable ch_graham_andrew($points_CGAL) samples=30 seconds=10000)
                runningTimeGrahamAndrew[k,:] = report(bmGrahamAndrew, 10)
            else
                # QuickHull
                exportFileQH = string(resultDirectory, instanceNames[k], "_QHull")
                callQHull(points, exportResult, exportFileQH)

                # Octagonal algorithm
                exportFileOct = string(resultDirectory, instanceNames[k], "_Octagonal")
                octagonalCH(points, exportResult, exportFileOct)

                # Hexadecagonal algorithm
                exportFileHex = string(resultDirectory, instanceNames[k], "_Hexadecagonal")
                hexadecagonalCH(points, exportResult, exportFileHex)
            end
        end
    end

    if benchmarking
        baseName = string(resultDirectory, dataTypeName[dataType], "_")
        # QuickHull
        exportReport(instanceNames, runningTimeQH, string(baseName,"qHull_running_time.csv"))
        # Octagonal algorithm
        exportReport(instanceNames, runningTimeOct, string(baseName,"octagonal_running_time.csv"))
        # Hexadecagonal algorithm
        exportReport(instanceNames, runningTimeHex, string(baseName,"hexadecagonal_running_time.csv"))
        # Akl Toussaint
        exportReport(instanceNames, runningTimeAklToussaint, string(baseName,"Akl_Toussaint_running_time.csv"))
        # Bykat
        exportReport(instanceNames, runningTimeBykat, string(baseName,"Bykat_running_time.csv"))
        # Eddy
        exportReport(instanceNames, runningTimeEddy, string(baseName,"Eddy_running_time.csv"))
        # Graham Andrew
        exportReport(instanceNames, runningTimeGrahamAndrew, string(baseName,"Graham_Andrew_running_time.csv"))
    end
end

main()

# Print the timings in the default way
# show(TimeOutputOct)
# show(TimeOutputQH)
