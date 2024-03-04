using JLD2

function arrayCompare(first, second)
    first_set = Set(first)
    second_set = Set(second)

    return issetequal(first_set, second_set)
end

function main()
    # sizes = [10000, 100000, 1000000, 10000000]
    sizes = [1000, 3000, 10000, 30000, 100000, 300000, 1000000, 3000000, 10000000, 30000000, 100000000]

    setNumbers = 10
    resultDirectory = "result/"

    # 1 for float random number in rectangle
    # 2 for float random number in hollow rectangle
    # 3 for int random number in rectangle
    # 4 for int random number in hollow rectangle
    dataType = 3
    dataTypeName = ["rectangle", "hollowRectangle", "rectangleInt", "hollowRectangleInt"]

    for i in 1:length(sizes)
        for j in 1:setNumbers
            println()
            instanceName = string(dataTypeName[dataType], "_", sizes[i], "_", j)
            octagonal_file = string(resultDirectory, instanceName, "_Octagonal.jld2")
            qHull_file = string(resultDirectory, instanceName, "_QHull.jld2")

            octagonal_vertices = load(octagonal_file, "vertices")
            qHull_vertices = load(qHull_file, "vertices")

            # note that this check only compares indices of the vertices
            # hence in case there are identical points having diffenrent indices
            # it may return false even the two set of points are identical
            println("Checking ", instanceName)
            println("\t OCT vs QH, results are identical: ", arrayCompare(octagonal_vertices, qHull_vertices))
            if (dataType==3 || dataType==4)
                hexadecagonal_file = string(resultDirectory, instanceName, "_Hexadecagonal.jld2")
                hexadecagonal_vertices = load(hexadecagonal_file, "vertices")
                println("\t HEX vs QH, results are identical: ", arrayCompare(hexadecagonal_vertices, qHull_vertices))
            end
        end
    end
end

main()
