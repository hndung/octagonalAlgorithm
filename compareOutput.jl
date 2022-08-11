using JLD2

function arrayCompare(first, second)
    first_set = Set(first)
    second_set = Set(second)

    return issetequal(first_set, second_set)
end

function main()
    sizes = [100000, 1000000, 10000000]
    # sizes = [1000, 3000, 10000, 30000, 100000, 300000, 1000000, 3000000, 10000000, 30000000, 100000000]

    setNumbers = 3
    resultDirectory = "result/"

    for i in 1:length(sizes)
        for j in 1:setNumbers
            println()
            instanceName = string(sizes[i], "_", j)
            octagonal_file = string(resultDirectory, instanceName, "_Octagonal.jld2")
            qHull_file = string(resultDirectory, instanceName, "_QHull.jld2")

            octagonal_vertices = load(octagonal_file, "vertices")
            qHull_vertices = load(qHull_file, "vertices")

            print("Checking ", instanceName)
            print(", results are identical: ", arrayCompare(octagonal_vertices, qHull_vertices))
        end
    end
end

main()
