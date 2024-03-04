# Octagonal cut algorithm
A fast algorithm for 2D convex hull problem

### Install libraries
- Open a terminal and go to the directory containing the code
- Run 
>        julia install.jl

### Run the programs
- Open a terminal and go to the directory containing the code 
- Run all algorithms in sequential mode for real data
>        julia main.jl
- Run all algorithms in parallel mode for real data
>        julia -t numberOfThreads main.jl
- Run all algorithms in sequential mode for integer data
>        julia mainIntData.jl
- Run all algorithms in parallel mode for integer data
>        julia -t numberOfThreads mainIntData.jl
One can choose the data type (rectangle or hollow rectangle) by setting dataType equal 1 or 2 in the main files main.jl and mainIntData.jl.

### Setting
- Benchmarking mode
> Set benchmarking = true in the main functions
- Export the convex hull to file
> Set benchmarking = false and exportResult = true in the main functions
