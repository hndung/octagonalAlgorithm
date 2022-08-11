# Octagonal cut algorithm
A fast algorithm for 2D convex hull problem

### Install libraries
- Open a terminal and go to the directory containing the code
- Run 
>        julia install.jl

### Run the programs
- Open a terminal and go to the directory containing the code 
- Run octogonal cut algorithm in sequential mode
>        julia mainOctagonal.jl
- Run octogonal cut algorithm in parallel mode
>        julia -t numberOfThreads mainOctagonal.jl
- Run quickhull
>        julia mainQH.jl

### Setting
- Benchmarking mode
> Set benchmarking = true in the main functions
- Export the convex hull to file
> Set benchmarking = false and exportResult = true in the main functions
