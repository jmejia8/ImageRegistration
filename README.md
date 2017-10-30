# Image Registration using Evolutionary Centers Algorithm

Experiments for finding affine transformations between two sets of points. 
I compare an Evolutionary Algorithm and a classic method RANSAC (written in MATLAB).

It is coded for Julia 0.6.0

## Instructions

First install dependencies:
- [Metaheuristics](https://github.com/jmejia8/Metaheuristics.jl) for evolutionary algorithms.
	```julia
		Pkg.clone("git@github.com:jmejia8/Metaheuristics.jl.git")
	```
- [PyPlot](https://github.com/JuliaPy/PyPlot.jl) for plotting in Julia.
	```julia
		Pkg.add("PyPlot")
	```

- [MATLAB](https://github.com/JuliaInterop/MATLAB.jl)  provides an interface for using [MATLAB™](https://www.mathworks.com/products/matlab.html) from Julia.
	```julia
		Pkg.add("MATLAB")
	```

## My First Time

After dependencies installation:

1. Open terminal.
2. cd to repo.
3. Write `include("main.jl")` in Julia REPL.
4. Imagine a better world.