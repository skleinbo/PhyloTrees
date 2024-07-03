# Usage

1. Clone this repository.
2. Start Julia in the repo folder `julia --project`
3. **Either** install `WeightedSampling` directly from its repo
   - `] add https://github.com/skleinbo/WeightedSampling.jl.git`
   - or subscribe to my package registry
     1. `] registry add https://github.com/skleinbo/JuliaRegistry.jl/`
     2. `] add WeightedSampling`

4. Install the unregistered dependencies

    ```julia
    julia> ]add https://github.com/skleinbo/BinaryTrees.jl https://github.com/skleinbo/TreeProcesses.jl
    ```

# Example

Include `run.jl` to sample data and `plot.jl` to generate this plot from them

![test.png](images/test.png)
