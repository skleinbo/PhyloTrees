include("moran_vs_coalescent.jl")

using GLMakie, GraphMakie, NetworkLayout

## -- Test -- ##
## Run a small simulation and plot the tree twice side by side 
## annotated with A and C respectively.
Pbd = birthdeath(1, 32^2, 0.9; N=0)[1]
Pbal = maximally_unbalanced(6)

Abd,Cbd = treevalues!(Pbd)
Abal,Cbal = treevalues!(Pbal)


fig = Figure()
graphplot(fig[1,1], SimpleDiGraph(adjacency_matrix(Pbd)), nlabels=repr.(Abd), layout=Buchheim())
graphplot(fig[1,2], (Pbd), nlabels=repr.(Cbd), layout=Buchheim())
graphplot(fig[2,1], reverse(Pbal), nlabels=repr.(Abal), layout=Buchheim())
graphplot(fig[2,2], reverse(Pbal), nlabels=repr.(Cbal), layout=Buchheim())
display(fig)
## ---------- ##
#
if !isinteractive()
    println("Press ENTER to quit.")
    readline()
end
