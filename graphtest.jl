include("moran_vs_coalescent.jl")

using GLMakie, GraphMakie, NetworkLayout

## -- Test -- ##
## Run a small simulation and plot the tree twice side by side 
## annotated with A and C respectively.
Pbd = birthdeath(1, 32^2, 0.9; N=0)
Pbal = maximally_unbalanced(6)

fig = Figure()
graphplot(fig[1,1], reverse(Pbd), nlabels=repr.(A(Pbd)), layout=Buchheim())
graphplot(fig[1,2], reverse(Pbd), nlabels=repr.(C(Pbd)), layout=Buchheim())
graphplot(fig[2,1], reverse(Pbal), nlabels=repr.(A(Pbal)), layout=Buchheim())
graphplot(fig[2,2], reverse(Pbal), nlabels=repr.(C(Pbal)), layout=Buchheim())
display(fig)
## ---------- ##
#
if !isinteractive()
    println("Press ENTER to quit.")
    readline()
end
