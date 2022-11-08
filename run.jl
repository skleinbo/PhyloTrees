includet("moran_vs_coalescent.jl")

## -- Test -- ##
@time P = birthdeath(32, 2*32^2, 1.0; N=0)
P = maximally_unbalanced(16)
fig = Figure()
graphplot(fig[1,1], reverse(P), nlabels=repr.(A(P)), layout=Buchheim())
graphplot(fig[1,2], reverse(P), nlabels=repr.(C(P)), layout=Buchheim())
## ---------- ##

n = 2^12
d = 0.975

Acoal = Int[]
Ccoal = Int[]
@time for _ in 1:100
  Pcoal = coalescent(n)
  @assert ne(Pcoal) == nv(Pcoal) - 1
  append!(Acoal, A(Pcoal))
  append!(Ccoal, C(Pcoal))
end
Pbalanced  = maximally_balanced(log2(n))
Punbalanced  = maximally_unbalanced(n÷2)
Abd = Int[]
Cbd = Int[]
@time for _ in 1:50
  global Pbd = birthdeath(1, ceil(Int, 1/(1-d))*n, d; N=0)
  @assert ne(Pbd) == nv(Pbd) - 1
  append!(Abd, A(Pbd))
  append!(Cbd, C(Pbd))
end
function to_mean_dataframe(A,C)
  df = DataFrame(;A,C)
  combine(groupby(df, :A), :A=>mean=>:a, [:A, :C]=>((a,c)->mean(c./a))=>:covera)
end
# graphplot(Pcoal, names=vertices(P), method=:tree, nodesize=0.1, fontsize=12, nodecolor=[begin v in cg ? :blue : :orange end for v in vertices(P)])
Pmoran = moran2(n, n^2)
Pyule = yule(n)

# Acoal = A(Pcoal)
# Ccoal = C(Pcoal)
dfcoal = to_mean_dataframe(Acoal, Ccoal)
Amoran = A(Pmoran)
Cmoran = C(Pmoran)
Ayule = A(Pyule)
Cyule = C(Pyule)
# Abd = A(Pbd)
# Cbd = C(Pbd)
dfbd = to_mean_dataframe(Abd,Cbd)
Aunb = A(Punbalanced)
Cunb = C(Punbalanced)
Ab = A(Pbalanced)
Cb = C(Pbalanced)

idxcoal = findall(x->1<x<1e6, Acoal)
Acoal = Acoal[idxcoal]
Ccoal = Ccoal[idxcoal]

idxmoran = findall(x->1<x<1e8, Amoran)
Amoran = Amoran[idxmoran]
Cmoran = Cmoran[idxmoran]

idxyule = findall(x->1<x<1e6, Ayule)
Ayule = Ayule[idxyule]
Cyule = Cyule[idxyule]

idxbd = findall(x->3<x<1e6, Abd)
Abd = Abd[idxbd]
Cbd = Cbd[idxbd]
## Plots
begin
  fig = Figure(resolution=72 .*(11,8), fontsize=12)
  colors = Makie.wong_colors()
  ax = Axis(fig[1,1], xscale=log10, yscale=log10,
    xlabel="A", ylabel="C/A")
  p1 = Makie.scatter!(ax, Aunb, Cunb./Aunb,
    label="Unbalanced", markersize=10, marker=:circle
    )
  ax2 = Axis(fig[1,2], xscale=log10, yscale=identity,
    xlabel="A", ylabel="C/A")
  ylims!(ax2, 1e0,30)
  xlims!(ax2, 1e0,1e4)
  p2 = Makie.scatter!(ax2, Ab, Cb./Ab,
    label="Balanced", markersize=10, marker=:star5,
    color=colors[3]
    )
  # Makie.scatter!(ax, Amoran, Cmoran./Amoran, label="Moran") 
  p3 = Makie.scatter!(ax2, dfbd.a, dfbd.covera,
    label="B/D", markersize=10, marker=:rect, color=(colors[1], 0.8))
  p32 = Makie.scatter!(ax2, dfcoal.a, dfcoal.covera,
    label="coalescent", markersize=10, marker=:rect,color=(colors[2], 0.5))

  # plot!(x->2*x^1.36, label="x^1.36")
  η = 1.50
  p4 = Makie.lines!(ax, 1..100,A->A/4+3/4/A, label="A")
  p5 = Makie.lines!(ax, 1..1e4,x->x^(η-1), label="A^$(η-1)", color=:black, linestyle=:dash)
  p6 = Makie.lines!(ax2, 1..1e4, x->x^0.5, label="A^0.5", color=:black, linestyle=:dash)
  p6 = Makie.lines!(ax2, 1..1e4, A->1/A+(A+1)/A*(log(A+1)/log(2)-1))
  # br = fig[2,1:2] = GridLayout()
  # linkyaxes!(ax,ax2)
  Legend(fig[2,1], ax, orientation=:horizontal, tellheight=true)
  Legend(fig[2,2], ax2, orientation=:horizontal, tellheight=true, tellwidth=true)
  colsize!(fig.layout, 1, Relative(1/2))
  colsize!(fig.layout, 2, Relative(1/2))
  # resize_to_layout!(fig)
  current_figure()
end