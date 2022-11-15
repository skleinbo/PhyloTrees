include("moran_vs_coalescent.jl")

using GLMakie
d = 0.975

## Aggregate A,C over multiple coalescent runs.
Acoal = Int[]
Ccoal = Int[]
@time for _ in 1:10
  Pcoal = coalescent(n)
  _A, _C = treevalues!(Pcoal)
  append!(Acoal, _A)
  append!(Ccoal, _C)
end
dfcoal = to_mean_dataframe(Acoal, Ccoal)

## Aggregate A,C over multiple Birth/Death runs.
Abd = Int[]
Cbd = Int[]
@time for _ in 1:50
  global Pbd = birthdeath(1, ceil(Int, 1/(1-d))*n, d; N=0)
  _A, _C = treevalues!(Pbd[1])
  append!(Abd, _A)
  append!(Cbd, _C)
end
dfbd = to_mean_dataframe(Abd,Cbd)

## Aggregate A,C over multiple fluctuating_coalescent runs.
Afc = Int[]
Cfc = Int[]
@time for _ in 1:100
  global Pfc = fluctuating_coalescent(n)
  _A, _C = treevalues!(Pfc)
  append!(Afc, _A)
  append!(Cfc, _C)
end
dffc = to_mean_dataframe(Afc,Cfc)

## Comparison: unbalanced tree C~A*lnA
Punbalanced  = maximally_unbalanced(n÷2)
Aunb, Cunb = treevalues!(Punbalanced)

## Comparison: balanced tree C~A^2
Pbalanced  = maximally_balanced(log2(n))
Ab, Cb = treevalues!(Pbalanced)

## Plots
begin
  fig = Figure(resolution=72 .*(11,8), fontsize=12)
  colors = Makie.wong_colors()

  ax = Axis(fig[1,1], xscale=log10, yscale=log10,
    xlabel="A", ylabel="C/A")
  ax2 = Axis(fig[1,2], xscale=log10, yscale=identity,
    xlabel="A", ylabel="C/A")
  ylims!(ax2, 1e0,30)
  xlims!(ax2, 1e0,1e4)
  
  p11 = Makie.scatter!(ax, Aunb, Cunb./Aunb,
    label="Unbalanced", markersize=10, marker=:circle
  )
  p14 = Makie.scatter!(ax, dffc.a, dffc.covera,
    label="Fluc. Coal.", markersize=10, marker=:circle
  )
  p21 = Makie.scatter!(ax2, Ab, Cb./Ab,
    label="Balanced", markersize=10, marker=:star5,
    color=colors[3]
  )

  p22 = Makie.scatter!(ax2, dfbd.a, dfbd.covera,
    label="B/D", markersize=10, marker=:rect, color=(colors[1], 0.8))
  p23 = Makie.scatter!(ax2, dfcoal.a, dfcoal.covera,
    label="coalescent", markersize=10, marker=:rect,color=(colors[2], 0.5))
  
  
  η = 1.50
  p12 = Makie.lines!(ax,  1..1e2, A->A/4+1-1/4/A, label="A")
  p13 = Makie.lines!(ax,  1..1e4, A->A^(η-1), label="A^$(η-1)", color=:black, linestyle=:dash)
  p24 = Makie.lines!(ax2, 1..1e4, A->A^0.5, label="A^0.5", color=:black, linestyle=:dash)
  p25 = Makie.lines!(ax2, 1..1e4, A->1/A+(A+1)/A*(log(A+1)/log(2)-1), label="logA")
  
  Legend(fig[2,1], ax,  nbanks=2, orientation=:horizontal, tellheight=true)
  Legend(fig[2,2], ax2, nbanks=2, orientation=:horizontal, tellheight=true, tellwidth=true)
  colsize!(fig.layout, 1, Relative(1/2))
  colsize!(fig.layout, 2, Relative(1/2))

  current_figure()
end
