push!(LOAD_PATH, "./TreeProcesses/src/")
push!(LOAD_PATH, "./BinaryTrees/src/")
using BinaryTrees, TreeProcesses
using CategoricalArrays, DataFrames, DataFramesMeta, StatsBase
using JLD2

using GLMakie

n = 2^13
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
  Pbd = birthdeath(1, ceil(Int, 1/(1-d))*n, d; N=0)
  _A, _C = treevalues!(Pbd[1])
  append!(Abd, _A)
  append!(Cbd, _C)
end
dfbd = to_mean_dataframe(Abd,Cbd)

## Aggregate A,C over multiple fluctuating_coalescent runs.
Afc = Int[]
Cfc = Int[]
@time for _ in 1:50
  Pfc = fluctuating_coalescent(n; fuse=max)
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

##--------------------------------------##


dffc_eta = load("data/03.jld2", "dffc_eta")

using ThreadsX, StatsBase
import Base.Iterators: product
etas = [0.0,2.0,3.0,6.0,12.0]
etas = 1.2 #1.2:0.005:(1.3-0.005)
etas = [1.0, 1.2, 1.3, 1.4, 2.0]
ns = 2 .^[8,10,14,16,18]
dffc_eta = DataFrame()
(map(product(etas, ns)) do (η,n)
  Afc = Int[]
  Cfc = Int[]
  if η==0.0
    fuse = max
  else
    fuse=(x,y)->(x+y)/η
  end
  lck = ReentrantLock()
  @time Threads.@threads for _ in 1:300
    Pfc = fluctuating_coalescent(n; fuse)
    _A, _C = treevalues!(Pfc)
    lock(lck)
    try
      append!(Afc, _A)
      append!(Cfc, _C)
    finally
      unlock(lck)
    end
  end
  flush(stdout)
  dffc = to_mean_dataframe(Afc,Cfc)
  dffc.eta .= η
  dffc.n .= n
  append!(dffc_eta, dffc, cols=:union)
end)


begin
  plot_etas = [collect(1.2:0.005:1.4)[1:1:end]; [0.0,0.5,1.0,2.0]]
  plot_etas = 1.37:0.005:1.38
  plot_etas = [1.37,1.375,1.2, 1.0]
  plot_etas = levels(dffc_eta.eta)
  plot_etas = [1.0, 2.0]
  etac = 1.375
  α = 0.0#1/3
  β = 1#1.54 #1.515
  @transform!(dffc_eta, :coverabeta = :covera ./ :a.^(β-1))
  if α===-0.0
    @transform!(dffc_eta, :x = @. -1/log(abs((:eta-etac)/etac))*:a )
  else
    @transform!(dffc_eta, :x = abs.((:eta.-etac)/etac).^α .* :a )
  end
  # @transform!(dffc_eta, :x = :a )
  dfsub = @subset(dffc_eta, :eta .∈ Ref(plot_etas), :x.>=0, :x.<Inf)
  @transform!(dfsub, :x = abs.(:x))
  bw = 1/200
  bins = 10.0.^range(log10.(extrema(dfsub.x).+eps())..., step=bw)
  dfsub.xbin = cut(dfsub.x, bins[2:end-1], labels=midpoints(bins), extend=true)|>Array
  cdfsub = combine(groupby(dfsub, [:n, :eta, :xbin]), :coverabeta => mean => :coverabeta)
end
begin
  fig2 = Figure(resolution=72 .*(11,8), fontsize=12)
  ax3 = Axis(fig2[1,1], xscale=log2, yscale=log2, backgroundcolor=:lightgray)

  etas = levels(cdfsub.eta)
  etamin, etamax = extrema(etas)
  colors = cgrad(:RdYlGn_11)
  for d in groupby(cdfsub, [:n, :eta])
    eta = d[1, :eta]
    n = d[1, :n]
    #d = @subset(cdfsub, :eta.==eta)
    if eta == etac
      marker=:rect
    else
      marker=:circle
    end
    Makie.scatter!(ax3, d.xbin, d.coverabeta, markersize=8, label="$eta",
     #color=get(colors, (eta-etac)/(1.5) + 1/2), marker=marker
     )
  end
  xlims!(ax3, 1, 1e5)
  ylims!(ax3, 1e-1, 1e2)
  Legend(fig2[1,2], ax3, nbanks=2, bgcolor=:lightgray)
  # begin
  #   x1, y1 = 4, 14.2
  #   x0, y0 = 3, 10.5
  #   m = y1-y0
  #   lines!(ax3, 1e1..1e5, x->m*log10(x)+((y1+y0)-m*(x1+x0))/2, color=:black, linewidth=4,linestyle=:dash)
  # end
  # begin
  #   x1, y1 = 4, 52
  #   x0, y0 = 3, 25
  #   m = y1-y0
  #   lines!(ax3, 2e2..1e5, x->m*log10(x)+((y1+y0)-m*(x1+x0))/2, color=:black, linewidth=4,linestyle=:dash)
  # end
  # lines!(ax3, 1e0..1e3, x->x^(0.50), linewidth=4, linestyle=:dot)
  # lines!(ax3, 1e2..1e4, x->10^(-1/2)*x^(1.0), linewidth=4, linestyle=:dot)
  # hlines!(ax3, 10^-(0.04))
  lines!(ax3, 1e0..1e5, x->x^(0.515))
  # lines!(ax3, 1e2..1e5, x->0.9*(1e2/x)^(0.25))
  fig2
end 