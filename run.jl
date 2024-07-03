using Base.Iterators: product
using BinaryTrees, TreeProcesses
using CategoricalArrays, DataFrames, DataFramesMeta, StatsBase
using JLD2
using EasyFit


power_log_law(x,p) = @. p[1]*x^p[2]*log(x)
power_law(x,p) = @. p[1]*x^p[2]
power_law_off(x,p) = @. p[1]*x^p[2] + p[3]

n = 2^20 #4_000_000

## Aggregate A,C,D over multiple weighted_coalescent runs.
fuse = (a,b)->a+b
stopats = [n÷2]
innerruns = 10
outerruns = 1

Afc = Int[]
Cfc = Int[]
Dfc = Int[]
dffc = DataFrame()
foreach(product(stopats, 1:outerruns)) do (stopat, run)
  empty!(Afc)
  empty!(Cfc)
  empty!(Dfc)
  @time for iter in 1:innerruns
    Pfc = preferential_coalescent(n, ones(n); fuse=fuse, stopat=stopat)
    foreach((ACD!(P) for P in Pfc)) do (k,a,c,d)
      push!(Afc,last(a))
      push!(Cfc,last(c))
      push!(Dfc,last(d))
    end
    @show iter
  end
  df = to_mean_dataframe(Afc, Cfc, Dfc)
  sort!(df, :A)
  df.run .= run
  df.stopat .= stopat
  select!(df, :stopat, :run, Not([:run, :stopat]))
  append!(dffc, df)
end

fc_params = combine(groupby(dffc, [:stopat, :run])) do dffc
  df = @subset dffc 1e0 .<= :A .<= 2e2
  fcfit_c_pow = fitlinear(log.(df.A), log.(df.covera))

  df = @subset dffc 5e1 .<= :A .<= 2e2
  fcfit_d_pow = fitlinear(log.(df.A), log.(df.dovera))

  (cfit=[exp(fcfit_c_pow.b) fcfit_c_pow.a], dfit=[exp(fcfit_d_pow.b) fcfit_d_pow.a])
end

## Aggregate A,C over multiple niche model runs.
n_niche = 2^20
maxiters = 1000
sigmas = [2.0, 3.0]
innerruns = 10
outerruns = 1

Anm = Int[]
Cnm = Int[]
Dnm = Int[]
TS = Int[]
dfnm = DataFrame()
foreach(product(sigmas, 1:outerruns)) do (σ, run)
  empty!(Anm)
  empty!(Cnm)
  empty!(Dnm)
  @time for _ in 1:innerruns
    iter = 0
    k = 0
    Pnm = nothing
    while iter < maxiters && k < n_niche
      Pnm, k = nichemodel(n_niche, σ; ϵ=0.0)
      iter += 1
    end
    iter == maxiters && @warn("Tree could not be generated") || @info "Tree found after $iter iterations"

    k,a,c,d = ACD!(Pnm)
    append!(Anm,a)
    append!(Cnm,c)
    append!(Dnm,d)

  end
  df = to_mean_dataframe(Anm, Cnm, Dnm)
  sort!(df, :A)
  df.run .= run
  df.σ .= σ
  select!(df, :σ, :run, Not([:run, :σ]))
  append!(dfnm, df)
end

nm_params = combine(groupby(dfnm, [:σ, :run])) do dffc
  df = @subset dffc 1e0 .<= :A .<= 2e2
  fcfit_c_pow = fitlinear(log.(df.A), log.(df.covera))
  df = @subset dffc 1e1 .<= :A .<= 2e2
  fcfit_d_pow = fitlinear(log.(df.A), log.(df.dovera))

  (cfit=[exp(fcfit_c_pow.b) fcfit_c_pow.a], dfit=[exp(fcfit_d_pow.b) fcfit_d_pow.a])
end


## Comparison: unbalanced tree C~A^2
Punbalanced  = maximally_unbalanced(n÷2)
_, Aunb, Cunb, Dunb = ACD!(Punbalanced)
dfunb = to_mean_dataframe(Aunb, Cunb, Dunb)

## Comparison: balanced tree C~A*log(A)
Pbalanced  = maximally_balanced(log2(n))
_, Ab, Cb, Db = ACD!(Pbalanced)
dfb = to_mean_dataframe(Ab, Cb, Db);
nothing
