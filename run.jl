using Base.Iterators: product
using BinaryTrees, TreeProcesses
using CategoricalArrays, DataFrames, DataFramesMeta, StatsBase
using JLD2
using EasyFit


n = 2^19 #4_000_000

## Aggregate A,C,D over multiple weighted_coalescent runs.
# fuse = (a,b)->max(0, (a+b)/1)
stopats = [n÷2]
αs = [1.0]
innerruns = 10
outerruns = 1

# Afc = Int[]
# Cfc = Int[]
# Dfc = Int[]
dffc = DataFrame()
@time foreach(product(1:outerruns, αs)) do (run, α)
  fuse = (a,b)->(a+b)*α
  # empty!(Afc)
  # empty!(Cfc)
  # empty!(Dfc)
  @time Pfc = map(1:innerruns) do iter
    P = preferential_coalescent(n, ones(n); fuse=fuse, stopat=1)[1]
    ACD!(P)
    @show iter
    P
  end
  for stopat in stopats
    Pslice = mapreduce(vcat, Pfc) do P
      TreeProcesses.time_slice(P, n-stopat)
    end
    obs = TreeProcesses.get_observables(Pslice)
    df = @views to_mean_dataframe(obs[:, 1], obs[:, 2], obs[:, 3])
    sort!(df, :A)
    df.run .= run
    df.stopat .= stopat
    df.α .= α
    select!(df, :stopat, :run, Not([:run, :stopat]))
    append!(dffc, df)
  end
end

## Full tree
Afc = Int[]
Cfc = Int[]
Dfc = Int[]
dffc2 = DataFrame()
@time foreach(product(1:outerruns, αs)) do (run, α)
  fuse = (a,b)->(a+b)*α
  empty!(Afc)
  empty!(Cfc)
  empty!(Dfc)
  @time foreach(1:innerruns) do iter
    P = preferential_coalescent(n, ones(n); fuse=fuse, stopat=1)[1]
    _,a,c,d = ACD!(P)
    append!(Afc, a)
    append!(Cfc, c)
    append!(Dfc, d)
  end
  df = @views to_mean_dataframe(Afc, Cfc, Dfc)
  sort!(df, :A)
  df.run .= run
  df.α .= α
  df.stopat .= 1
  select!(df, :run, Not([:run]))
  append!(dffc2, df)
end

fc_params = combine(groupby(dffc, [:stopat, :run, :α])) do dffc
  df = @subset dffc 1e0 .<= :A .<= 2e2
  fcfit_c_pow = fitlinear(log.(df.A), log.(df.covera))

  df = @subset dffc min(maximum(:A), 1e2) .<= :A .<= 2e2
  fcfit_d_pow = fitlinear(log.(df.A), log.(df.dovera))

  (cfit=[exp(fcfit_c_pow.b) fcfit_c_pow.a], dfit=[exp(fcfit_d_pow.b) fcfit_d_pow.a])
end

fc_params2 = combine(groupby(dffc2, [:stopat, :run, :α])) do dffc
  df = @subset dffc 1e0 .<= :A .<= 2e2
  fcfit_c_pow = fitlinear(log.(df.A), log.(df.covera))

  df = @subset dffc min(maximum(:A), 1e2) .<= :A .<= 2e2
  fcfit_d_pow = fitlinear(log.(df.A), log.(df.dovera))

  (cfit=[exp(fcfit_c_pow.b) fcfit_c_pow.a], dfit=[exp(fcfit_d_pow.b) fcfit_d_pow.a])
end

## Aggregate A,C over multiple niche model runs.
n_niche = n÷2
maxiters = 1_000
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
    # @info σ
    while iter < maxiters && k < n_niche
      Pnm, _, k = nichemodel(n_niche, σ; ϵ=0.0, excludepruned=true)
      iter += 1
    end
    if iter == maxiters
      @warn("Tree could not be generated")
    else
      @info "Tree found after $iter iterations"
    end

    _,a,c,d = ACD!(Pnm)
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
