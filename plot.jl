using AlgebraOfGraphics, CairoMakie

power_log_law(x,p) = @. p[1]*x^p[2]*log(x)
power_law(x,p) = @. p[1]*x^p[2]
power_law_off(x,p) = @. p[1]*x^p[2] + p[3]

include("plotheme.jl")

function setmarkersize!(leg, ms)
  for i in eachindex(leg.entrygroups[][1][2])
    leg.entrygroups[][1][2][i].elements[1].attributes[:markersize] = Observable(ms)
  end
  notify(leg.entrygroups)
end

begin
_stopat = 1#2^18
_σ = 3.0

# df_pc = @subset(dffc, :stopat .==_stopat, :run .== 1, :α.==1.0)
df_pc = @subset(dffc2, :A .<= 1e4, :run .== 1, :α.==1.0)
p_pc = @subset(fc_params2, :stopat .==_stopat, :run .== 1, :α.==1.0)

df_nm = @subset(dfnm, :A .< 1e4, :σ .==_σ, :run .==(1))
p_nm = @subset(nm_params, :σ .==_σ, :run .==(1))
end

begin
  fig = Figure(size=(510,510), fontsize=20)
  colors = Makie.wong_colors()

  ax = Axis(fig[1, 1], xscale=log10, yscale=Makie.log10,
    xminorticksvisible=true, xminorticks=IntervalsBetween(9),
    yminorticksvisible=true, yminorticks=IntervalsBetween(9),
    xlabel=L"$A$", ylabel=L"$C/A$")
  xlims!(ax, 1, 1e4)
  ylims!(ax, 2 / 3, 1e2)

  # p11 = Makie.lines!(ax, Aunb, Cunb./Aunb,
  #   label="Unbalanced", markersize=10, marker=:circle
  # )
  params = reshape(p_pc[1, :cfit], :)
  # p14 = Makie.lines!(ax, 1..2e2, x->power_law(x, (params)),
  #   color=:gray37
  # )
  p12 = Makie.scatter!(ax, df_pc.A, df_pc.covera,
  label="preferential coalescence", markersize=6
  )
  p16 = Makie.scatter!(ax, df_nm.A, df_nm.covera,
    label=L"Niche model ($σ = %$_σ$)",
    alpha=0.1, markersize=6
  )
  
  lines!(ax, 1..10^4, x->x^(1/2), color=:red)

  # p14 = Makie.lines!(ax, df_pc.A, power_law.(df_pc.A, Ref(params)),
  #   color=:gray37, linestyle=:solid, label=L"\sim A^{%$(round(params[2],digits=2))}"
  # )

  # leg = axislegend(ax, position=:lt, framevisible=false)
  # setmarkersize!(leg, 6)


  # params = copy(params)
  # x0 = log(10^2)
  # beta = (params[2] -= 1/x0)
  # params[1] += 1 - log(x0)
  # push!(params, 0)

  p21 = Makie.lines!(ax, dfb.A, dfb.covera,
  linestyle=:dash,
  linewidth=1,
  color=colors[3]
  )
  p22 = Makie.lines!(ax, dfunb.A, dfunb.covera,
  linestyle=:dash,
  linewidth=1,
  color=colors[3]
  )
  save("covera_vs_a.pdf", fig)
  fig
end

begin
  fig = Figure(size=(230,230), fontsize=20)
  colors = Makie.wong_colors()

  ax = Axis(fig[1, 1], xscale=log10, yscale=Makie.log10,
    xminorticksvisible=true, xminorticks=IntervalsBetween(9),
    yminorticksvisible=true, yminorticks=IntervalsBetween(9),
    xlabel=L"$A$", ylabel=L"$C/A$")
  xlims!(ax, 1, 1e4)
  ylims!(ax, 2 / 3, 1e2)

  # p11 = Makie.lines!(ax, Aunb, Cunb./Aunb,
  #   label="Unbalanced", markersize=10, marker=:circle
  # )
  params = reshape(p_pc[1, :cfit], :)
  # p14 = Makie.lines!(ax, 1..2e2, x->power_law(x, (params)),
  #   color=:gray37
  # )
  p12 = Makie.scatter!(ax, df_pc.A, df_pc.covera,
  label="preferential coalescence"
  )
  p16 = Makie.scatter!(ax, df_nm.A, df_nm.covera,
    label=L"Niche model ($σ = %$_σ$)",
    alpha=0.1
  )
  
  # p14 = Makie.lines!(ax, df_pc.A, power_law.(df_pc.A, Ref(params)),
  #   color=:gray37, linestyle=:solid, label=L"\sim A^{%$(round(params[2],digits=2))}"
  # )

  # leg = axislegend(ax, position=:lt, framevisible=false)
  setmarkersize!(leg, 6)


  # params = copy(params)
  # x0 = log(10^2)
  # beta = (params[2] -= 1/x0)
  # params[1] += 1 - log(x0)
  # push!(params, 0)

  p21 = Makie.lines!(ax, dfb.A, dfb.covera,
  linestyle=:dash,
  linewidth=1,
  color=colors[3]
  )
  p22 = Makie.lines!(ax, dfunb.A, dfunb.covera,
  linestyle=:dash,
  linewidth=1,
  color=colors[3]
  )
  save("CoverA.pdf", fig)
  fig
end

begin
  fig2 = Figure(size=(5, 5) .* 72)
  colors = Makie.wong_colors()

  ax2 = Axis(fig2[1, 1], xscale=log10, yscale=Makie.log10,
  xminorticksvisible=true, xminorticks=IntervalsBetween(9),
  yminorticksvisible=true, yminorticks=IntervalsBetween(9),
  xlabel=L"$A$", ylabel=L"$D/A$")
  xlims!(ax2, 1, 1e4)
  ylims!(ax2, 1 / 2, 1e3)

  params = reshape(p_pc[1, :dfit], :)
  p15 = Makie.lines!(ax2, df_pc.A, power_law.(df_pc.A, Ref(params)),
    color=:gray37,
    label=L"\sim A^{%$(round(params[2],digits=2))}"
  )
  p22 = Makie.scatter!(ax2, df_pc.A, df_pc.dovera,
    label="preferential coalescence"
  )
  p26 = Makie.scatter!(ax2, df_nm.A, df_nm.dovera,
    label=L"Niche model ($σ = %$_σ$)",
    alpha=0.1
  )
  
  p21 = Makie.lines!(ax2, dfb.A, dfb.dovera,
  linestyle=:dash,
  linewidth=1,
  color=colors[3]
  )
  p22 = Makie.lines!(ax2, dfunb.A, dfunb.dovera,
  linestyle=:dash,
  linewidth=1,
  color=colors[3]
  )
  leg2 = axislegend(ax2, position=:lt, framevisible=false)
  setmarkersize!(leg2, 6)

  # Legend(fig[2, 1:2], ax, nbanks=1, orientation=:horizontal, tellheight=true, position=:lt)
  save("DoverA.pdf", fig2)
  fig2
end

begin
  p = 0.0
  p2 = 0.5
  fig3 = Figure()
  ax3 = Axis(fig3[1,1], xscale=log10, yscale=log10, xlabel="A", ylabel=L"f_A",
    xminorticksvisible=true, xminorticks=IntervalsBetween(9),
    yminorticksvisible=true, yminorticks=IntervalsBetween(9)
  )
  scatter!(df_pc.A, df_pc.nrow./sum(df_pc.nrow), label="preferential coalescence")
  scatter!(df_nm.A, df_nm.nrow./sum(df_nm.nrow), label="niche model")
  lines!(1..1e4, x->2*x^(-2.0))
  lines!(1..1e4, x->1*x^(-1.6))
  xlims!(ax3, 1, 1e3)
  leg3 = axislegend(ax3)
  setmarkersize!(leg3, 6)

  save("AvsN.pdf", fig3)
  fig3
end