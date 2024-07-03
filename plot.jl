using AlgebraOfGraphics, CairoMakie

CairoMakie.activate!(type="svg", pt_per_unit=1/0.75)

set_theme!(merge(Theme(linewidth=2,
  fontsize=12,
  color=:seaborn_deep,
  Scatter=(; alpha=0.5, marker=:rect, markersize=6),
  Legend=(; alpha=0.0)
), Makie.theme_latexfonts()))

begin
  _stopat = 2^19
  _σ = 3.0

  df_pc = @subset(dffc, :stopat .==_stopat, :run .== 1)
  p_pc = @subset(fc_params, :stopat .==_stopat, :run .== 1)

  df_nm = @subset(dfnm, :σ .==_σ, :run .==(1))
  p_nm = @subset(nm_params, :σ .==_σ, :run .==(1))


  fig = Figure(size=(7, 5) .* 72)
  colors = Makie.wong_colors()

  ax = Axis(fig[1, 1], xscale=log10, yscale=Makie.log10,
    xminorticksvisible=true, xminorticks=IntervalsBetween(9),
    yminorticksvisible=true, yminorticks=IntervalsBetween(9),
    xlabel=L"$A$", ylabel=L"$C/A$")
  ax2 = Axis(fig[1, 2], xscale=log10, yscale=Makie.log10,
    xminorticksvisible=true, xminorticks=IntervalsBetween(9),
    yminorticksvisible=true, yminorticks=IntervalsBetween(9),
    xlabel=L"$A$", ylabel=L"$D/A$")
  xlims!(ax, 1, 1e4)
  ylims!(ax, 2 / 3, 1e3)
  xlims!(ax2, 1, 1e4)
  ylims!(ax2, 1 / 2, 1e3)

  # p11 = Makie.lines!(ax, Aunb, Cunb./Aunb,
  #   label="Unbalanced", markersize=10, marker=:circle
  # )
  p16 = Makie.scatter!(ax, df_nm.A, df_nm.covera,
    label=L"Niche model ($σ = %$_σ$)"
  )
  p12 = Makie.scatter!(ax, df_pc.A, df_pc.covera,
    label="preferential coalescence"
  )

  p26 = Makie.scatter!(ax2, df_nm.A, df_nm.dovera,
    label=L"Niche model ($σ = 2.0$)"
  )
  p22 = Makie.scatter!(ax2, df_pc.A, df_pc.dovera,
    label="preferential coalescence"
  )

  params = reshape(fc_params[1, :cfit], :)
  p14 = Makie.lines!(ax, df_pc.A, power_law.(df_pc.A, Ref(params)),
    color=:gray37
  )
  axislegend(ax, [p14], [L"\sim A^{%$(round(params[2],digits=2))}"], position=:lt, framevisible=false)

  params = copy(params)
  x0 = log(10^2)
  beta = (params[2] -= 1/x0)
  params[1] += 1 - log(x0)
  push!(params, 0)
  params = reshape(p_pc[1, :dfit], :)
  p15 = Makie.lines!(ax2, df_pc.A, power_law.(df_pc.A, Ref(params)),
    color=:gray37
  )
  axislegend(ax2, [p15], [L"\sim A^{%$(round(params[2],digits=2))}"], position=:lt, framevisible=false)

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

  Legend(fig[2, 1:2], ax, nbanks=1, orientation=:horizontal, tellheight=true, position=:lt)

  current_figure()
end