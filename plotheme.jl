CairoMakie.activate!(type="svg", pt_per_unit=1/0.75)

my_default_thm = merge(Theme(linewidth=2,
fontsize=12,
color=:seaborn_deep,
Scatter=(; alpha=0.5, marker=:rect, markersize=3),
Legend=(; alpha=0.0)
), Makie.theme_latexfonts())

set_theme!(my_default_thm)