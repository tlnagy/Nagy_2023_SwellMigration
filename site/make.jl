using Documenter
using DemoCards

figpage, postprocess_cb, fig_assets = makedemos(joinpath("..", "figures"))
assets = []
isnothing(fig_assets) || (push!(assets, fig_assets))

makedocs(
    sitename="Neutrophils Actively Swell to Potentiate Migration",
    format = Documenter.HTML(
        size_threshold_ignore = ["figures/Notebooks/density_sim.md"],
        assets = assets
    ),
    warnonly = true,
    pages = [
        "Home" => "index.md",
        figpage
    ]
)

postprocess_cb()