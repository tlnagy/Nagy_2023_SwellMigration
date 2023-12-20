using Documenter
using DemoCards

# Download necessary files from DataDryad
include(normpath(joinpath(@__DIR__, "..", "export", "dryad.jl")))
dst = normpath(joinpath(@__DIR__, "..", "data"))
files = ["fxm_uncaging_augmented.csv", "coulter_data.zip", "coulter_datasets.csv"]
for file in files
    if !isfile(joinpath(dst, file))
        @info "$file missing from disk, downloading..."
        download(root, versionpath, dst, file, headers)
    else
        @info "$file already present on disk"
    end
end

figpage, postprocess_cb, fig_assets = makedemos(joinpath("..", "figures"))
assets = []
isnothing(fig_assets) || (push!(assets, fig_assets))

makedocs(
    sitename="Neutrophils Actively Swell to Potentiate Migration",
    format = Documenter.HTML(
        size_threshold_ignore = ["figures/Notebooks/density_sim.md", "figures/Notebooks/latb_volume.md", "figures/Notebooks/motility.md"],
        assets = assets
    ),
    warnonly = true,
    pages = [
        "Home" => "index.md",
        figpage
    ]
)

postprocess_cb()

# deploydocs(
#     repo = "github.com/tlnagy/Nagy_2023_SwellMigration.git",
# )