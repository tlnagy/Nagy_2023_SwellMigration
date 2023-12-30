using Documenter
using DemoCards

# Download necessary files from DataDryad
include(normpath(joinpath(@__DIR__, "..", "export", "dryad.jl")))
dst = normpath(joinpath(@__DIR__, "..", "data"))
files = ["fxm_uncaging_augmented.csv",
         "coulter_data.zip",
         "coulter_datasets.csv",
         "fxm_highres_datasets.csv",
         "fxm_highres_augmented_w_localities.csv",
         "fxm_highres_ds1_230310_V8_WT_fxmcorr.tif",
         "fxm_highres_ds2_230310_V8_BIX_fxmcorr.tif",
         "fxm_highres_ds3_231004_V8_Duvel_fxmcorr.tif"
]
for file in files
    if !isfile(joinpath(dst, file))
        @info "$file missing from disk, downloading..."
        download(root, versionpath, dst, file, headers)
        open(joinpath(dst, "dataset_version.txt"), "w") do io
            write(io, "version=$(version["versionNumber"])")
        end
    else
        @info "$file already present on disk"
    end
end

figpage, postprocess_cb, fig_assets = makedemos(joinpath("..", "figures"))
assets = []
isnothing(fig_assets) || (push!(assets, fig_assets))

makedocs(
    sitename="Neutrophils Actively Swell\nto Potentiate Migration",
    format = Documenter.HTML(
        size_threshold_ignore = ["figures/Notebooks/density_sim.md", "figures/Notebooks/latb_volume.md", "figures/Notebooks/motility.md", "figures/Notebooks/volume.md", "figures/Notebooks/singlecell.md"],
        assets =vcat([
                asset("https://analytics.tamasnagy.com/js/script.js", class=:js, attributes=Dict(Symbol("data-domain") => "tamasnagy.com", :defer => ""))
            ], assets),
    ),
    warnonly = true,
    pages = [
        "Home" => "index.md",
        figpage
    ]
)

postprocess_cb()

deploydocs(
    repo = "github.com/tlnagy/Nagy_2023_SwellMigration.git",
    versions=nothing
)