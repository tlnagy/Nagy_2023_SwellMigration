using CSV
using DataFrames
using TiffImages
using OMETIFF
using FileIO
using AxisArrays
using ImageCore
using ImageContrastAdjustment
using ProgressMeter
using Measures

# datasets = CSV.File(normpath(@__DIR__, "../data/fxm/fxm_dataset_info.csv")) |> DataFrame
datasets = CSV.File(normpath(@__DIR__, "../data/fxm/fxm_highres_datasets.csv")) |> DataFrame
# filter!(x->x.used == "Y", datasets)

# ds = datasets[1, :]
function prepdata(ds)
    newdspath = normpath(@__DIR__, joinpath("..", "data", "fxm", "highres", ds.new_dataset))
    if isdir(newdspath)
        println("Dataset $(ds.new_dataset) already exists, skipping")
        return
    end

    rawomedir = normpath(@__DIR__, joinpath(ds.old_dataset, "data_raw"))
    rawomepath = filter(x->endswith(x, ".ome.tif"), readdir(rawomedir)) |> first

    rawome = load(joinpath(rawomedir, rawomepath))

    if :position in axisnames(rawome)
        rawimg = view(rawome, Axis{:position}(Symbol(ds.lane)))
    else
        rawimg = rawome
    end
    # fix ordering
    chandim = axisdim(rawimg, Axis{:channel}) 
    (chandim !== 3) && @warn "Weird channel config"
    channames = axisvalues(rawimg)[chandim]

    fchan = y ->findfirst(x->occursin(y, string(x)), channames)

    rawimgord = rawimg[:, :, fchan.(["BFP", "RFP", "DeepRed"]), :]
    rawimgclr = colorview(RGB, PermutedDimsArray(rawimgord, (3, 1, 2, 4)))

    newdata = mkpath(newdspath)
    outimg = empty(LazyBufferedTIFF, RGB{N0f16}, joinpath(newdata, "raw_img.tif"))
    @showprogress for slice in eachslice(rawimgclr, dims = 3)
        push!(outimg, slice)
    end
    close(outimg)

    procdir = normpath(@__DIR__, ds.old_dataset)

    fxmcorr = TiffImages.load(joinpath(procdir, "fxmcorr.tif"))

    TiffImages.save(joinpath(newdata, "fxmcorr.tif"), convert.(Gray{N0f16}, adjust_histogram(fxmcorr, LinearStretching(nothing => (0, 1)))))

    # cp(joinpath(procdir, "augmented.csv"), joinpath(newdata, "augmented.csv"))
    # footprints = CSV.File(joinpath(procdir, "augmented_w_footprints.csv")) |> DataFrame
    # footprints |> CSV.write(joinpath(newdata, "augmented_w_footprints_gzipped.csv"), compress = true)
    # try
    #     localities = CSV.File(joinpath(procdir, "augmented_w_localities.csv")) |> DataFrame
    #     localities |> CSV.write(joinpath(newdata, "augmented_w_localities_gzipped.csv"), compress = true)
    # catch
    #     println("No localities for $(ds.new_dataset)")
    # end
end

# ds = datasets[13, :]
for ds in eachrow(datasets)
    prepdata(ds)
    GC.gc()
end

include("../figures/utils.jl")

linked = vcat(map(eachrow(datasets)) do ds
    procdir = normpath(@__DIR__, ds.old_dataset)
    types = Dict([:rel_volume, :abs_volume_um3] .=> Measurement{Float64})
    types[:footprint_cart] = EmbedVector{CartesianIndex{2}, ';'}
    types[:locality_cart] = EmbedVector{CartesianIndex{2}, ';'}

    footprints = CSV.File(joinpath(procdir, "augmented_w_localities.csv"), types = types) |> DataFrame

    rename!(footprints, [s => join(split(s, "_")[1:end-1], "_") for s in filter(x->occursin("_387LP", x), names(footprints))])

    insertcols!(footprints, 1, :date => ds.date, 
                               :volunteer => ds.volunteer, 
                               :condition => ds.condition)
    insertcols!(footprints,
                               :lane => ds.lane,
                               :dish => ds.dish,
                               :scriptver => ds.scriptver,
                               :dataset => joinpath("movies", ds.new_dataset))
end...)
linked |> CSV.write(joinpath("data", "fxm", "fxm_highres_augmented_w_localities.csv"), compressed = false, missingstring="n/a")

types = Dict([:rel_volume, :abs_volume_um3] .=> Measurement{Float64})
types[:footprint_cart] = EmbedVector{CartesianIndex{2}, ';'}

footprints = CSV.File(joinpath("data", "fxm", "augmented_w_footprints.csv"), types = types) |> DataFrame

sort!(datasets, :date)
datasets.prefix .= ""
@showprogress for (i, ds) in enumerate(eachrow(datasets))
    src = joinpath(normpath(@__DIR__, "..", "data", "fxm", "fxm_uncaging_movies", ds.new_dataset))
    dst = joinpath(normpath(@__DIR__, "..", "data", "fxm", "flat"))
    prefix = "fxm_uncaging_ds$(i)_$(ds.new_dataset)"
    cp(joinpath(src, "fxmcorr.tif"), joinpath(dst, "$(prefix)_fxmcorr.tif"))
    cp(joinpath(src, "raw_img.tif"), joinpath(dst, "$(prefix)_rawrgb.tif"))
    ds.prefix = prefix
end
