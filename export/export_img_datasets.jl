include("../../loader.jl")

# df = loadfxmcsv("augmented.csv")
datasets = CSV.File("data/fxm/fxm_highres_datasets.csv") |> DataFrame

meta_columns = [:date, :condition, :lane, :dish, :scriptver, :old_datapath]

datasets = combine(groupby(df, meta_columns),
    meta_columns .=> first)

datasets = CSV.File("image_datasets.csv") |> DataFrame
datasets.new_dataset .= @. Dates.format(datasets.date, "yymmdd") * "_" * datasets.volunteer * "_" * datasets.condition
CSV.write("data/fxm/fxm_highres_datasets.csv", datasets)