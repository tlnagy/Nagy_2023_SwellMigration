# https://github.com/JuliaData/CSV.jl/issues/740
struct EmbedVector{T, delim}
    values::Vector{T}
    
    EmbedVector{delim}(values::Vector{T}) where {T, delim} = new{T, delim}(values)
end

Base.zero(::Type{EmbedVector{T, delim}}) where {T, delim} = EmbedVector{delim}(T[]) 
Base.string(x::EmbedVector{T, delim}) where {T, delim} = join(map(x->x.I, x.values), delim)

Base.tryparse(::Type{EmbedVector{T, delim}}, str) where {T, delim} = EmbedVector{delim}(map(x -> tryparse(T, x), split(str, delim)))

function Base.tryparse(::Type{CartesianIndex{2}}, str)
    i = findfirst(',', str)
    if length(str) > 16
        CartesianIndex(parse(Int, str[16:i-1]), parse(Int, str[i+1:end-1]))
    else
        CartesianIndex(parse(Int, str[2:i-1]), parse(Int, str[i+1:end-1]))
    end
end

Base.length(x::EmbedVector) = length(x.values)

# https://mikemol.github.io/technique/colorblind/2018/02/11/color-safe-palette.html
okabe_ito = parse.(Colorant, ["#56B4E9", "#E69F00", "#CC79A7", "#009E73", "#000000", "#F0E442", "#0072B2", "#D55E00"])
default(
    guidefontsize=10, tickfontsize=10, titlefontsize=12,
    gridlinewidth=2, minorgridlinewidth=2, gridstyle=:dash, right_margin = 10mm,
    thickness_scaling=1.5, framestyle = :box, linewidth=1, tick_orientation = :out,
    fontfamily = "helvetica", palette = okabe_ito)