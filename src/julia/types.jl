abstract Protocells

type Protocell <: Protocells
    composition :: Vector{Int}
    devtime :: Float64
    metabolism::Function
end

type Population <: Protocells
    individuals::Array{Protocell}
end

