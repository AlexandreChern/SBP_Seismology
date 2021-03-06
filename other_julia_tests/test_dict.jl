struct a
    a1::Int
    a2::Int
    a3::Dict{Int64,Array{Any,1}}
end

b = a(1,1,Dict(i=>[] for i=1:10))