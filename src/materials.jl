module materials

export Material,Constant,get_permittivity

abstract type Material end

struct Constant <: Material
    ε::Complex{Float64}
end

function get_permittivity(mat::Constant,λ)
    return mat.ε
end

end
