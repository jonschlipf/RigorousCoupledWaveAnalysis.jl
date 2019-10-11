module materials

export Material,ConstantPerm,get_permittivity

abstract type Material end

struct ConstantPerm <: Material
    ε::Complex{Float64}
end

function get_permittivity(mat::ConstantPerm,λ)
    return mat.ε
end

end
