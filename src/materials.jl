module materials
using Interpolations

export Material,ConstantPerm,get_permittivity,InterpolPerm

abstract type Material end

struct ConstantPerm <: Material
    ε::Complex{Float64}
end

function get_permittivity(mat::ConstantPerm,λ)
    return mat.ε
end

struct InterpolPerm <: Material
    ε
end
function get_permittivity(mat::InterpolPerm,λ)
    return mat.ε(λ)
end

end
