
using Interpolations

export Material,ConstantPerm,get_permittivity,InterpolPerm,ModelPerm

abstract type Material end

struct ConstantPerm <: Material
    ε::Complex{Float64}
end

struct ModelPerm<: Material
    f::Function
end

function get_permittivity(mat::ConstantPerm,λ)
    return Complex(mat.ε)
end

struct InterpolPerm <: Material
    ε
end
function get_permittivity(mat::InterpolPerm,λ)
    return Complex(mat.ε(λ))
end
function get_permittivity(mat::ModelPerm,λ)
    return Complex(mat.f(λ))
end
