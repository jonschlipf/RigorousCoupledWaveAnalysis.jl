module materials
using Interpolations

<<<<<<< HEAD
export Material,ConstantPerm,get_permittivity,InterpolPerm,ModelPerm
=======
export Material,ConstantPerm,get_permittivity,InterpolPerm
>>>>>>> 00f14cd568c1e9c26d1fe34db6812f725b5ada19

abstract type Material end

struct ConstantPerm <: Material
    ε::Complex{Float64}
end

<<<<<<< HEAD
struct ModelPerm<: Material
    f::Function
end

=======
>>>>>>> 00f14cd568c1e9c26d1fe34db6812f725b5ada19
function get_permittivity(mat::ConstantPerm,λ)
    return mat.ε
end

struct InterpolPerm <: Material
    ε
end
function get_permittivity(mat::InterpolPerm,λ)
    return mat.ε(λ)
end
<<<<<<< HEAD
function get_permittivity(mat::ModelPerm,λ)
    return mat.f(λ)
end
=======

>>>>>>> 00f14cd568c1e9c26d1fe34db6812f725b5ada19
end
