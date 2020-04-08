
using Interpolations

export Material,ConstantPerm,get_permittivity,InterpolPerm,InterpolPermA,ModelPerm,ConstantPermA

abstract type Material end
abstract type Isotropic <: Material end
abstract type Anisotropic <: Material end

struct ConstantPerm <: Isotropic
    ε::Complex{Float64}
end
function get_permittivity(mat::ConstantPerm,λ,index=1)
    return Complex(mat.ε)
end
struct ModelPerm<: Isotropic
    f::Function
end
function get_permittivity(mat::ModelPerm,λ,index=1)
    return Complex(mat.f(λ))
end
struct InterpolPerm <: Isotropic
    ε
end
function get_permittivity(mat::InterpolPerm,λ,index=1)
    return Complex(mat.ε(λ))
end
struct InterpolPermA<: Anisotropic
    ε
end
function get_permittivity(mat::InterpolPermA,λ,index=1)
    return Complex(mat.ε[index](λ))
end
struct ConstantPermA<: Anisotropic
    ε
end
function get_permittivity(mat::ConstantPermA,λ,index=1)
    return mat.ε[index]
end
