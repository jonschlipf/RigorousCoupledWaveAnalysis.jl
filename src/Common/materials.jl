
using Interpolations

export Material,ConstantPerm,get_permittivity,InterpolPerm,InterpolPermA,ModelPerm,ConstantPermA,ModelPermA

abstract type Material end
abstract type Isotropic <: Material end
abstract type Anisotropic <: Material end

"""
    ConstantPerm(ε)
RCWA isotropic material model with wavelength-independent permittivity
# Attributes
* `ε` : Complex scalar denoting the permittivity
"""
struct ConstantPerm <: Isotropic
    ε::Complex{Float64}
end
"""
    ConstantPermA(ε)
RCWA anisotropic material model with wavelength-independent permittivity
# Attributes
* `ε` : Complex 5-element vector ([εxx,εxy,εyx,εyy,εzz]) denoting the permittivity
"""
struct ConstantPermA<: Anisotropic
    ε
end
"""
    ModelPerm(ε)
RCWA isotropic material model with the permittivity defined by an analytical model
# Attributes
* `f` : Function to compute the complex permittivity value for a single wavelength
"""
struct ModelPerm<: Isotropic
    f::Function
end
"""
    ModelPermA(ε)
RCWA anisotropic material model with the permittivity defined by an analytical model
# Attributes
* `f` : Vector containing five functions to compute the complex permittivity tensor for a single wavelength ([εxx,εxy,εyx,εyy,εzz])
"""
struct ModelPermA<: Anisotropic
    f::Vector{Function}
end
"""
    InterpolPerm(ε)
RCWA isotropic material model with permittivity given by interpolated sampled data
# Attributes
* `ε` : Interpolation object (see the Interpolations.jl package)
"""
struct InterpolPerm <: Isotropic
    ε
end
"""
    InterpolPermA(ε)
RCWA anisotropic material model with wavelength-independent permittivity
# Attributes
* `ε` : Interpolation object (see the Interpolations.jl package), interpolates 5-element vectors ([εxx,εxy,εyx,εyy,εzz])
"""
struct InterpolPermA<: Anisotropic
    ε
end
"""
    get_permittivity(mat::Material,λ)
    get_permittivity(mat::Material,λ,index)
Returns the permittivity of a material model object for given wavelength and coordinate direction
# Arguments
* `mat` :  material object
* `λ` :  free space wavelength
* `index` :  optional, specifies the coordinate direction for anisotropic materials
# Outputs
* `ε` : Complex scalar permittivity value for the given wavelength and coordinate direction
"""
function get_permittivity(mat::ConstantPerm,λ,index=1)
    return Complex(mat.ε)*(index==1||index==4||index==5)
end
function get_permittivity(mat::ModelPerm,λ,index=1)
    return Complex(mat.f(λ))*(index==1||index==4||index==5)
end
function get_permittivity(mat::ModelPermA,λ,index=1)
    return Complex(mat.f[index](λ))
end
function get_permittivity(mat::InterpolPerm,λ,index=1)
    return Complex(mat.ε(λ))*(index==1||index==4||index==5)
end
function get_permittivity(mat::InterpolPermA,λ,index=1)
    return Complex(mat.ε[index](λ))
end
function get_permittivity(mat::ConstantPermA,λ,index=1)
    return mat.ε[index]
end
