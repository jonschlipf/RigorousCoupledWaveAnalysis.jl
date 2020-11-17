export Layer,SimpleLayer,PatternedLayer,AnisotropicLayer,RCWAModel
export Geometry,drawable,reciprocal

"""
    Geometry
Structure to store complex geometries for layer patterning
"""
abstract type Geometry end
#Geometries are specified here:
include("shapes.jl")
include("transforms.jl")
"""
    Layer
Structure to store the layers of the RCWA model
"""
abstract type Layer end
"""
    PatternedLayer(t,materials,geometries)

An inhomogenous layer with two or more materials forming a pattern
# Attributes
* `t` : Layer thickness
* `materials` : n-element array for the n materials making up the layer
* `geometries` : (n-1)-element array specifying the patterning
"""
struct PatternedLayer <: Layer
    thickness::Float64
    materials::Array{Material,1}
    geometries::Array{Geometry,1}
end
"""
    SimpleLayer(t,material)

A homogenous, isotropic layer
# Attributes
* `t` : Layer thickness
* `material` : material model for the permittivity of the layer
"""
struct SimpleLayer <: Layer
    thickness::Float64
    material::Material
end
"""
    AnisotropicLayer(t,material)

A homogenous, in-plane anisotropic layer
# Attributes
* `t` : Layer thickness
* `material` : material model for the permittivity of the layer
"""
struct AnisotropicLayer <: Layer
    thickness::Float64
    material::Material
end
"""
    RCWAModel(layers,εsup,εsub)

A structure to store the full RCWA stack including superstrate and substrate
# Attributes
* `layers` : Array of layer objects
* `εsup` : material model for the superstrate
* `εsub` : material model for the substrate
"""
struct RCWAModel
    layers::Array{Layer,1}
    εsup::Material
    εsub::Material
end

