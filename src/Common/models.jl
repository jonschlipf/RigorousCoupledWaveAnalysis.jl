









export Custom,Layer,SimpleLayer,PatternedLayer,RCWAModel,Circle,Rectangle,Ellipse,reciprocal
export Combination,Rotation,Shift,drawable,Geometry,AnisotropicLayer

abstract type Geometry end

include("shapes.jl")
include("transforms.jl")




abstract type Layer end
struct PatternedLayer <: Layer
    thickness::Float64
    materials::Array{Material,1}
    geometries::Array{Geometry,1}
end
struct SimpleLayer <: Layer
    thickness::Float64
    material::Material
end
struct AnisotropicLayer <: Layer
    thickness::Float64
    material::Material
end

struct RCWAModel
    layers::AbstractArray{Layer,1}
    εsup::Material
    εsub::Material
end

#function SimpleLayer(thi,mat)
#    return SimpleLayer(thi,mat)
#end
