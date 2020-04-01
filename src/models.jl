

module models

using ..ft2d
using ..materials




export Custom,Layer,PlainLayer,PatternedLayer,RCWAModel,Circle,Rectangle,Ellipse,reciprocal
export Combination,Rotation,Shift,drawable,Geometry

abstract type Geometry end

include("shapes.jl")
include("transforms.jl")




abstract type Layer end
struct PatternedLayer <: Layer
    thickness::Float64
    materials::Array{Material,1}
    geometries::Array{Geometry,1}
end
struct PlainLayer <: Layer
    thickness::Float64
    material::Material
end

struct RCWAModel
    layers::AbstractArray{Layer,1}
    εsup::Material
    εsub::Material
end

#function plainLayer(thi,mat)
#    return PlainLayer(thi,mat)
#end

end
