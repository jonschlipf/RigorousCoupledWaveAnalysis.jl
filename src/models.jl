module models

using ..materials
using ..ft2d

export Layer,PlainLayer,plainLayer,PatternedLayer,Model,Circle,Rectangle,Ellipse,reciprocal

abstract type Geometry end
struct Circle <: Geometry
    radius::Float64
end
function reciprocal(c::Circle,dnx,dny)
    return circft(c.radius,dnx,dny)
end
struct Rectangle <: Geometry
    dx::Float64
    dy::Float64
end
function reciprocal(r::Rectangle,dnx,dny)
    return rectft(r.dx,r.dy,dnx,dny)
end
struct Ellipse <: Geometry
    rx::Float64
    ry::Float64
end
function reciprocal(c::Ellipse,dnx,dny)
    return ellipft(e.rx,e.ry,dnx,dny)
end
struct arbitrary <: Geometry
    mask::Array{Float64,2}
end


abstract type Layer end
struct PatternedLayer <: Layer
    thickness::Float64
    materials::AbstractArray{Material,1}
    geometries::AbstractArray{Geometry,1}
end
struct PlainLayer <: Layer
    thickness::Float64
    material::Material
end

struct Model
    layers::AbstractArray{Layer,1}
    εsup::Material
    εsub::Material
end

function plainLayer(thi,mat)
    return PlainLayer(thi,mat)
end

end
