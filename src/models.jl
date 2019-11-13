module models

using ..materials
using ..ft2d

export Custom,Layer,PlainLayer,PatternedLayer,Model,Circle,Rectangle,Ellipse,reciprocal
export Combination,Rotation,Shift
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
struct Custom <: Geometry
    F::AbstractArray{Complex{Float64},2}
end
function reciprocal(c::Custom,dnx,dny)
    return c.F
end
struct Rotation <: Geometry
    G::Geometry
    θ::Float64
end
function reciprocal(c::Rotation,dnx,dny)
    dnx2=dnx*cos(c.θ)+dny*sin(c.θ)
    dny2=-dnx*sin(c.θ)+dny*cos(c.θ)
    return reciprocal(c.G,dnx2,dny2)
end
struct Shift <: Geometry
    G::Geometry
    ax::Float64
    ay::Float64
end
function reciprocal(c::Shift,dnx,dny)
    return reciprocal(c.G,dnx,dny).*exp.(-2im*pi*c.ax.*dnx).*exp.(-2im*pi*c.ay.*dny)
end
struct Combination <: Geometry
    G::Array{Geometry,1}
end
function reciprocal(c::Combination,dnx,dny)
    ret=dnx*0
    for i=1:length(c.G)
        ret+=reciprocal(c.G[i],dnx,dny)
    end
    return ret
end

struct Ellipse <: Geometry
    rx::Float64
    ry::Float64
end
function reciprocal(e::Ellipse,dnx,dny)
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

#function plainLayer(thi,mat)
#    return PlainLayer(thi,mat)
#end

end
