module models

using ..materials
using ..ft2d

export Custom,Layer,PlainLayer,PatternedLayer,RCWAModel,Circle,Rectangle,Ellipse,reciprocal
export Combination,Rotation,Shift,drawable
abstract type Geometry end
struct Circle <: Geometry
    radius::Float64
end
function reciprocal(c::Circle,dnx,dny)
    return circft(c.radius,dnx,dny)
end
function drawable(c::Circle)
    phivals=2pi*(0:.01:1)
    return c.radius*.5cos.(phivals),c.radius*.5sin.(phivals)
end
struct Rectangle <: Geometry
    dx::Float64
    dy::Float64
end
function reciprocal(r::Rectangle,dnx,dny)
    return rectft(r.dx,r.dy,dnx,dny)
end
function drawable(r::Rectangle)
    return [.5*r.dx,-.5*r.dx,-.5*r.dx,.5*r.dx],[.5*r.dy,.5*r.dy,-.5*r.dy,-.5*r.dy]
end
struct Custom <: Geometry
    F::Array{Complex{Float64},2}
end
function reciprocal(c::Custom,dnx,dny)
    return c.F
end
function drawable(c::Custom)
    return 0,0
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
function drawable(c::Rotation)
    x,y=drawable(c.G)
    u=x*cos(-c.θ)+y*sin(-c.θ)
    v=-x*sin(-c.θ)+y*cos(-c.θ)
    return u,v
end
struct Shift <: Geometry
    G::Geometry
    ax::Float64
    ay::Float64
end
function reciprocal(c::Shift,dnx,dny)
    return reciprocal(c.G,dnx,dny).*exp.(-2im*pi*c.ax.*dnx).*exp.(-2im*pi*c.ay.*dny)
end
function drawable(c::Shift)
    x,y=drawable(c.G)
    return x.+c.ax,y.+c.ay
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
function drawable(c::Combination)
    x,y=drawable(c.G[1])
    for i=2:length(c.G)
        xon,yon=drawable(c.G[i])
        x=cat(x,xon,dims=1)
        y=cat(y,yon,dims=1)
    end
    return x,y
end
struct Ellipse <: Geometry
    rx::Float64
    ry::Float64
end
function reciprocal(e::Ellipse,dnx,dny)
    return ellipft(e.rx,e.ry,dnx,dny)
end

function drawable(c::Ellipse)
    phivals=2pi*(0:.01:1)
    return c.rx*.5cos.(phivals),c.ry*.5sin.(phivals)
end

struct arbitrary <: Geometry
    mask::Array{Float64,2}
end


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
