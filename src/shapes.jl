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
    return [.5*r.dx,-.5*r.dx,-.5*r.dx,.5*r.dx,.5*r.dx],[.5*r.dy,.5*r.dy,-.5*r.dy,-.5*r.dy,.5*r.dy]
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
