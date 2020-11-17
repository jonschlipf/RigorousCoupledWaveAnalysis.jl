export Custom,Circle,Rectangle,Ellipse
using SpecialFunctions
"""
    Circle(d)
A circular inclusion
# Attributes
* `d` : diameter (relative to unit cell)
"""
struct Circle <: Geometry
    d::Float64
end
"""
    Rectangle(dx,dy)
A rectangular inclusion
# Attributes
* `dx` : width (relative to unit cell)
* `dy` : height (relative to unit cell)
"""
struct Rectangle <: Geometry
    dx::Float64
    dy::Float64
end
"""
    Custom(F)
A custom inclusion
# Attributes
* `F` : predefinded custom reciprocal space matrix
"""
struct Custom <: Geometry
    F::Array{Complex{Float64},2}
end
"""
    Ellipse(dx,dy)
An elliptic inclusion
# Attributes
* `dx` : x-axis diameter (relative to unit cell)
* `dy` : y-axis diameter (relative to unit cell)
"""
struct Ellipse <: Geometry
    dx::Float64
    dy::Float64
end
function reciprocal(c::Circle,dnx,dny)
	radix=.5*sqrt.((dnx*c.d).^2+(dny*c.d).^2)
    result=besselj.(1,2*pi*radix)./radix
    #asymptotic value for going towards 0
    result[radix.==0].=pi
    #scale amplitude
    return .25result*c.d^2
end
function reciprocal(r::Rectangle,dnx,dny)
	return r.dx*r.dy*sinc.(r.dx*dnx).*sinc.(r.dy*dny)
end
function reciprocal(c::Custom,dnx,dny)
    return c.F
end
function reciprocal(e::Ellipse,dnx,dny)
	radix=.5*sqrt.((dnx*e.dx).^2+(dny*e.dy).^2)
    result=besselj.(1,2*pi*radix)./radix
    #asymptotic value for going towards 0
    result[radix.==0].=pi
    #scale amplitude
    return .25result*e.dx*e.dy
end
function drawable(r::Rectangle)
    return [.5*r.dx,-.5*r.dx,-.5*r.dx,.5*r.dx,.5*r.dx],[.5*r.dy,.5*r.dy,-.5*r.dy,-.5*r.dy,.5*r.dy]
end
function drawable(c::Custom)
    #not yet implemented...
	return 0,0
end
function drawable(c::Circle)
    phivals=2pi*(0:.01:1)
    return c.d*.5cos.(phivals),c.d*.5sin.(phivals)
end

function drawable(c::Ellipse)
    phivals=2pi*(0:.01:1)
    return c.dx*.5cos.(phivals),c.dy*.5sin.(phivals)
end

