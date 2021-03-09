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
#Already defined in transforms.jl
#"""
#    reciprocal(geo::Geometry,dnx,dny)
#Compute the reciprocal space representation of a geometry object
## Arguments
#* `geo` :  Geometry object
#* `dnx` : reciprocal space grid in x
#* `dny` : reciprocal space grid in y
## Outputs
#* `F` : reciprocal space representation of the geometry
#"""
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
#Already defined in transforms.jl
#"""
#    drawable(geo::Geometry)
#computes a vector of points for a drawable real-space representation of a geometry object
## Arguments
#* `geo` :  Geometry object
## Outputs
#* `x` : x points of the real-space representation 
#* `y` : y points of the real-space representation 
#"""
function drawable(r::Rectangle)
    return [.5*r.dx,-.5*r.dx,-.5*r.dx,.5*r.dx,.5*r.dx],[.5*r.dy,.5*r.dy,-.5*r.dy,-.5*r.dy,.5*r.dy]
end
function drawable(c::Custom)
    prinln("This method is not yet implemented for Custom geometries")
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

