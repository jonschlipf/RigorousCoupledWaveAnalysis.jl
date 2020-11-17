export Combination,Rotation,Shift
"""
    Rotation(G,θ)
A geometry object for the description of the rotation of another geometry object G by the angle θ
# Attributes
* `G` : Rotated geometry
* `θ` : Rotation angle
"""
struct Rotation <: Geometry
    G::Geometry
    θ::Float64
end
"""
    Shift(G,ax,ay)
A geometry object for the description of the shift of another geometry object G by ax (relative to the unit cell size) in x and ay in y
# Attributes
* `G` : Rotated geometry
* `ax` : x shift
* `ay` : y shift
"""
struct Shift <: Geometry
    G::Geometry
    ax::Float64
    ay::Float64
end
"""
    Combination(G)
A geometry object for the combination of several geometry objects into one
# Attributes
* `G` : Rotated geometry
"""
struct Combination <: Geometry
    G::Array{Geometry,1}
end
"""
    reciprocal(G,dnx,dny)
Compute the reciprocal space convolution matrix for a given geometry
# Arguments
* `G` :  geometry object
* `dnx` : reciprocal space grid in x
* `dny` : reciprocal space grid in y
# Output
* `ret` : reciprocal space convolution matrix
"""
function reciprocal(c::Rotation,dnx,dny)
	#rotated reciprocal coordinate system
    dnx2=dnx*cos(c.θ)+dny*sin(c.θ)
    dny2=-dnx*sin(c.θ)+dny*cos(c.θ)
    return reciprocal(c.G,dnx2,dny2)
end
function reciprocal(c::Shift,dnx,dny)
	#multiplication by a complex amplitude in reciprocal space
    return reciprocal(c.G,dnx,dny).*exp.(-2im*pi*c.ax.*dnx).*exp.(-2im*pi*c.ay.*dny)
end
function reciprocal(c::Combination,dnx,dny)
	#linearity of Fourier transform
    ret=dnx*0
    for i=1:length(c.G)
        ret+=reciprocal(c.G[i],dnx,dny)
    end
    return ret
end
"""
    drawable(G)
Compute a simple plottable real-space representation a given geometry
# Arguments
* `G` :  geometry object
# Output
* `x` : x coordinates of drawing points
* `y` : y coordinates of drawing points
"""
function drawable(c::Rotation)
	#parent object
    x,y=drawable(c.G)
	#rotated coordinate system
    u=x*cos(-c.θ)+y*sin(-c.θ)
    v=-x*sin(-c.θ)+y*cos(-c.θ)
    return u,v
end
function drawable(c::Shift)
	#parent object
    x,y=drawable(c.G)
	#simple shift
    return x.+c.ax,y.+c.ay
end
function drawable(c::Combination)
	#first object
    x,y=drawable(c.G[1])
    for i=2:length(c.G)
		#concatenate the others
        xon,yon=drawable(c.G[i])
        x=cat(x,xon,dims=1)
        y=cat(y,yon,dims=1)
    end
    return x,y
end
