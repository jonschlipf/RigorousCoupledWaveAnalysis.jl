
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
