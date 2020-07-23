#Common methodes for computations with eigenmodes and transforms

module  Common
using LinearAlgebra
include("ft2d.jl")
include("materials.jl")
include("models.jl")
include("grids.jl")
export Eigenmodes,Halfspace
export eigenmodes,halfspace
export a2e,a2e2d,a2p,e2p,slicehalf
"""
Structure to store the eigenmodes of a layer
"""
struct Eigenmodes
    V::AbstractArray{Complex{Float64},2} #Transform towards magnetic fields
    W::AbstractArray{Complex{Float64},2} #Transform towards electric fields
    X::AbstractArray{Complex{Float64},2} #Propagation of the amplitude vector through the layer
    q::AbstractArray{Complex{Float64},2} #pseudo wave vector

end
"""
Structure to store eigenmodes of the upper and lower halfspace (substrate and superstrate)
"""
struct Halfspace
    Kz::AbstractArray{Complex{Float64},2} #z-component of the wave vector in the medium
    V::AbstractArray{Complex{Float64},2}  #Transform towards magnetic fields. W is identity anyway.
end
"""
    eigenmodes(g::RcwaGrid,λ,l::Layer)
    eigenmodes(dnx,dny,Kx,Ky,k0,λ,l::Layer)
    eigenmodes(g::RcwaGrid,λ,l::Array{Layer,1})
Compute the eigenmodes of a layer
# Arguments
* `g` :  grid object
* `λ` :  free space wavelength
* `l` :  layer object
* `dnx` : reciprocal space grid in x
* `dny` : reciprocal space grid in y
* `Kx` : kx component of the wave vector in reciprocal space
* `Ky` : ky component of the wave vector in reciprocal space
* `k0` : free space wavevector (for grid computation)
"""
function eigenmodes(dnx,dny,Kx,Ky,k0,λ,l::PatternedLayer)
    #get the base permittivity
    εxx=Kx*0+get_permittivity(l.materials[1],λ,1)*I
    if typeof(l.materials[1])<:Isotropic
        εzz=εyy=εxx
        εxy=εyx=0εxx
    else
        εxy=Kx*0+get_permittivity(l.materials[1],λ,2)*I
        εyx=Kx*0+get_permittivity(l.materials[1],λ,3)*I
        εyy=Kx*0+get_permittivity(l.materials[1],λ,4)*I
        εzz=Kx*0+get_permittivity(l.materials[1],λ,5)*I
    end
    #add the permittivity for all inclusions
    for ct=1:length(l.geometries)
        rec=reciprocal(l.geometries[ct],dnx,dny)
        εxx+=rec*(get_permittivity(l.materials[ct+1],λ,1)-get_permittivity(l.materials[ct],λ,1))
        εxy+=rec*(get_permittivity(l.materials[ct+1],λ,2)-get_permittivity(l.materials[ct],λ,2))
        εyx+=rec*(get_permittivity(l.materials[ct+1],λ,3)-get_permittivity(l.materials[ct],λ,3))
        εyy+=rec*(get_permittivity(l.materials[ct+1],λ,4)-get_permittivity(l.materials[ct],λ,4))
        εzz+=rec*(get_permittivity(l.materials[ct+1],λ,5)-get_permittivity(l.materials[ct],λ,5))
    end
    #reciprocal of permittivity
    η=I/εzz
    # η=I/εxx
    #Maxwell equations transformed
    P=[Kx*η*Ky I-Kx*η*Kx;Ky*η*Ky-I -Ky*η*Kx]
    Q=[Kx*Ky+εyx εyy-Kx*Kx;Ky*Ky-εxx -εxy-Ky*Kx]
    #Q=[Kx*Ky εxx-Kx*Kx;Ky*Ky-εxx -Ky*Kx]
    #eigenmodes
    ev=eigen(Matrix(P*Q))
    q=Diagonal(sqrt.(Complex.(ev.values)))
    q[real.(q).>0].*=-1
    #W is transform between amplitude vector and E-Field
    W=ev.vectors
    #V is transform between amplitude vector and H-Field
    V=Q*W/Diagonal(q)
    #X the factor applied to the amplitudes when propagatin through the layer
    X=exp(q*k0*l.thickness)
    #create struct
    return Eigenmodes(Matrix(V),Matrix(W),X,q)
end
function eigenmodes(dnx,dny,Kx,Ky,k0,λ,l::SimpleLayer)
    ε=0Kx+get_permittivity(l.material,λ)*I
    Kz=sqrt.(Complex.(ε-Kx*Kx-Ky*Ky))
    Q=[Kx*Ky ε-Kx*Kx;Ky*Ky-ε -Ky*Kx]
    q=[1im*Kz zeros(size(Kz));zeros(size(Kz)) 1im*Kz]
    q[real.(q).>0].*=-1
    #W is identity
    V=Q/Diagonal(q)
    W=I+0*V
    X=exp(Matrix(q*k0*l.thickness))
    return Eigenmodes(Matrix(V),Matrix(W),X,q)
end
function eigenmodes(dnx,dny,Kx,Ky,k0,λ,l::AnisotropicLayer)
    εxx=Kx*0+get_permittivity(l.material,λ,1)*I
    if typeof(l.material)<:Isotropic
        εzz=εyy=εxx
        εxy=εyx=0εxx
    else
        εxy=Kx*0+get_permittivity(l.material,λ,2)*I
        εyx=Kx*0+get_permittivity(l.material,λ,3)*I
        εyy=Kx*0+get_permittivity(l.material,λ,4)*I
        εzz=Kx*0+get_permittivity(l.material,λ,5)*I
    end
    η=I/εzz
    #Maxwell equations transformed
    P=[Kx*η*Ky I-Kx*η*Kx;Ky*η*Ky-I -Ky*η*Kx]
    Q=[Kx*Ky+εyx εyy-Kx*Kx;Ky*Ky-εxx -εxy-Ky*Kx]
    #eigenmodes
    ev=eigen(Matrix(P*Q))
    q=Diagonal(sqrt.(Complex.(ev.values)))
    q[real.(q).>0].*=-1
    #W is transform between amplitude vector and E-Field
    W=ev.vectors
    #V is transform between amplitude vector and H-Field
    V=Q*W/Diagonal(q)
    #X the factor applied to the amplitudes when propagatin through the layer
    X=exp(q*k0*l.thickness)
    #create struct
    return Eigenmodes(Matrix(V),Matrix(W),X,q)
end
function eigenmodes(g::RcwaGrid,λ,l::Layer)
    return eigenmodes(g.dnx,g.dny,g.Kx,g.Ky,g.k0,λ,l)
end
function eigenmodes(g::RcwaGrid,λ,l::Array{Layer,1})
    #initialize array
    rt=Array{Eigenmodes,1}(undef,length(l))
    #iterate
    for cnt=1:length(l)
        rt[cnt]=eigenmodes(g,λ,l[cnt])
    end
    return rt
end
"""
    halfspace(Kx,Ky,ε)

Compute the eigenmodes of a halfspace
# Arguments
* `Kx` : kx component of the wave vector in reciprocal space
* `Ky` : ky component of the wave vector in reciprocal space
* `ε` : complex permittivity of the halfspace material
"""
function halfspace(Kx,Ky,material,λ)
    εxx=Kx*0+get_permittivity(material,λ,1)*I
    if typeof(material)<:Isotropic
        εzz=εyy=εxx
        εxy=εyx=0εxx
    else
        εxy=Kx*0+get_permittivity(material,λ,2)*I
        εyx=Kx*0+get_permittivity(material,λ,3)*I
        εyy=Kx*0+get_permittivity(material,λ,4)*I
        εzz=Kx*0+get_permittivity(material,λ,5)*I
    end
    Kz=sqrt.(Complex.(εzz-Kx*Kx-Ky*Ky))
    Q=[Kx*Ky+εyx εyy-Kx*Kx;Ky*Ky-εxx -εxy-Ky*Kx]

    q0=[1im*Kz zeros(size(Kz));zeros(size(Kz)) 1im*Kz]
    V=Q/Diagonal(q0)
    return Halfspace(Kz,V)
end

"""
    slicehalf(e)

Just slices a vector in two vectors of half length
# Arguments
* `e` : vector to be halfed
"""
function slicehalf(e)
    mylength=convert(Int64,size(e,1)/2)
    return e[1:mylength,:],e[mylength+1:end,:]
end

"""
    e2p(ex,ey,ez,Kz,k0)

Converts the reciprocal-space electric field (in a substrate or superstrate) into a power flow
# Arguments
* `ex` : x component of the electric field in reciprocal space
* `ey` : y component of the electric field in reciprocal space
* `ez` : z component of the electric field in reciprocal space
* `Kz` : z component of the wavevector in the medium
* `kz0` : z component of the plane wave wavevector in the superstrate
"""

function e2p(ex,ey,ez,Kz,kz0)
    #amplitudes squared
    P=abs.(ex).^2+abs.(ey).^2+abs.(ez).^2
    #weight by z component
    P=sum(real.(Kz)*P/real(kz0))
    return P
end


"""
    a2p(a,W,Kx,Ky,Kz,k0)

converts an amplitude vector (in substrate or superstrate) to power flow
# Arguments
* `a` : x amplitude vector
* `W` : eigenmodes of the halfspace
* `Kx` : x component of the wavevector in the medium
* `Ky` : y component of the wavevector in the medium
* `Kz` : z component of the wavevector in the medium
* `kz0` : z component of the plane wave wavevector in the superstrate
"""

function a2p(a,W,Kx,Ky,Kz,kz0)
    ex,ey,ez=a2e(a,W,Kx,Ky,Kz)
    return e2p(ex,ey,ez,Kz,kz0)
end
"""
    a2e(a,W,Kx,Ky,Kz)

converts an amplitude vector (in substrate or superstrate) to reciprocal-space electric fields
# Arguments
* `a` : x amplitude vector
* `W` : eigenmodes of the halfspace
* `Kx` : x component of the wavevector in the medium
* `Ky` : y component of the wavevector in the medium
* `Kz` : z component of the wavevector in the medium
"""

function a2e(a,W,Kx,Ky,Kz)
    e=W*a
    ex,ey=slicehalf(e)
    ez=-Kz\(Kx*ex+Ky*ey)
    return ex,ey,ez
end
"""
    a2e2d(a,W)

Converts an amplitude vector to reciprocal-space electric fields Ex and Ey.
This is a "light" version of the "a2e" method.
# Arguments
* `a` : x amplitude vector
* `W` : eigenmodes of the halfspace
"""
function a2e2d(a,W)
    e=W*a
    ex,ey=slicehalf(e)
    return ex,ey
end



end  # module  eigenmodes
