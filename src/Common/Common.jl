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
	Eigenmodes(V,W,X,q)

Structure to store the eigenmodes of a layer
# Attributes
* `V` : Magnetic field eigenmode
* `W` : Electric field eigenmode
* `X` : Factor for propagating amplitude vector
* `q` : Pseudo wave vector
"""
struct Eigenmodes
    V::AbstractArray{Complex{Float64},2} #Transform towards magnetic fields
    W::AbstractArray{Complex{Float64},2} #Transform towards electric fields
    X::AbstractArray{Complex{Float64},2} #Propagation of the amplitude vector through the layer
    q::AbstractArray{Complex{Float64},2} #pseudo wave vector

end
"""
	Halfspace(V,Kz)

Structure to store the eigenmodes of a layer
# Attributes
* `V` : Magnetic field eigenstate
* `Kz` : z component of the wave vector in the medium
"""
struct Halfspace
    Kz::AbstractArray{Complex{Float64},2} #z-component of the wave vector in the medium
    V::AbstractArray{Complex{Float64},2}  #Transform towards magnetic fields. W is identity anyway.
end
"""
    eigenmodes(dnx,dny,Kx,Ky,λ,l::Layer)
    eigenmodes(g::RcwaGrid,λ,l::Layer)
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
# Outputs
* `em` : Eigenmode object
"""
function eigenmodes(dnx,dny,Kx,Ky,λ,l::PatternedLayer)
    k0=2π/real(λ)
	#get the base permittivity
    εxx=get_permittivity(l.materials[1],λ,1)*I
    if typeof(l.materials[1])<:Isotropic
    else
    end
    #add the permittivity for all inclusions
	if minimum([typeof(m)<:Isotropic for m in l.materials])
		εxx=get_permittivity(l.materials[1],λ)*I
	    for ct=1:length(l.geometries)
    	    rec=reciprocal(l.geometries[ct],dnx,dny)
        	εxx+=rec*(get_permittivity(l.materials[ct+1],λ)-get_permittivity(l.materials[ct],λ))
		end
        εzz=εyy=εxx
        εxy=εyx=0I
	else
		εxx=get_permittivity(l.materials[1],λ,1)*I
        εxy=get_permittivity(l.materials[1],λ,2)*I
        εyx=get_permittivity(l.materials[1],λ,3)*I
        εyy=get_permittivity(l.materials[1],λ,4)*I
        εzz=get_permittivity(l.materials[1],λ,5)*I
	    for ct=1:length(l.geometries)
    	    rec=reciprocal(l.geometries[ct],dnx,dny)
        	εxx+=rec*(get_permittivity(l.materials[ct+1],λ,1)-get_permittivity(l.materials[ct],λ,1))
			εxy+=rec*(get_permittivity(l.materials[ct+1],λ,2)-get_permittivity(l.materials[ct],λ,2))
        	εyx+=rec*(get_permittivity(l.materials[ct+1],λ,3)-get_permittivity(l.materials[ct],λ,3))
        	εyy+=rec*(get_permittivity(l.materials[ct+1],λ,4)-get_permittivity(l.materials[ct],λ,4))
			εzz+=rec*(get_permittivity(l.materials[ct+1],λ,5)-get_permittivity(l.materials[ct],λ,5))
    	end
	end	 	
    #reciprocal of permittivity
    η=I/εzz
    # η=I/εxx
    #Maxwell equations transformed
    #this is old code
	#P=[Kx*η*Ky I-Kx*η*Kx;Ky*η*Ky-I -Ky*η*Kx]
    Q=[Kx*Ky+εyx εyy-Kx*Kx;Ky*Ky-εxx -εxy-Ky*Kx]
    #M=Matrix(P*Q)
	#analytic multiplication can speed things up:
	A=η*(Ky*εyx+Kx*εxx)
	B=η*(Ky*εyy+Kx*εxy)
	M=[Ky.^2-εxx+Kx*A -Ky*Kx-εxy+Kx*B;-Kx*Ky-εyx+Ky*A Kx.^2-εyy+Ky*B]
	#eigenmodes
    ev=eigen(M)
    q=Diagonal(sqrt.(Complex.(ev.values)))
	#select negative root
    q[real.(q).>0].*=-1
    #W is transform between amplitude vector and E-Field
    W=ev.vectors
    #V is transform between amplitude vector and H-Field
    V=Q*W/Diagonal(q)
    #X the factor applied to the amplitudes when propagatin through the layer
    X=exp(Diagonal(q*k0*l.thickness))
    #create struct
    return Eigenmodes(V,W,X,q)
end
function eigenmodes(dnx,dny,Kx,Ky,λ,l::SimpleLayer)
	k0=2π/real(λ)
	#permittivity tensor
    ε=get_permittivity(l.material,λ)*I
	#z component
    Kz=sqrt.(Complex.(ε-Kx*Kx-Ky*Ky))
    q=[1im*Kz 0I;0I 1im*Kz]
    q[real.(q).>0].*=-1
	#magnetic field eigenmodes
    Q=[Kx*Ky ε-Kx*Kx;Ky*Ky-ε -Ky*Kx]
    V=Q/Diagonal(q)
    #W is identity
    W=I+0*V
	#amplitude propagation
    X=exp(Diagonal(q*k0*l.thickness))
    return Eigenmodes(V,W,X,q)
end
function eigenmodes(dnx,dny,Kx,Ky,λ,l::AnisotropicLayer)
    k0=2π/real(λ)
	#permittivity tensor
	εxx=get_permittivity(l.material,λ,1)*I
	#anisotropy
    if typeof(l.material)<:Isotropic
        εzz=εyy=εxx
        εxy=εyx=0I
    else
        εxy=get_permittivity(l.material,λ,2)*I
        εyx=get_permittivity(l.material,λ,3)*I
        εyy=get_permittivity(l.material,λ,4)*I
        εzz=get_permittivity(l.material,λ,5)*I
    end
    η=I/εzz
    #Maxwell equations transformed
    P=[Kx*η*Ky I-Kx*η*Kx;Ky*η*Ky-I -Ky*η*Kx]
    Q=[Kx*Ky+εyx εyy-Kx*Kx;Ky*Ky-εxx -εxy-Ky*Kx]
    #eigenmodes
    ev=eigen(Matrix(P*Q))
    q=Diagonal(sqrt.(Complex.(ev.values)))
	#select negative root
    q[real.(q).>0].*=-1
    #W is transform between amplitude vector and E-Field
    W=ev.vectors
    #V is transform between amplitude vector and H-Field
    V=Q*W/Diagonal(q)
    #X the factor applied to the amplitudes when propagatin through the layer
    X=exp(Diagonal(q*k0*l.thickness))
    #create struct
    return Eigenmodes(V,W,X,q)
end
function eigenmodes(g::RcwaGrid,λ,l::Layer)
    return eigenmodes(g.dnx,g.dny,g.Kx,g.Ky,λ,l)
end
function eigenmodes(g::RcwaGrid,λ,l::Array{Layer,1})
    #initialize array
    rt=Array{Eigenmodes,1}(undef,length(l))
    #iterate through layers
    for cnt=1:length(l)
        rt[cnt]=eigenmodes(g,λ,l[cnt])
    end
    return rt
end
"""
    halfspace(Kx,Ky,material,λ)

Compute the eigenmodes of a halfspace
# Arguments
* `Kx` : kx component of the wave vector in reciprocal space
* `Ky` : ky component of the wave vector in reciprocal space
* `material` : medium
* `λ` : wavelength
# Outputs
* `em` : halfspace eigenmode object
"""
function halfspace(Kx,Ky,material,λ)
    #Base value
	εxx=get_permittivity(material,λ,1)*I
    #probably add off-diagonal elements
	if typeof(material)<:Isotropic
        εzz=εyy=εxx
        εxy=εyx=0I
    else
        εxy=get_permittivity(material,λ,2)*I
        εyx=get_permittivity(material,λ,3)*I
        εyy=get_permittivity(material,λ,4)*I
        εzz=get_permittivity(material,λ,5)*I
    end
	#z component of wave vector
	Kz=sqrt.(Complex.(εzz-Kx*Kx-Ky*Ky))
	#simple solution like free space
    q0=[1im*Kz 0I;0I 1im*Kz]
	#Magnetic field eigenvalues
    Q=[Kx*Ky+εyx εyy-Kx*Kx;Ky*Ky-εxx -εxy-Ky*Kx]
    V=Q/Diagonal(q0)
    return Halfspace(Kz,V)
end

"""
    slicehalf(e)

Utility function, just slices a vector in two vectors of half length
# Arguments
* `v` : vector to be halfed
# Outputs
* `v1` : upper half of v
* `v2` : lower half of v
"""
function slicehalf(v)
    mylength=convert(Int64,size(v,1)/2)
    return v[1:mylength,:],v[mylength+1:end,:]
end

"""
    e2p(ex,ey,ez,Kz,kz0)

Converts the reciprocal-space electric field (in a substrate or superstrate) into a Poynting power flow in z direction
# Arguments
* `ex` : x component of the electric field in reciprocal space
* `ey` : y component of the electric field in reciprocal space
* `ez` : z component of the electric field in reciprocal space
* `Kz` : z component of the wavevector in the medium
* `kz0` : z component of the 0th-order wavevector in the superstrate
# Outputs
* `P` : power flow
"""

function e2p(ex,ey,ez,Kz,kz0)
    #amplitudes squared
    P=abs.(ex).^2+abs.(ey).^2+abs.(ez).^2
    #weight by z component
    P=sum(real.(Kz)*P/real(kz0))
    return P
end


"""
    a2p(a,W,Kx,Ky,Kz,kz0)

Converts an amplitude vector (in substrate or superstrate) to Poynting power flow in z direction
# Arguments
* `a` : amplitude vector
* `W` : eigenmodes of the halfspace
* `Kx` : x component of the wavevector in the medium
* `Ky` : y component of the wavevector in the medium
* `Kz` : z component of the wavevector in the medium
* `kz0` : z component of the plane wave wavevector in the superstrate
# Outputs
* `P` : power flow
"""

function a2p(a,W,Kx,Ky,Kz,kz0)
    ex,ey,ez=a2e(a,W,Kx,Ky,Kz)
    return e2p(ex,ey,ez,Kz,kz0)
end

"""
    a2e2d(a,W)

Converts an amplitude vector to reciprocal-space electric fields Ex and Ey.
This is a "light" version of the "a2e" method.
# Arguments
* `a` : amplitude vector
* `W` : eigenmodes of the halfspace
# Outputs
* `ex` : x-component of the electric field
* `ey` : y-component of the electric field
"""
function a2e2d(a,W)
    e=W*a
	#e contains both components
    ex,ey=slicehalf(e)
    return ex,ey
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
# Outputs
* `ex` : x-component of the electric field
* `ey` : y-component of the electric field
* `ez` : z-component of the electric field
"""

function a2e(a,W,Kx,Ky,Kz)
	ex,ey=a2e2d(a,W)
	#Plane wave, E⊥k, E*k=0
    ez=-Kz\(Kx*ex+Ky*ey)
    return ex,ey,ez
end



end  # module  eigenmodes
