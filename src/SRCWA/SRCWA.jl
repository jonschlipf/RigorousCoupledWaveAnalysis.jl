module SRCWA
using LinearAlgebra
export srcwa_reftra,a2p,slicehalf,srcwa_matrices,Srcwa_matrices,srcwa_amplitudes,srcwa_flow

export innerSource,dipoleRad,pointDipole

using ..Common
include("scatterMatrices.jl")
include("emission.jl")

"""
    srcwa_reftra(ψin,m,grd,λ)
    srcwa_reftra(ψin,εsup,εsub,S,grd,λ)
Computes reflection and transmission according to the S-matrix method
# Arguments
* `ψin` :  incoming (source) amplitude vector
* `m` :  RCWA model object
* `grd` : RCWA grid object
* `λ` : free-space wavelength
* `S` :  total scattering matrix of the device (computed from m if not given)
* `εsup` :  superstrate material object (taken from m if not given)
* `εsub` :  substrate material object (taken from m if not given)
lation)
# Outputs
* `R` : reflection by the device
* `T` : transmission by the device
"""
function srcwa_reftra(ψin,model::RCWAModel,grd::RCWAGrid,λ)
    ref=halfspace(grd.Kx,grd.Ky,model.εsup,λ) #superstrate and substrate halfspace
    tra=halfspace(grd.Kx,grd.Ky,model.εsub,λ)
    S=scattermatrix_ref(ref,grd.V0) #total scattering matrix of device
    for ct=1:length(model.layers)
        S=concatenate(S,scattermatrix_layer(eigenmodes(grd.dnx,grd.dny,grd.Kx,grd.Ky,λ,model.layers[ct]),grd.V0))
    end
    S=concatenate(S,scattermatrix_tra(tra,grd.V0))
	kzin=grd.k0[3] # impinging wave vector z component
    R=a2p(0ψin,S.S11*ψin,ref.V,I,kzin) # get reflected power from reflected amplitude vector
    T=-a2p(S.S21*ψin,0ψin,tra.V,I,kzin) # get transmitted power form transmitted amplitude vector
    return R,T
end
function srcwa_reftra(ψin,εsup,εsub,S,grd::RCWAGrid,λ)
	ref=halfspace(grd.Kx,grd.Ky,εsup,λ) #superstrate and substrate halfspaces
    tra=halfspace(grd.Kx,grd.Ky,εsub,λ)
	kzin=grd.k0[3] # impinging wave vector z component
	R=a2p(0ψin,S.S11*ψin,ref.V,I,kzin)  # get reflected power from reflected amplitude vector
	T=-a2p(S.S21*ψin,0ψin,tra.V,I,kzin) # get transmitted power form transmitted amplitude vector
    return R,T
end
"""
    srcwa_amplitudes(ψin,m,grd,λ)
    srcwa_amplitudes(ψin,grd,mtr)
Computes the propagating wave amplituded vectors between layer interfaces according to the S-matrix method
# Arguments
* `ψin` :  incoming (source) amplitude vector
* `m` :  RCWA model object
* `grd` : RCWA grid object
* `λ` : free-space wavelength
* `mtr` :  array with scattering matrices includeing superstrate and substrate
lation)  
# Outputs
* `a` : array of forward propagating wave amplitudes
* `b` : array of backward propagating wave amplitudes
"""
function srcwa_amplitudes(ψin,m::RCWAModel,grd::RCWAGrid,λ)
	mtr=[scattermatrix_layer(eigenmodes(grd,λ,l),grd.V0) for l in m.layers] #compute scattering matrices of layers
	ref=scattermatrix_ref(halfspace(grd.Kx,grd.Ky,m.εsup,λ),grd.V0) #superstrate and substrate halfspaces
	tra=scattermatrix_ref(halfspace(grd.Kx,grd.Ky,m.εsub,λ),grd.V0)
	mtr=cat(ref,mtr,tra,dims=1) # store all matrices in an array
	a,b=srcwa_amplitudes(ψin,grd,mtr)
	return a,b
end
function srcwa_amplitudes(ψin,grd::RCWAGrid,mtr::Array{ScatterMatrix,1})
    a=Array{Array{Complex{Float64},1},1}(undef,length(mtr)) #preallocate
    b=Array{Array{Complex{Float64},1},1}(undef,length(mtr))
    a[1]=ψin #incoming wave
	b[1]=concatenate(mtr).S11*a[1] #reflected wave
	for ct=1:length(a)-1
        Sbefore=concatenate(mtr[1:ct]) #total scatter matrix before and after the interface in question
        Safter=concatenate(mtr[ct+1:end])
        a[ct+1]=(I-Sbefore.S22*Safter.S11)\(Sbefore.S21*ψin) #compute intermediate amplitude vectors
        b[ct+1]=Safter.S11*a[ct+1]
    end
    return a,b
end
"""
    srcwa_flows(a,b,em,kz0)
Computes the powers flow in z direction at multiple interfaces
# Arguments
* `a` :  array of forward wave amplitude vectors at all interfaces
* `b` :  array of backward wave amplitude vectors at all interfaces
* `V0` :  magnetic field eigenmodes of free space
* `kz0` :  z-component of the impinging wave vector (for normalization)
# Outputs
* `flow` : array of power flows in the z direction at all interfaces
"""
function srcwa_flows(a,b,V0,kzin)
    p=zeros(size(a,2)) #preallocate
    for ct=1:size(a,2)
        ex,ey=a2e2d(a[:,ct]+b[:,ct],I) #compute electric field
        hx,hy=a2e2d(a[:,ct]-b[:,ct],V0) #compute magnetic field
        p[ct]=imag(sum(ex.*conj.(hy)-ey.*conj.(hx)))/kzin #compute z-component of poynting vector
    end
    return p
end



end
