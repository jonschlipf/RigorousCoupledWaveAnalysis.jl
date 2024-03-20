#Enhanced transmission matrix algorithm by Moharam et al.
module ETM
using LinearAlgebra
using ..Common
export etm_reftra,etm_amplitudes,etm_propagate,etm_flow,etm_reftra_flows,etm_getfields_stack
#utility method for dealing with eigenmodes
function F(em)
    return Matrix([em.W em.W;em.V -em.V])
end
"""
    etm_propagate(sup,sub,em,ψin,grd)
    etm_propagate(sup,sub,em,ψin,grd,get_r)
The propagation of waves according to the ETM method
# Arguments
* `sup` :  superstrate halfspace object
* `sub` :  substrate halfspace object
* `em` :  array with eigenmodes for all layers
* `ψin` :  incoming (source) amplitude vector
* `grd` : RCWA grid object
* `get_r` : set false if you do not need the internal backward propagating waves (only required for absorption calculation)
# Outputs
* `ψref` : reflected amplitude vector
* `ψtra` : transmitted amplitude vector
* `ψp` : array of internal back propagating wave vectors in each layer
* `ψm` : array of internal forward propagating wave vectors in each layer
"""
function etm_propagate(sup,sub,em,ψin,grd,get_r=true)
    ψm=Array{Array{Complex{Float64},1},1}(undef,length(em)) #preallocate
    ψp=Array{Array{Complex{Float64},1},1}(undef,length(em)) 
if length(em)>0
    #backward iteration
    a=Array{Array{Complex{Float64},2},1}(undef,length(em)) #preallocate
    b=Array{Array{Complex{Float64},2},1}(undef,length(em))
    a[end],b[end]=slicehalf(F(em[end])\Matrix([I;-sub.V])) #transmission matrix for a wave in the last layer into the substrate
    for cnt=length(em):-1:2
        a[cnt-1],b[cnt-1]=slicehalf(F(em[cnt-1])\F(em[cnt])*[em[cnt].X*(a[cnt]/b[cnt])*em[cnt].X ;I]) #successively compute the transmission matrix from the i-th layer into the substrate
    end
	#forward iteration
    ψref,ψm1=slicehalf(-cat([I;sup.V],F(em[1])*[em[1].X*(a[1]/b[1])*em[1].X;I],dims=2)\([I;-sup.V]*ψin)) #compute the reflected wave and the forward wave in the first layer
    ψm[1]=vec(ψm1) 
    for cnt=1:length(em)-1
        ψm[cnt+1]=b[cnt]\I*(em[cnt].X*ψm[cnt]) #successively compute the forward wave in the i-th layer
		get_r&&(ψp[cnt]=em[cnt].X*a[cnt]*ψm[cnt+1]) #if required, compute the backward wave in the i-th layer
    end
	get_r&&(ψp[1]=em[1].X*(a[1]/b[1])*em[1].X*ψm[1]) #if required, computethe backward wave in the first layer
    ψtra=b[end]\I*em[end].X*ψm[end] #compute the transmitted wave from the forward wave in the last layer
	get_r&&(ψp[end]=em[end].X*a[end]*ψtra) #compute the backward wave in the last layer, if required
	else
	ψref,ψtra=slicehalf([-I I;sup.V sub.V]\[I*ψin;sup.V*ψin]) #the case for "empty" model (single interface, no layers)
	end
    return ψref,ψtra,ψp,ψm
end
"""
    etm_reftra(ψin,m,grd,λ)
    etm_reftra(ψin,m,grd,λ,em,sup,sub)
Computes reflection and transmission according to the ETM method
# Arguments
* `ψin` :  incoming (source) amplitude vector
* `m` :  RCWA model object
* `grd` : RCWA grid object
* `λ` : free-space wavelength
* `em` :  array with eigenmodes for all layers (computed from m if not given)
* `sup` :  superstrate halfspace object (computed from m if not given)
* `sub` :  substrate halfspace object (computed from m if not given)
lation)
# Outputs
* `R` : reflection by the device 
* `T` : transmission by the device 
"""
function etm_reftra(ψin,m::RCWAModel,grd::RCWAGrid,λ,em,sup,sub)
	kzin=grd.k0[3]#*real(sqrt(get_permittivity(m.εsup,λ)))
	ro,to,r,t=etm_propagate(sup,sub,em,ψin,grd,false) #propagate amplitudes
    R=a2p(0ro,ro,sup.V,I,kzin) #compute amplitudes to powers
    T=-a2p(to,0to,sub.V,I,kzin)
    return R,T
end
function etm_reftra(s,m::RCWAModel,grd::RCWAGrid,λ)
	ems=eigenmodes(grd,λ,m.layers) #eigenmodes of all layers
	ref=halfspace(grd.Kx,grd.Ky,m.εsup,λ) #superstrate and substrate
	tra=halfspace(grd.Kx,grd.Ky,m.εsub,λ)
	R,T=etm_reftra(s,m,grd,λ,ems,ref,tra)
	return R,T
end
"""
    etm_reftra_flows(ψin,m,grd,λ)
    etm_reftra_flows(ψin,m,grd,λ,em,sup,sub)
Computes reflection and transmission, as well as the net power flows between layers according to the ETM method
# Arguments
* `ψin` :  incoming (source) amplitude vector
* `m` :  RCWA model object
* `grd` : RCWA grid object
* `λ` : free-space wavelength
* `em` :  array with eigenmodes for all layers (computed from m if not given)
* `sup` :  superstrate halfspace object (computed from m if not given)
* `sub` :  substrate halfspace object (computed from m if not given)
lation)
# Outputs
* `R` : reflection by the device 
* `T` : transmission by the device 
* `flw` : array containing the power flows between layers 
"""
function etm_reftra_flows(s,m::RCWAModel,grd::RCWAGrid,λ,ems,sup,sub)
	kzin=grd.k0[3]#*real(sqrt(get_permittivity(m.εsup,λ)))
    ro,to,b,a=etm_propagate(sup,sub,ems,s,grd) #propagate waves
    R=a2p(0ro,ro,sup.V,I,kzin) #reflected and transmitted power
    T=-a2p(to,0to,sub.V,I,kzin)
	flw=[etm_flow(a[i],b[i],ems[i],kzin) for i=1:length(a)] #intermediate power flows
    return R,T,flw
end
function etm_reftra_flows(s,m::RCWAModel,grd::RCWAGrid,λ)
	ems=eigenmodes(grd,λ,m.layers) # layer eigenmodes
	ref=halfspace(grd.Kx,grd.Ky,m.εsup,λ) #superstrate and substrate eigenmodes
	tra=halfspace(grd.Kx,grd.Ky,m.εsub,λ)
	R,T,flw=etm_reftra_flows(s,m,grd,λ,ems,ref,tra)
	return R,T,flw
end
"""
    etm_amplitudes(ψin,m,grd,λ)
    etm_amplitudes(ψin,m,grd,λ,em,sup,sub)
Computes the amplitude vectors of forward and backward propagating waves throught the ETM method
# Arguments
* `ψin` :  incoming (source) amplitude vector
* `m` :  RCWA model object
* `grd` : RCWA grid object
* `λ` : free-space wavelength
* `em` :  array with eigenmodes for all layers (computed from m if not given)
* `sup` :  superstrate halfspace object (computed from m if not given)
* `sub` :  substrate halfspace object (computed from m if not given)  
lation)
# Outputs
* `a` : amplitude vectors of forward waves
* `b` : amplitude vectors of backward waves
"""
function etm_amplitudes(ψin,m::RCWAModel,grd::RCWAGrid,λ,em,sup,sub)
    ro,to,r,t=etm_propagate(sup,sub,em,ψin,grd) #propagate wave
	return cat([ψin],t,[to],dims=1),cat([ro],r,[0ro],dims=1) #put in order
end	
function etm_amplitudes(ψin,m::RCWAModel,grd::RCWAGrid,λ)
	ems=eigenmodes(grd,λ,m.layers) 	#layer eigenmodes
	ref=halfspace(grd.Kx,grd.Ky,m.εsup,λ) #superstrate
	tra=halfspace(grd.Kx,grd.Ky,m.εsub,λ) #substrate
	a,b=etm_amplitudes(ψin,m,grd,λ,ems,ref,tra)
	return a,b
end
"""
    etm_flow(a,b,em,kz0)
Computes the power flow in z direction at a single location
# Arguments
* `a` :  forward wave amplitude vector
* `b` :  backward wave amplitude vector
* `em` :  eigenmode object
* `kz0` :  z-component of the impinging wave vector (for normalization)
# Outputs
* `flow` : power flow in the z direction
"""
function etm_flow(a,b,em,kz0)
	ex,ey=a2e2d(a+b,em.W) #electric field
	hx,hy=a2e2d(a-b,em.V) #magnetic field
	return imag(sum(ex.*conj.(hy)-ey.*conj.(hx)))/kz0 #poynting vector z coordinate
end
"""
    etm_getfields_stack(mdl,grd::RCWAGrid,xypoints,zpoints,λ,a,b,em::Eigenmodes,window,padding)

computes the electric and magnetic fields within the whole stack
# Arguments
* `mdl` : RCWA model object
* `grd` : reciprocal space grid object
* `xypoints` : two-element vector specifying the number of points in x, y, and z for which the fields are to be computed
* `zpoints` : array of desired z-axis points relative to the top interface of the layer
* `λ` : wavelength
* `a` : forward amplitude vectors
* `b` : backward amplitude vectors
* `em` : eigenmodes of the layer
* `window` : type of Window to use against Gibbs' phenomenon, available: Hann (default), None
* `padding` : padding of the window in x and y, default is [0,0]
# Outputs
* `efield` : 4D tensor for the electric field (dimensions are x, y, z, and the component (E_x or E_y or E_z)
* `hfield` : 4D tensor for the magnetic field (dimensions are x, y, z, and the component (E_x or E_y or E_z)
"""
function etm_getfields_stack(mdl::RCWAModel,grd::RCWAGrid,xypoints,zpoints,λ,a,b,em::Array{Eigenmodes,1},window="Hann",padding=[0,0])
    zstack=0
    zind=1
    E=zeros(xypoints[1],xypoints[2],length(zpoints),3)*1im
    H=zeros(xypoints[1],xypoints[2],length(zpoints),3)*1im
    for i in eachindex(mdl.layers)
        zpoints_separated=zpoints[zstack.<=zpoints.<zstack+mdl.layers[i].thickness]
        Ei,Hi=getfields(a[i],b[i],em[i],grd,xypoints,zpoints_separated.-zstack,λ,window,padding)
        E[:,:,zind:zind+length(zpoints_separated)-1,:]=Ei
        H[:,:,zind:zind+length(zpoints_separated)-1,:]=Hi
        zstack=zstack+mdl.layers[i].thickness
        zind=zind+length(zpoints_separated)
    end
    return E,H
end
function etm_getfields_stack(mdl::RCWAModel,grd::RCWAGrid,xypoints,zpoints,λ,source,window="Hann",padding=[0,0])
    em=eigenmodes(grd,λ,mdl.layers)
    ref=halfspace(grd.Kx,grd.Ky,mdl.εsup,λ)
    tra=halfspace(grd.Kx,grd.Ky,mdl.εsub,λ)
    b0,a0,b,a=etm_propagate(ref,tra,em,source,grd)
    E,H=etm_getfields_stack(mdl,grd,xypoints,zpoints,λ,a,b,em,window,padding)
    return E,H
end

end
