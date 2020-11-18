
#=
Lion Augel, Yuma Kawaguchi, Stefan Bechler, Roman Körner, Jörg Schulze, Hironaga Uchida and Inga A. Fischer,
"Integrated Collinear Refractive Index Sensor with Ge PIN Photodiodes,"
ACS Photonics 2018, 5, 11, 4586-4593 (2018)
=#
using RCWA,LinearAlgebra
si=InterpolPerm(RCWA.si_schinke)
ge=InterpolPerm(RCWA.ge_nunley) #Ge from interpolated measured values
ox=ModelPerm(RCWA.sio2_malitson) #SiO2 from dispersion formula
air=ConstantPerm(1.0)
etoh=ConstantPerm(1.353^2)
al=ModelPerm(RCWA.al_rakic)

N=4 #one needs much larger N here for accurate results
wls=1100:5:1600
p=950
d=480

nha=PatternedLayer(100,[al,etoh],[Circle(d/p)])
spa=SimpleLayer(50,ox)
nsi=SimpleLayer(20,si)
nge=SimpleLayer(20,ge)
ige=SimpleLayer(480,ge)
mdl=RCWAModel([nha,spa,nsi,nge,ige],etoh,si)

nmax=4
A=zeros(length(wls),nmax+1)
R=zeros(length(wls),nmax+1)
T=zeros(length(wls),nmax+1)
function flowabs(a,b,V,W,kzin)
    p=zeros(size(a,2))
    for ct=1:size(a,2)
        ex,ey=RCWA.a2e2d(a[:,ct]+b[:,ct],W)
        hx,hy=RCWA.a2e2d(-a[:,ct]+b[:,ct],V)
        p[ct]=imag(sum(ex.*conj.(hy)-ey.*conj.(hx)))/kzin
    end
    return p
end

nsup=1.353
for N=0:nmax
	println(N)
for i=1:length(wls)
    λ=wls[i] #get wavelength from array
    grd=rcwagrid(N,N,p,p,1E-5,0,λ)
    ste,stm=etmSource(grd,1.353)
println(i)
	R[i,N+1],T[i,N+1],flw=etm_reftra_flows(ste,mdl,grd,λ)
	#println(flw)
	A[i,N+1]=flw[end]-T[i,N+1]
	
end
end
con=[sum((A[:,i+1]-A[:,i]).^2) for i=1:nmax]
