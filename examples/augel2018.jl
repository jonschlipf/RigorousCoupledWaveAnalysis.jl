
#=
Lion Augel, Yuma Kawaguchi, Stefan Bechler, Roman Körner, Jörg Schulze, Hironaga Uchida and Inga A. Fischer,
"Integrated Collinear Refractive Index Sensor with Ge PIN Photodiodes,"
ACS Photonics 2018, 5, 11, 4586-4593 (2018)
=#
using RCWA
si=InterpolPerm(si_schinke)
ge=InterpolPerm(ge_nunley) #Ge from interpolated measured values
ox=ModelPerm(sio2_malitson) #SiO2 from dispersion formula
air=ConstantPerm(1.0)
etoh=ConstantPerm(1.353^2)
al=ModelPerm(al_rakic)

N=3 #one needs much larger N here for accurate results
wls=1100:5:1600
p=950
d=480

nha=PatternedLayer(100,[al,etoh],[Circle(d/p)])
spa=SimpleLayer(50,ox)
nsi=SimpleLayer(20,si)
nge=SimpleLayer(20,ge)
ige=SimpleLayer(480,ge)
mdl=RCWAModel([nha,spa,nsi,nge,ige],etoh,si)

A=zeros(size(wls))
R=zeros(size(wls))
T=zeros(size(wls))
R2=zeros(size(wls))
T2=zeros(size(wls))
R3=zeros(size(wls))
T3=zeros(size(wls))
nsup=1.353
for i=1:length(wls)

    λ=wls[i] #get wavelength from array
    grd=rcwagrid(N,N,λ,1E-5,0,p,p)
    ate,atm=scatterSource(grd,nsup)
    mtr=scatMatrices(mdl,grd,λ)
    a,b=srcwa_amplitudes(ate,grd,mtr)
    flw=srcwa_abs(a,b,grd.V0,grd.k0[3]*1.353)
    R2[i],T2[i]=srcwa_reftra(ate,mdl,grd,λ)
    R[i]=1-flw[1]
    T[i]=flw[end]
	ste,stm=etmSource(grd,nsup)
	R3[i],T3[i]=etm_reftra(ste,mdl,grd,λ)
	
    #println(flw)
    #println(R)
    #println(T)
    A[i]=flw[end-1]-flw[end]
end
