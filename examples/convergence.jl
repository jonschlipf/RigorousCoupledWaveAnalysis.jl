
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

nmax=5

A=zeros(length(wls),nmax+1)
R=zeros(length(wls),nmax+1)
T=zeros(length(wls),nmax+1)

nsup=1.353
for N=0:nmax
println(N)
for i=1:length(wls)

    λ=wls[i] #get wavelength from array
    grd=rcwagrid(N,N,p,p,1E-5,0,λ)
    ate,atm=scatterSource(grd,nsup)
    mtr=scatMatrices(mdl,grd,λ)
    a,b=srcwa_amplitudes(atm,grd,mtr)
    flw=srcwa_abs(a,b,grd.V0,grd.k0[3]*1.353)
    R[i,N+1]=1-flw[1]
    T[i,N+1]=flw[end]
    A[i,N+1]=flw[end-1]-flw[end]
end
end
con=[sum((A[:,i+1]-A[:,i]).^2) for i=1:nmax]
