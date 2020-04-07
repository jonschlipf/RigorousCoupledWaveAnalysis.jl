
"""
Lion Augel, Yuma Kawaguchi, Stefan Bechler, Roman Körner, Jörg Schulze, Hironaga Uchida and Inga A. Fischer,
"Integrated Collinear Refractive Index Sensor with Ge PIN Photodiodes,"
ACS Photonics 2018, 5, 11, 4586-4593 (2018)
"""
using RCWA
si=InterpolPerm(si_schinke)
ge=InterpolPerm(ge_nunley) #Ge from interpolated measured values
ox=ModelPerm(sio2_malitson) #SiO2 from dispersion formula
air=ConstantPerm(1.0)
etoh=ConstantPerm(1.353^2)
al=ModelPerm(al_rakic)

N=4
wls=1100:5:1600
p=950
d=480

nha=PatternedLayer(100,[al,etoh],[Circle(d/p)])
spa=PlainLayer(50,ox)
nsi=PlainLayer(20,si)
nge=PlainLayer(20,ge)
ige=PlainLayer(480,ge)
mdl=RCWAModel([nha,spa,nsi,nge,ige],etoh,si)

for i=1:1#length(wls)
    λ=wls[i] #get wavelength from array
    grd=rcwagrid(mdl,N,N,λ,1E-5,0,p,p)
    ate,atm=scatterSource(grd.kin,N,N)
    mtr=scatMatrices(mdl,grd,λ)
    a,b=srcwa_amplitudes(ate,grd,mtr)
    p=srcwa_abs(a,b,grd::RcwaGrid)
    println(length(p))
end
