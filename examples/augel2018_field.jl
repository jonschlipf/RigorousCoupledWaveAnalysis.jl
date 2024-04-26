
#=
Lion Augel, Yuma Kawaguchi, Stefan Bechler, Roman Körner, Jörg Schulze, Hironaga Uchida and Inga A. Fischer,
"Integrated Collinear Refractive Index Sensor with Ge PIN Photodiodes,"
ACS Photonics 2018, 5, 11, 4586-4593 (2018)
=#
using RigorousCoupledWaveAnalysis
Si=InterpolPerm(RigorousCoupledWaveAnalysis.si_schinke) #Si from interpolated literature values
Ge=InterpolPerm(RigorousCoupledWaveAnalysis.ge_nunley) #Ge from interpolated literature values
SiO2=ModelPerm(RigorousCoupledWaveAnalysis.sio2_malitson) #SiO2 from literature dispersion formula

n_H2O=1.321
n_CH3COOH=1.353 #Constant refractive indices
Al=ModelPerm(RigorousCoupledWaveAnalysis.al_rakic) #Al from dispersion formula

N=6 #one needs much larger N (~11 is good, 15 is better) here for accurate results
λ=1300
p=950 #pitch
d=480 #hole diameter

function build_model(n_sup)
	sup=ConstantPerm(n_sup^2)
	top=SimpleLayer(100,sup)
	nha=PatternedLayer(100,[Al,ConstantPerm(n_sup^2)],[Circle(d/p)])#patterned NHA layer
	spa=SimpleLayer(50,SiO2)
	nsi=SimpleLayer(20,Si)
	nge=SimpleLayer(20,Ge)
	ige=SimpleLayer(480,Ge)
	return RCWAModel([top,nha,spa,nsi,nge,ige],sup,Si)
end


grd=rcwagrid(N,N,p,p,1E-5,0,λ,ConstantPerm(n_H2O^2),true) #reciprocal space grid
ste,stm=rcwasource(grd) #te and tm source amplitudes
mdl=build_model(n_H2O)
#number of desired data points
xypoints=[100,100]
#z axis
zpoints=0:4.5:500
#compute fields
E,H=RigorousCoupledWaveAnalysis.etm_getfields_stack(mdl,grd,xypoints,zpoints,λ,stm)
