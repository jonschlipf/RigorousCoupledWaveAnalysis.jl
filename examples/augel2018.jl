
#=
Lion Augel, Yuma Kawaguchi, Stefan Bechler, Roman Körner, Jörg Schulze, Hironaga Uchida and Inga A. Fischer,
"Integrated Collinear Refractive Index Sensor with Ge PIN Photodiodes,"
ACS Photonics 2018, 5, 11, 4586-4593 (2018)
=#
using RCWA
Si=InterpolPerm(RCWA.si_schinke) #Si from interpolated literature values
Ge=InterpolPerm(RCWA.ge_nunley) #Ge from interpolated literature values
SiO2=ModelPerm(RCWA.sio2_malitson) #SiO2 from literature dispersion formula

n_H2O=1.321
n_CH3COOH=1.353 #Constant refractive indices
Al=ModelPerm(RCWA.al_rakic) #Al from dispersion formula

N=6 #one needs much larger N (~11 is good, 15 is better) here for accurate results
wls=1100:5:1600 #wavelength array
p=950 #pitch
d=480 #hole diameter
function build_model(n_sup)
	nha=PatternedLayer(100,[Al,ConstantPerm(n_sup^2)],[Circle(d/p)])#patterned NHA layer
	spa=SimpleLayer(50,SiO2)
	nsi=SimpleLayer(20,Si)
	nge=SimpleLayer(20,Ge)
	ige=SimpleLayer(480,Ge)
	return RCWAModel([nha,spa,nsi,nge,ige],ConstantPerm(n_sup^2),Si)
end

A_H2O=zeros(size(wls)) #array for absorption
R_H2O=zeros(size(wls)) #array for reflection
T_H2O=zeros(size(wls)) #array for transmission
R_CH3COOH=zeros(size(wls))
T_CH3COOH=zeros(size(wls))
A_CH3COOH=zeros(size(wls))
for i=1:length(wls)
    λ=wls[i] #get wavelength from array
    grd=rcwagrid(N,N,p,p,1E-5,0,λ) #reciprocal space grid
    #compute for H2O
	ste,stm=rcwasource(grd,n_H2O) #te and tm source amplitudes
	R_H2O[i],T_H2O[i],flw=etm_reftra_flows(ste,build_model(n_H2O),grd,λ) #compute ref, tra and power flows for te
	A_H2O[i]=-flw[end-1]-T_H2O[i] #absorption is the power entering the last layer minus the power leaving the device

	#now same for CH3COOH
    ste,stm=rcwasource(grd,n_CH3COOH)
	R_CH3COOH[i],T_CH3COOH[i],flw=etm_reftra_flows(ste,build_model(n_CH3COOH),grd,λ)
	A_CH3COOH[i]=-flw[end-1]-T_CH3COOH[i]
end
