
using RigorousCoupledWaveAnalysis
Si=InterpolPerm(RigorousCoupledWaveAnalysis.si_schinke)
air=ConstantPerm(1.0)



N=6
mdl=RCWAModel([PatternedLayer(100,[air,Si],[Circle(.1)])],air,air)
grd=rcwagrid_pml(N,N,10000,10000,1e-5,0,1000,mdl.εsup,.5)
ste,stm=rcwasource(grd)
xypoints=[100,100]
zpoints=0:4.5:500
E,H=RigorousCoupledWaveAnalysis.etm_getfields_stack(mdl,grd,xypoints,zpoints,1000.0,stm)
