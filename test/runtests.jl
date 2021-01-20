using RCWA
using Test
@testset "VerticalIncidence" begin
n1=rand()*10
n2=rand()*10

mdl=RCWAModel([],ConstantPerm(n1^2),ConstantPerm(n2^2))
grd=rcwagrid(0,0,100,100,1E-3,rand()*360,1000)
ste,stm=rcwasource(grd,n1)
R,T=etm_reftra(ste,mdl,grd,1000)
@test ((n1-n2)/(n1+n2))^2≈R
@test 1-((n1-n2)/(n1+n2))^2≈T
end
