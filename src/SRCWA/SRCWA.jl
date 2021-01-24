module SRCWA
using LinearAlgebra
export srcwa_reftra,a2p,slicehalf,srcwa_matrices,Srcwa_matrices,srcwa_amplitudes,srcwa_flow

export innerSource,dipoleRad,pointDipole

using ..Common
include("scatterMatrices.jl")
include("emission.jl")



function srcwa_reftra(a0,model::RCWAModel,grd::RCWAGrid,λ)
    ref=halfspace(grd.Kx,grd.Ky,model.εsup,λ)
    tra=halfspace(grd.Kx,grd.Ky,model.εsub,λ)
    S=scattermatrix_ref(ref,grd.V0)
    for ct=1:length(model.layers)
        S=concatenate(S,scattermatrix_layer(eigenmodes(grd.dnx,grd.dny,grd.Kx,grd.Ky,λ,model.layers[ct]),grd.V0))
    end
    S=concatenate(S,scattermatrix_tra(tra,grd.V0))
	kzin=grd.k0[3]
	R=a2p(0a0,S.S11*a0,ref.V,I,kzin)
	T=-a2p(S.S21*a0,0a0,tra.V,I,kzin)
    return R,T
end

function srcwa_reftra(a0,εsup,εsub,S,grd::RCWAGrid,λ)
	ref=halfspace(grd.Kx,grd.Ky,εsup,λ)
    tra=halfspace(grd.Kx,grd.Ky,εsub,λ)
	kzin=grd.k0[3]
	R=a2p(0a0,S.S11*a0,ref.V,I,kzin)
	T=-a2p(S.S21*a0,0a0,tra.V,I,kzin)
    return R,T
end

function srcwa_amplitudes(source,grd::RCWAGrid,mtr::Array{ScatterMatrix,1})
    a=Array{Array{Complex{Float64},1},1}(undef,length(mtr))
    b=Array{Array{Complex{Float64},1},1}(undef,length(mtr))
    a[1]=source
	b[1]=concatenate(mtr).S11*a[1]
	for ct=1:length(a)-1
        Sbefore=concatenate(mtr[1:ct])
        Safter=concatenate(mtr[ct+1:end])
        a[ct+1]=(I-Sbefore.S22*Safter.S11)\(Sbefore.S21*source)
        b[ct+1]=Safter.S11*a[ct+1]
    end
    return a,b
end
function srcwa_amplitudes(source,m::RCWAModel,grd::RCWAGrid,λ)
	mtr=[scattermatrix_layer(eigenmodes(grd,λ,l),grd.V0) for l in m.layers]
	ref=scattermatrix_ref(halfspace(grd.Kx,grd.Ky,m.εsup,λ),grd.V0)
	tra=scattermatrix_ref(halfspace(grd.Kx,grd.Ky,m.εsub,λ),grd.V0)
	mtr=cat(ref,mtr,tra,dims=1)
	a,b=srcwa_amplitudes(source,grd,mtr)
	return a,b
end
function srcwa_flows(a,b,V0,kzin)
    p=zeros(size(a,2))
    for ct=1:size(a,2)
        ex,ey=a2e2d(a[:,ct]+b[:,ct],I)
        hx,hy=a2e2d(a[:,ct]-b[:,ct],V0)
        p[ct]=imag(sum(ex.*conj.(hy)-ey.*conj.(hx)))/kzin
    end
    return p
end



end
