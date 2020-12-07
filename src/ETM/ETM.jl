#Enhanced transmission matrix algorithm by Moharam et al.
module ETM
using LinearAlgebra
using ..Common
export etm_reftra,etm_amplitudes,etm_propagate,etm_flow,etm_reftra_flows
function F(em)
    return Matrix([em.W em.W;em.V -em.V])
end
function etm_propagate(ref,tra,em,ψin,grd,get_r=true)
    #backward iteration
    a=Array{Array{Complex{Float64},2},1}(undef,length(em))
    b=Array{Array{Complex{Float64},2},1}(undef,length(em))
    a[end],b[end]=slicehalf(F(em[end])\Matrix([I;-tra.V]))
    for cnt=length(em):-1:2
        a[cnt-1],b[cnt-1]=slicehalf(F(em[cnt-1])\F(em[cnt])*[em[cnt].X*(a[cnt]/b[cnt])*em[cnt].X ;I])
    end
	#forward iteratio
    ψm=Array{Array{Complex{Float64},1},1}(undef,length(em))
    ψp=Array{Array{Complex{Float64},1},1}(undef,length(em))
    ψref,ψm1=slicehalf(-cat([I;ref.V],F(em[1])*[em[1].X*(a[1]/b[1])*em[1].X;I],dims=2)\([I;-ref.V]*ψin))
    ψm[1]=vec(ψm1)
    for cnt=1:length(em)-1
        ψm[cnt+1]=b[cnt]\I*(em[cnt].X*ψm[cnt])
		get_r&&ψp[cnt]=em[cnt].X*a[cnt]*ψm[cnt+1]
    end
	get_r&&ψp[1]=em[1].X*(a[1]/b[1])*em[1].X*ψm[1]
    ψtra=b[end]\I*em[end].X*ψm[end]
	ψp[end]=em[end].X*a[end]*ψtra
    return ψref,ψtra,ψp,ψm
end
	
function etm_reftra(s,m::RCWAModel,grd::RcwaGrid,λ,em,ref,tra)
	kzin=grd.k0[3]*real(sqrt(get_permittivity(m.εsup,λ)))
	ro,to,r,t=etm_propagate(ref,tra,em,s,grd,false)
    R=a2p(ro,I,grd.Kx,grd.Ky,ref.Kz,kzin)
    T=a2p(to,I,grd.Kx,grd.Ky,tra.Kz,kzin)
    return R,T
end
function etm_reftra(s,m::RCWAModel,grd::RcwaGrid,λ)
	ems=eigenmodes(grd,λ,m.layers)
	ref=halfspace(grd.Kx,grd.Ky,m.εsup,λ)
	tra=halfspace(grd.Kx,grd.Ky,m.εsub,λ)
	R,T=etm_reftra(s,m,grd,λ,ems,ref,tra)
	return R,T
end
	
function etm_reftra_flows(s,m::RCWAModel,grd::RcwaGrid,λ,ems,ref,tra)
	kzin=grd.k0[3]*real(sqrt(get_permittivity(m.εsup,λ)))
    ro,to,b,a=etm_propagate(ref,tra,ems,s,grd)
    R=a2p(ro,I,grd.Kx,grd.Ky,ref.Kz,kzin)
    T=a2p(to,I,grd.Kx,grd.Ky,tra.Kz,kzin)
	flw=[etm_flow(a[i],b[i],ems[i],kzin) for i=1:length(a)]
    return R,T,flw
end
function etm_reftra_flows(s,m::RCWAModel,grd::RcwaGrid,λ)
	ems=eigenmodes(grd,λ,m.layers)
	ref=halfspace(grd.Kx,grd.Ky,m.εsup,λ)
	tra=halfspace(grd.Kx,grd.Ky,m.εsub,λ)
	R,T,flw=etm_reftra_flows(s,m,grd,λ,ems,ref,tra)
	return R,T,flw
end
function etm_amplitudes(s,m::RCWAModel,grd::RcwaGrid,λ,em,ref,tra)
    ro,to,r,t=etm_propagate(ref,tra,em,s,grd)
	return cat(ro,r,0ro,dims=1),cat(0to,t,to,dims=1)
end	
function etm_amplitudes(s,m::RCWAModel,grd::RcwaGrid,λ)
	ems=eigenmodes(grd,λ,m.layers)
	ref=halfspace(grd.Kx,grd.Ky,m.εsup,λ)
	tra=halfspace(grd.Kx,grd.Ky,m.εsub,λ)
	a,b=etm_amplitudes(s,m,grd,λ,ems,ref,tra)
	return a,b
end
function etm_flow(a,b,em,kzin)
	ex,ey=a2e2d(a+b,em.W)
	hx,hy=a2e2d(a-b,em.V)
	return -imag(sum(ex.*conj.(hy)-ey.*conj.(hx)))/kzin
end

end
