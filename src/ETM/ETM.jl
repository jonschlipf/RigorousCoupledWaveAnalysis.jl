#Enhanced transmission matrix algorithm by Moharam et al.
module ETM
using LinearAlgebra
using ..Common
export etmSource,etm_reftra,etm_amplitudes,etm_propagate,etm_flow,etm_reftra_flows
function F(em)
    return Matrix([em.W em.W;-em.V em.V])
end
function etm_propagate(ref,tra,em,s,grd,get_r=true)
    B=[I I*0;0*I I; 1im*grd.Kx*grd.Ky/tra.Kz 1im*(grd.Ky^2+tra.Kz^2)/tra.Kz;-1im*(grd.Kx^2+tra.Kz^2)/tra.Kz -1im*grd.Kx*grd.Ky/tra.Kz]
    #backward iteration
    a=Array{Array{Complex{Float64},2},1}(undef,length(em))
    b=Array{Array{Complex{Float64},2},1}(undef,length(em))
    a[end],b[end]=slicehalf(F(em[end])\Matrix(B))
    for cnt=length(em):-1:2
        a[cnt-1],b[cnt-1]=slicehalf(F(em[cnt-1])\F(em[cnt])*[I I*0;I*0 em[cnt].X]*[I ;(b[cnt]/a[cnt])*em[cnt].X])
    end
    A=Matrix([ I 0*I;0*I I;-1im*grd.Kx*grd.Ky/ref.Kz -1im*(ref.Kz^2+grd.Ky*grd.Ky)/ref.Kz;1im*(grd.Kx*grd.Kx+ref.Kz^2)/ref.Kz 1im*grd.Kx*grd.Ky/ref.Kz])
    Bprime=Matrix(F(em[1])*[I I*0;I*0 em[1].X]*[I;(b[1]/a[1])*em[1].X])
    #forward iteratio
    t=Array{Array{Complex{Float64},1},1}(undef,length(em))
    r=Array{Array{Complex{Float64},1},1}(undef,length(em))
    ro,tu=slicehalf(cat(-A,Bprime,dims=2)\s)
    t[1]=vec(tu)
    for cnt=1:length(em)-1
        t[cnt+1]=a[cnt]\I*(em[cnt].X*t[cnt])
		r[cnt]=em[cnt].X*b[cnt]*t[cnt+1]
    end
	r[1]=em[1].X*(b[1]/a[1])*em[1].X*t[1]
    to=a[end]\I*em[end].X*t[end]
	r[end]=em[end].X*b[end]*to
    return ro,to,r,t
end
	
function etm_reftra(s,m::RCWAModel,grd::RcwaGrid,λ,ems,ref,tra)
	kzin=grd.k0[3]*real(sqrt(get_permittivity(m.εsup,λ))
    ro,to,r,t=etm_propagate(ref,tra,ems,s,grd,false)
    R=a2p(ro,I,grd.Kx,grd.Ky,ref.Kz,kzin)
    T=a2p(to,I,grd.Kx,grd.Ky,tra.Kz,kzin)
    return R,T
end
function etm_reftra(s,m::RCWAModel,grd::RcwaGrid,λ)
	ems=eigenmodes(grd,λ,m.layers)
	ref=halfspace(grd.Kx,grd.Ky,m.εsup,λ)
	tra=halfspace(grd.Kx,grd.Ky,m.εsub,λ))
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
	tra=halfspace(grd.Kx,grd.Ky,m.εsub,λ))
	R,T,flw=etm_reftra_flows(s,m,grd,λ,ems,ref,tra)
	return R,T,flw
end
function etm_amplitudes(s,m::RCWAModel,grd::RcwaGrid,λ,ems,ref,tra)
    ro,to,r,t=etm_propagate(ref,tra,ems,s,grd)
	return cat(ro,r,0ro,dims=1),cat(0to,t,to,dims=1)
end	
function etm_amplitudes(s,m::RCWAModel,grd::RcwaGrid,λ)
	ems=eigenmodes(grd,λ,m.layers)
	ref=halfspace(grd.Kx,grd.Ky,m.εsup,λ)
	tra=halfspace(grd.Kx,grd.Ky,m.εsub,λ))
	a,b=etm_amplitudes(s,m,grd,λ,ems,ref,tra)
	return a,b
end
function etm_flow(a,b,em,kzin)
	ex,ey=a2e2d(a+b,em.W)
	hx,hy=a2e2d(-a+b,em.V)
	return imag(sum(ex.*conj.(hy)-ey.*conj.(hx)))/kzin
end

function etmSource(grd,nsup)
	kin=grd.k0*nsup
    width=(grd.Nx*2+1)*(grd.Ny*2+1)
    #vertical
    normal=[0,0,1]
    #te polarization E-field is perpendicular with z-axis and propagation direction (so, parallel with surface)
    kte=cross(normal,kin)/norm(cross(normal,kin))
    #tm polarization E-field is perpendicular with te and propagation direction (so, not necessarily parallel with surface)
    ktm=cross(kin,kte)/norm(cross(kin,kte))
    esource=zeros(width*4)*1im
    esource[convert(Int64,(width+1)/2)]=kte[1]
    esource[convert(Int64,(width+1)/2)+width]=kte[2]
    esource[convert(Int64,(width+1)/2)+2*width]=(kte[2]*kin[3]-kte[3]*kin[2])*1im
    esource[convert(Int64,(width+1)/2)+3*width]=(kte[3]*kin[1]-kte[1]*kin[3])*1im
    ste=esource
    esource=zeros(width*4)*1im
    esource[convert(Int64,(width+1)/2)]=ktm[1]
    esource[convert(Int64,(width+1)/2)+width]=ktm[2]
    esource[convert(Int64,(width+1)/2)+2*width]=(ktm[2]*kin[3]-ktm[3]*kin[2])*1im
    esource[convert(Int64,(width+1)/2)+3*width]=(ktm[3]*kin[1]-ktm[1]*kin[3])*1im
    stm=esource#/sqrt(epsref)
    return ste,stm
end

end
