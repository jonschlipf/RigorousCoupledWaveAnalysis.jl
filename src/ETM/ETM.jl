#Enhanced transmission matrix algorithm by Moharam et al.
module ETM
using LinearAlgebra
using ..Common
export etmSource,etm_reftra
function F(em)
    return [em.W em.W;-em.V em.V]
end
function  etm_reftra(ref,tra,em,s,grd,kzin)
    B=[I I*0;0*I I; 1im*grd.Kx*grd.Ky/tra.Kz 1im*(grd.Ky^2+tra.Kz^2)/tra.Kz;-1im*(grd.Kx^2+tra.Kz^2)/tra.Kz -1im*grd.Kx*grd.Ky/tra.Kz]
    #backward iteration
    a=Array{Array{Complex{Float64},2},1}(undef,length(em))
    a[end],b=slicehalf(F(em[end])\Matrix(B))
    for cnt=length(em):-1:2
        a[cnt-1],b=slicehalf(F(em[cnt-1])\F(em[cnt])*[I I*0;I*0 em[cnt].X]*[I ;(b/a[cnt])*em[cnt].X])
    end
    A=Matrix([ I 0*I;0*I I;-1im*grd.Kx*grd.Ky/ref.Kz -1im*(ref.Kz^2+grd.Ky*grd.Ky)/ref.Kz;1im*(grd.Kx*grd.Kx+ref.Kz^2)/ref.Kz 1im*grd.Kx*grd.Ky/ref.Kz])
    Bprime=Matrix(F(em[1])*[I I*0;I*0 em[1].X]*[I;(b/a[1])*em[1].X])#*t1
    #forward iteratio
    t=Array{Array{Complex{Float64},1},1}(undef,length(em))
    r,tu=slicehalf(cat(-A,Bprime,dims=2)\s)
    t[1]=vec(tu)
    for cnt=1:length(em)-1
        t[cnt+1]=a[cnt]\(em[cnt].X*t[cnt])
    end
    t=a[end]\em[end].X*t[end]
    R=a2p(r,I,grd.Kx,grd.Ky,ref.Kz,kzin)
    T=a2p(t,I,grd.Kx,grd.Ky,tra.Kz,kzin)
    return R,T
end

function  etm_reftra(s,m::RCWAModel,grd::RcwaGrid,λ)
    ref=halfspace(grd.Kx,grd.Ky,m.εsup,λ)
    tra=halfspace(grd.Kx,grd.Ky,m.εsub,λ)
    ems=eigenmodes(grd,λ,m.layers)
	kzin=grd.k0[3]*real(sqrt(get_permittivity(m.εsup,λ)))
    R,T=etm_reftra(ref,tra,ems,s,grd,kzin)
    return R,T
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
