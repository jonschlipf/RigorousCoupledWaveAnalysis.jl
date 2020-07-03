function pointDipole(grd,px0,py0,pz0)
	fc=(2pi/grd.k0)^2/grd.px/grd.py#prefactor from fourier transform
	jx=1im*ones(length(grd.nx))*px0*sqrt(fc)#this is the normalized j, not the normalized one
	jy=1im*ones(length(grd.nx))*py0*sqrt(fc)
	jz=1im*ones(length(grd.nx))*pz0*sqrt(fc)
	return jx,jy,jz
end
function innerSource(grd,jx,jy,jz,Sin,Sout,eps)
	#P=2pi/grd.k0/sqrt(grd.px*grd.py)
	#P=P*
	P=[-1im*grd.Kx*(eps\jz);-1im*grd.Ky*(eps\jz);jy;-jx]
	M=[I+Sout.S11 -Sin.S22-I;grd.V0*(I-Sout.S11) grd.V0*(-Sin.S22+I)]
	prvec=M\P
	a0,b0=slicehalf(prvec)
	return a0,b0
end
function dipoleRad(a0,b0,Sout,grd,sub,jx,jy,jz)
	exb,eyb,ezb=a2e((I+Sout.S11)*a0,I,grd.Kx,grd.Ky,sub.Kz)
	#sz=[sqrt(length(grd.nx)),sqrt(length(grd.nx))]
	#x=[r  for r in -sz[1]/2+.5:sz[1]/2-.5, c in -sz[2]/2+.5:sz[2]/2-.5]/sz[1]
	#y=[c  for r in -sz[1]/2+.5:sz[1]/2-.5, c in -sz[2]/2+.5:sz[2]/2-.5]/sz[2]
	#pxr=recipvec2real(grd.nx,grd.ny,px,x,y)
	#pyr=recipvec2real(grd.nx,grd.ny,py,x,y)
	#pzr=recipvec2real(grd.nx,grd.ny,pz,x,y)
	#exr=recipvec2real(grd.nx,grd.ny,exb,x,y)
	#eyr=recipvec2real(grd.nx,grd.ny,eyb,x,y)
	#ezr=recipvec2real(grd.nx,grd.ny,ezb,x,y)
	#ru=(conj.(pxr).*exr+conj.(pyr).*eyr+conj.(pzr).*ezr)#/sqrt(get_permittivity(m1,2pi/grd.k0))
	#au=.5*2pi/grd.k0/sqrt(grd.px*grd.py)*sum(ru)/length(ru)
	#au=.5sum(ru)/length(ru)
	au=1im*.5sum(conj.(jx).*exb+conj.(jy).*eyb+conj.(jz).*ezb)#the 1im factor is to denormalize j
	return au
end
