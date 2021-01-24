function fresnel(n1,n2,θi)
	θt=asind.(n1./n2.*sind.(θi))
	rs=(n1.*cosd.(θi)-n2.*cosd.(θt))./(n1.*cosd.(θi)+n2.*cosd.(θt))
	ts=2n1.*cosd.(θi)./(n1.*cosd.(θi)+n2.*cosd.(θt))
	rp=(n2.*cosd.(θi)-n1.*cosd.(θt))./(n2.*cosd.(θi)+n1.*cosd.(θt))
	tp=2n1.*cosd.(θi)./(n2.*cosd.(θi)+n1.*cosd.(θt))
	return rs,rp,ts,tp
end
function fresnelPower(n1,n2,θi)
	rs,rp,ts,tp=fresnel(real(n1),real(n2),θi)
	θt=asind.(real.(n1)./real.(n2).*sind.(θi))
	Rs=abs(rs)^2	
	Rp=abs(rp)^2	
	Ts=abs(ts)^2*real(n2*cosd(θt))/real(n1*cosd(θi))	
	Tp=abs(tp)^2*real(n2*cosd(θt))/real(n1*cosd(θi))	
	return Rs,Rp,Ts,Tp
end
function fabryperot(n1,n2,n3,k,d,θ1)
	n1=real(n1)
	n3=real(n3)
	θ2=asind.(n1./n2.*sind.(θ1))
	θ3=asind.(n1./n3.*sind.(θ1))
	p=exp(1im*n2*d*k*cosd(θ2))
	r12s,r12p,t12s,t12p=fresnel(n1,n2,θ1)
	r21s,r21p,t21s,t21p=fresnel(n2,n1,θ2)
	r23s,r23p,t23s,t23p=fresnel(n2,n3,θ2)
	r32s,r32p,t32s,t32p=fresnel(n3,n2,θ3)
	
	Cs=t12s./(1 .-r21s.*r23s.*p.^2)
	Bs=r12s+t21s.*r23s.*p.^2 .*Cs
	Fs=t23s.*p.*Cs
	Rs=abs(Bs)^2
	Ts=abs(Fs^2)*real(n3*cosd(θ3))/real(n1*cosd(θ1))
	Cp=t12p./(1 .-r21p.*r23p.*p.^2)
	Bp=r12p+t21p.*r23p.*p.^2 .*Cp
	Fp=t23p.*p.*Cp
	Rp=abs(Bp)^2
	Tp=abs(Fp^2)*real(n3*cosd(θ3))/real(n1*cosd(θ1))
	return Rs,Rp,Ts,Tp
end
