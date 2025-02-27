using LinearAlgebra
using LinearAlgebra,CUDA
using ..Common

export squarescale_reftra

function taylorx(dnx,dny,Kx,Ky,λ,l::PatternedLayer)
    IMa=Diagonal(Kx*0 .+1)
    IM=[IMa 0*IMa;0*IMa IMa]
    k0=2π/real(λ)
    #get the base permittivity
    εxx=get_permittivity(l.materials[1],λ,1)*I
    #add the permittivity for all inclusions
    if minimum([typeof(m)<:Common.Isotropic for m in l.materials])
        #all isotropic
        εxx=get_permittivity(l.materials[1],λ)*I
        for ct=1:length(l.geometries)
            rec=reciprocal(l.geometries[ct],dnx,dny)
            εxx+=rec*(get_permittivity(l.materials[ct+1],λ)-get_permittivity(l.materials[ct],λ))
        end
        εzz=εyy=εxx
        εxy=εyx=0I
    else
        #anisotropic
        εxx=get_permittivity(l.materials[1],λ,1)*I
        εxy=get_permittivity(l.materials[1],λ,2)*I
        εyx=get_permittivity(l.materials[1],λ,3)*I
        εyy=get_permittivity(l.materials[1],λ,4)*I
        εzz=get_permittivity(l.materials[1],λ,5)*I
        for ct=1:length(l.geometries)
            rec=reciprocal(l.geometries[ct],dnx,dny)
            εxx+=rec*(get_permittivity(l.materials[ct+1],λ,1)-get_permittivity(l.materials[ct],λ,1))
            εxy+=rec*(get_permittivity(l.materials[ct+1],λ,2)-get_permittivity(l.materials[ct],λ,2))
            εyx+=rec*(get_permittivity(l.materials[ct+1],λ,3)-get_permittivity(l.materials[ct],λ,3))
            εyy+=rec*(get_permittivity(l.materials[ct+1],λ,4)-get_permittivity(l.materials[ct],λ,4))
            εzz+=rec*(get_permittivity(l.materials[ct+1],λ,5)-get_permittivity(l.materials[ct],λ,5))
        end
    end	 	
    a01=0
    a11=-0.10036558103014462001
    a21=-0.00802924648241156960
    a31=-0.00089213849804572995

    b01=0
    b11=0.39784974949964507614
    b21=1.36783778460411719922
    b31=0.49828962252538267755
    b41=-0.00063789819459472330

    b02=-10.9676396052962062593
    b12=1.68015813878906197182
    b22=0.05717798464788655127
    b32=-0.00698210122488052084
    b42=0.00003349750170860705

    b03=-0.09043168323908105619
    b13=-0.06764045190713819075
    b23=0.06759613017704596460
    b33=0.02955525704293155274
    b43=-0.00001391802575160607

    b04=0
    b14=0
    b24=-0.09233646193671185927
    b34=-.01693649390020817171
    b44=-0.00001400867981820361
    η=inv(εzz)
    P=[Kx*η*Ky I-Kx*η*Kx;Ky*η*Ky-I -Ky*η*Kx]
    Q=[Kx*Ky+εyx εyy-Kx*Kx;Ky*Ky-εxx -εxy-Ky*Kx]
    A0=[0IM P;Q 0IM]*k0*l.thickness
    nrm=maximum(sum(abs.(A0),dims=1))
    m=Int(ceil(log2(nrm)))
    A=A0*2.0^-m
    PQ=P*Q
    PQP=PQ*P
    QP=Q*P
    QPQ=QP*Q
    A2=[PQ 0IM;0IM QP]*(k0*l.thickness*2.0^-m)^2
    A3=[0IM PQP;QP*Q 0IM]*(k0*l.thickness*2.0^-m)^3
    A6=[PQP*QPQ 0IM;0IM QPQ*PQP]*(k0*l.thickness*2.0^-m)^6
    #A2=A*A
    #A3=A2*A
    #A6=A3*A3
    B1=a01*I+a11*A+a21*A2+a31*A3
    B2=b01*I+b11*A+b21*A2+b31*A3+b41*A6
    B3=b02*I+b12*A+b22*A2+b32*A3+b42*A6
    B4=b03*I+b13*A+b23*A2+b33*A3+b43*A6
    B5=b04*I+b14*A+b24*A2+b34*A3+b44*A6
    A9=B1*B5+B4
    X=B2+(B3+A9)*A9
    return X^(2^m)
end
function squarescalex(dnx,dny,Kx,Ky,λ,l::PatternedLayer)
    IMa=Diagonal(Kx*0 .+1)
    IM=[IMa 0*IMa;0*IMa IMa]
    k0=2π/real(λ)
    #get the base permittivity
    εxx=get_permittivity(l.materials[1],λ,1)*I
    #add the permittivity for all inclusions
    if minimum([typeof(m)<:Common.Isotropic for m in l.materials])
        #all isotropic
        εxx=get_permittivity(l.materials[1],λ)*I
        for ct=1:length(l.geometries)
            rec=reciprocal(l.geometries[ct],dnx,dny)
            εxx+=rec*(get_permittivity(l.materials[ct+1],λ)-get_permittivity(l.materials[ct],λ))
        end
        εzz=εyy=εxx
        εxy=εyx=0I
    else
        #anisotropic
        εxx=get_permittivity(l.materials[1],λ,1)*I
        εxy=get_permittivity(l.materials[1],λ,2)*I
        εyx=get_permittivity(l.materials[1],λ,3)*I
        εyy=get_permittivity(l.materials[1],λ,4)*I
        εzz=get_permittivity(l.materials[1],λ,5)*I
        for ct=1:length(l.geometries)
            rec=reciprocal(l.geometries[ct],dnx,dny)
            εxx+=rec*(get_permittivity(l.materials[ct+1],λ,1)-get_permittivity(l.materials[ct],λ,1))
            εxy+=rec*(get_permittivity(l.materials[ct+1],λ,2)-get_permittivity(l.materials[ct],λ,2))
            εyx+=rec*(get_permittivity(l.materials[ct+1],λ,3)-get_permittivity(l.materials[ct],λ,3))
            εyy+=rec*(get_permittivity(l.materials[ct+1],λ,4)-get_permittivity(l.materials[ct],λ,4))
            εzz+=rec*(get_permittivity(l.materials[ct+1],λ,5)-get_permittivity(l.materials[ct],λ,5))
        end
    end	 	
    η=inv(εzz)
    P=[Kx*η*Ky I-Kx*η*Kx;Ky*η*Ky-I -Ky*η*Kx]
    Q=[Kx*Ky+εyx εyy-Kx*Kx;Ky*Ky-εxx -εxy-Ky*Kx]
    A0=[0IM P;Q 0IM]*k0*l.thickness
    nrm=maximum(sum(abs.(A0),dims=1))
    m=Int(ceil(log2(nrm)))
    m=0
    A=A0*2.0^-m
    X=exp(A)
    return X^(2^m)
end
function taylor_reftra(ψin,m::RCWAModel,grd::RCWAGrid,λ)
    IMa=Diagonal(grd.Kx*0 .+1)
    IM=[IMa 0*IMa;0*IMa IMa]
    IMb=[IM 0*IM;0*IM IM]
    X=[taylorx(grd.dnx,grd.dny,grd.Kx,grd.Ky,λ,l) for l in m.layers]
    Xp=IMb
    for X in X
        Xp*=X
    end
    ref=halfspace(grd.Kx,grd.Ky,m.εsup,λ) #superstrate and substrate
    tra=halfspace(grd.Kx,grd.Ky,m.εsub,λ)
    Y=Xp*[IM;-tra.V]
    S=[Y [-IM;-ref.V]]\[IM;-ref.V]*ψin
    to,ro=slicehalf(S)

    kzin=grd.k0[3]#*real(sqrt(get_permittivity(m.εsup,λ)))
    R=a2p(0ro,ro,ref.V,IM,kzin) #compute amplitudes to powers
    T=-a2p(to,0to,tra.V,IM,kzin)

    return R,T

end
function squarescale_reftra(ψin,m::RCWAModel,grd::RCWAGrid,λ)
    IMa=Diagonal(grd.Kx*0 .+1)
    IM=[IMa 0*IMa;0*IMa IMa]
    IMb=[IM 0*IM;0*IM IM]
    X=[squarescalex(grd.dnx,grd.dny,grd.Kx,grd.Ky,λ,l) for l in m.layers]
    Xp=IMb
    for X in X
        Xp*=X
    end
    ref=halfspace(grd.Kx,grd.Ky,m.εsup,λ) #superstrate and substrate
    tra=halfspace(grd.Kx,grd.Ky,m.εsub,λ)
    Y=Xp*[IM;-tra.V]
    S=[Y [-IM;-ref.V]]\[IM;-ref.V]*ψin
    to,ro=slicehalf(S)

    kzin=grd.k0[3]#*real(sqrt(get_permittivity(m.εsup,λ)))
    R=a2p(0ro,ro,ref.V,IM,kzin) #compute amplitudes to powers
    T=-a2p(to,0to,tra.V,IM,kzin)

    return R,T

end
