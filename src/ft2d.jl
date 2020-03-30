module ft2d
using SpecialFunctions,FFTW
export rectft,circft,ellipft,recip2real,recipvec2real,real2recip

function rectft(frax,fray,dnx,dny)
    #rectangle transforms to sinc
    return frax*fray*sinc.(frax*dnx).*sinc.(fray*dny)
end

function ellipft(frax,fray,dnx,dny)
    #argument of the bessel function
    radix=.5*sqrt.((dnx*frax).^2+(dny*fray).^2)
    result=besselj.(1,2*pi*radix)./radix
    #asymptotic value for going towards 0
    result[radix.==0].=pi
    #scale
    return result*frax*fray/4
end

function circft(fra,dnx,dny)
    return ellipft(fra,fra,dnx,dny)
end
function real2recip(dnx,dny,f)
    dnx=(dnx.+size(f,1)).%size(f,1).+1
    dny=(dny.+size(f,2)).%size(f,2).+1
    F2=fft(ifftshift(f))/length(f)
    F3=0.0im*dnx
    for i=1:size(F3,1)
        for j=1:size(F3,2)
            F3[i,j]=F2[dnx[i,j],dny[i,j]]
        end
    end
    return F3
end

#function real2recip(dnx,dny,f,x,y)
#    M=Int(sqrt(size(dnx,1)))
#    M2=2M-1
#    Fs=zeros(M2,M2)*1im
#    for i=1:M2
#        for j=1:M2
#            Fs[i,j]=length(x)^-1*sum(f.*exp.(-2im*π*(x*(i-M)+y*(j-M))))
#        end
#    end
#    F=0.0im*dnx
#    for i=1:size(F,1)
#        for j=1:size(F,2)
#            F[i,j]=Fs[dnx[i,j]+M,dny[i,j]+M]
#        end
#    end
#    return F
#end

function recip2real(dnx,dny,F,x,y)
    f=zeros(size(x))*1im
    for i=1:size(x,1)
        for j=1:size(y,2)
            f[i,j]=size(dnx,1)^-1*sum(F.*exp.(1im*x[i,j]*dnx*2*pi+1im*y[i,j]*dny*2*pi))
        end
    end
    return f
end

function recipvec2real(nx,ny,F,x,y)
    f=zeros(size(x))*1im
    for i=1:size(x,1)
        for j=1:size(y,2)
            f[i,j]=sum(F.*exp.(1im*x[i,j]*nx*2*pi+1im*y[i,j]*ny*2*pi))
        end
    end
    return f
end
end
