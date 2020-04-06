#2D fourier transform (inefficient, but does the job)

using SpecialFunctions,FFTW
export rectft,circft,ellipft,recip2real,recipvec2real,real2recip
#2D Fourier transform of a rectangle
#Input:
#frax,fray: x any y width of the rectangle with respect to the simulation region
#dnx,dny: reciprocal lattice vectors for which the Fourier transform is to be computed
#Output: rectangle transformed to reciprocal space
function rectft(frax,fray,dnx,dny)
    #rectangle transforms to sinc
    return frax*fray*sinc.(frax*dnx).*sinc.(fray*dny)
end
#2D Fourier transform of an ellipse
#Input:
#frax,fray: x and y radius of the ellipse with respect to the simulation region
#dnx,dny: reciprocal lattice vectors for which the Fourier transform is to be computed
#Output: ellipse transformed to reciprocal space
function ellipft(frax,fray,dnx,dny)
    #argument of the bessel function
    radix=.5*sqrt.((dnx*frax).^2+(dny*fray).^2)
    result=besselj.(1,2*pi*radix)./radix
    #asymptotic value for going towards 0
    result[radix.==0].=pi
    #scale amplitude
    return result*frax*fray/4
end

#2D Fourier transform of a circle
#Input:
#fra radius of the circle with respect to the simulation region
#dnx,dny: reciprocal lattice vectors for which the Fourier transform is to be computed
#Output: circle transformed to reciprocal space
function circft(fra,dnx,dny)
    return ellipft(fra,fra,dnx,dny)
end

#transforms an arbitrary bitmask f to reciprocal space
#Input:
#dnx,dny: reciprocal lattice grid
#f: bitmask to be transformed
#x,y: 2D arrays of the x and y coordinate in real space
#Output: reciprocal space image of the bitmask


function real2recip(dnx,dny,f)
    F0=fftshift(fft(f))/length(f)
    a=Int64(ceil(size(f,1)/2))
    b=Int64(ceil(size(f,2)/2))
    F=0.0im*dnx
    for i=1:size(F,1)
        for j=1:size(F,2)
            F[i,j]=F0[dnx[i,j]+a,dny[i,j]+b]
        end
    end
    return F
end

#transforms a reciprocal space 2D map (for example permittivity) to real space
#Input:
#dnx,dny: reciprocal lattice grid
#F: reciprocal space map to be transformed
#x,y: 2D arrays of the x and y coordinate in real space
#Output: real space image


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

function recip2real(dnx,dny,F)
    a=maximum(abs.(dnx)) #the maximum dnx value, =2N
    b=maximum(abs.(dny))
    F0=zeros(2a+1,2b+1)*1im #reduced set with unique values of F
    a=Int64(ceil(size(F0,1)/2))
    b=Int64(ceil(size(F0,2)/2))
    for i=1:size(F0,1) #iterate
        for j=1:size(F0,2)
            indices=findall((dnx.==i-a).&(dny.==j-b)) #find the element with desired dnx and dny
            F0[i,j]=F[indices[1]] #put into new array
        end
    end
    return ifft(ifftshift(F0))*length(F0) #transform to real space domain
end


#function recip2real(dnx,dny,F,x,y)
    #initialize
#    f=zeros(size(x))*1im
    #iterate over real space
#    for i=1:size(x,1)
#        for j=1:size(y,2)
            #inverse Fourier transform
#            f[i,j]=size(dnx,1)^-1*sum(F.*exp.(1im*x[i,j]*dnx*2*pi+1im*y[i,j]*dny*2*pi))
#        end
#    end
#    return f
#end

#transforms a reciprocal space 1D vector (for example electric field) to real space
#Input:
#nx,ny: reciprocal space vector grid
#F: reciprocal space map to be transformed
#x,y: 2D arrays of the x and y coordinate in real space
#Output: real space image
function recipvec2real(nx,ny,F,x,y)
    f=zeros(size(x))*1im
    for i=1:size(x,1)
        for j=1:size(y,2)
            f[i,j]=sum(F.*exp.(1im*x[i,j]*nx*2*pi+1im*y[i,j]*ny*2*pi))
        end
    end
    return f
end
