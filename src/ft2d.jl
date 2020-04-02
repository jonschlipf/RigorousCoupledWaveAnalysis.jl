#2D fourier transform (inefficient, but does the job)
module ft2d
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
<<<<<<< HEAD
#transforms an arbitrary bitmask f to reciprocal space
#Input:
#dnx,dny: reciprocal lattice grid
#f: bitmask to be transformed
#x,y: 2D arrays of the x and y coordinate in real space
#Output: reciprocal space image of the bitmask
function real2recip(dnx,dny,f,x,y)
    #number of reciprocal lattice vectors
    M=Int(sqrt(size(dnx,1)))
    #max difference from zero
    M2=2M-1
    #initialize array of unique values
    Fs=zeros(M2,M2)*1im
    #iterate over reciprocal space
    for i=1:M2
        for j=1:M2
            Fs[i,j]=length(x)^-1*sum(f.*exp.(-2im*π*(x*(i-M)+y*(j-M))))
        end
    end
    #distribute over full reciprocal space map
    F=0.0im*dnx
    for i=1:size(F,1)
        for j=1:size(F,2)
            F[i,j]=Fs[dnx[i,j]+M,dny[i,j]+M]
=======
function real2recip(dnx,dny,f)
    dnx=(dnx.+size(f,1)).%size(f,1).+1
    dny=(dny.+size(f,2)).%size(f,2).+1
    F2=fft(ifftshift(f))/length(f)
    F3=0.0im*dnx
    for i=1:size(F3,1)
        for j=1:size(F3,2)
            F3[i,j]=F2[dnx[i,j],dny[i,j]]
>>>>>>> c69b66c5c2d065e408d1ab7ff2665d81f229fad8
        end
    end
    return F3
end
<<<<<<< HEAD
#transforms a reciprocal space 2D map (for example permittivity) to real space
#Input:
#dnx,dny: reciprocal lattice grid
#F: reciprocal space map to be transformed
#x,y: 2D arrays of the x and y coordinate in real space
#Output: real space image
=======

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

>>>>>>> c69b66c5c2d065e408d1ab7ff2665d81f229fad8
function recip2real(dnx,dny,F,x,y)
    #initialize
    f=zeros(size(x))*1im
    #iterate over real space
    for i=1:size(x,1)
        for j=1:size(y,2)
            #inverse Fourier transform
            f[i,j]=size(dnx,1)^-1*sum(F.*exp.(1im*x[i,j]*dnx*2*pi+1im*y[i,j]*dny*2*pi))
        end
    end
    return f
end

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
end
