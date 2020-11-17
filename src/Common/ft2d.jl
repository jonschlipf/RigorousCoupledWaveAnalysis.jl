#2D fourier transform (inefficient, but does the job)

using SpecialFunctions,FFTW
export recip2real,recipvec2real,real2recip
"""
	real2recip(dnx,dny,f)
Converts a 2D real space pattern into a reciprocal space convolution matrix
# Arguments
* `dnx` : reciprocal space grid in x
* `dny` : reciprocal space grid in y
* `f` : binary real-space representation of a pattern
#Output
* `F`: reciprocal space convolution of f
"""
function real2recip(dnx,dny,f)
    F0=fftshift(fft(f))/length(f)
    a=Int64(ceil(size(f,1)/2+.5))
    b=Int64(ceil(size(f,2)/2+.5))
    F=0.0im*dnx
    for i=1:size(F,1)
        for j=1:size(F,2)
            F[i,j]=F0[dnx[i,j]+a,dny[i,j]+b]
        end
    end
    return F
end
"""
	recip2real(dnx,dny,F)
	recip2real(dnx,dny,F,sx,sy)
Converts a reciprocal space convolution matrix into a 2D real space pattern. sx and sy can be specified to increase the image resolution. If not specified, an fft is carried out directly, and the result has the same size as the reciprocal space representation.
# Arguments
* `dnx` : reciprocal space grid in x
* `dny` : reciprocal space grid in y
* `F`: reciprocal space convolution matrix
* `sx`: number of pixels in the result along the x axis
* `sy`: number of pixels in the result along the y axis
#Output
* `f` : binary real-space representation corresponding to F
"""
function recip2real(dnx,dny,F)
    a=maximum(abs.(dnx)) #the maximum dnx value, =2N
    b=maximum(abs.(dny))
    F0=zeros(2a+1,2b+1)*1im #reduced set with unique values of F
    a=Int64(ceil(size(F0,1)/2+.5))
    b=Int64(ceil(size(F0,2)/2+.5))
    for i=1:size(F0,1) #iterate
        for j=1:size(F0,2)
            indices=findall((dnx.==i-a).&(dny.==j-b)) #find the element with desired dnx and dny
            F0[i,j]=F[indices[1]] #put into new array
        end
    end
    return ifft(ifftshift(F0))*length(F0) #transform to real space domain
end

function recip2real(dnx,dny,F,sx,sy)
    a=maximum(abs.(dnx)) #the maximum dnx value, =2N
    b=maximum(abs.(dny))
    F0=zeros(sx,sy)*1im #reduced set with unique values of F
    a=Int64(ceil(size(F0,1)/2+.5))
    b=Int64(ceil(size(F0,2)/2+.5))
    for i=1:size(F0,1) #iterate
        for j=1:size(F0,2)
            indices=findall((dnx.==i-a).&(dny.==j-b)) #find the element with desired dnx and dny
            if length(indices)>0
                F0[i,j]=F[indices[1]] #put into new array
            end
        end
    end
    return ifft(ifftshift(F0))*length(F0) #transform to real space domain
end
#legacy
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
"""
	recipvec2real(nx,ny,F,x,y)
Converts a reciprocal space amplitude vector into a 2D real space map.
# Arguments
* `nx` : reciprocal space grid in x
* `ny` : reciprocal space grid in y
* `F`: reciprocal space amplitude vector
* `x`: x coordinates of desired points
* `y`: y coordinates of desired points
#Output
* `f` : 2D real space map
"""

function recipvec2real(nx,ny,F,x,y)
    f=zeros(size(x))*1im
    for i=1:size(x,1)
        for j=1:size(y,2)
			#this is a manual and awful FT, should be replaced by FFTW
            f[i,j]=sum(F.*exp.(1im*x[i,j]*nx*2*pi+1im*y[i,j]*ny*2*pi))
        end
    end
    return f
end
