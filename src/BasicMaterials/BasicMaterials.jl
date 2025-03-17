module BasicMaterials

export ge_nunley,si_schinke,sio2_malitson,zno_bond,al_rakic,ag_rakic,au_rakic,ti_rakic,si3n4_luke

using CSV,Interpolations,DataFrames
include("rakic.jl")

ge_nunley_raw=Matrix(CSV.read(string(@__DIR__,"/ge_nunley.txt"),DataFrame,delim=" "))
ge_nunley_n=cat(ge_nunley_raw[:,1],ge_nunley_raw[:,2]+1im*ge_nunley_raw[:,3],dims=2)
ge_nunley=extrapolate(interpolate((1000ge_nunley_raw[:,1],),ge_nunley_n[:,2].^2,Gridded(Linear())),Flat())


sio2_malitson_n(l)=sqrt(Complex(1+(0.6961663l^2)/(l^2-0.0684043^2)+(0.4079426l^2)/(l^2-0.1162414^2)+(0.8974794l^2)/(l^2-9.896161^2)))
sio2_malitson(l)=sio2_malitson_n(.001l)^2 .+0im
#I. H. Malitson. Interspecimen comparison of the refractive index of fused silica, J. Opt. Soc. Am. 55, 1205-1208 (1965)

si3n4_luke_n(l)=sqrt(Complex(1+(3.0249l^2)/(l^2-0.1353406^2)+(40314l^2)/(l^2-1239.842^2)))
si3n4_luke(l)=si3n4_luke_n(.001l)^2 .+0im
#K. Luke, Y. Okawachi, M. R. E. Lamont, A. L. Gaeta, M. Lipson. Broadband mid-infrared frequency comb generation in a Si3N4 microresonator, Opt. Lett. 40, 4823-4826 (2015)

zno_bond_n(l)=sqrt(2.81418+0.87968l^2 /(l^2-0.3042^2)-0.00711.*l.^2)
zno_bond(l)=zno_bond_n(.001l)^2 .+0im
#W. L. Bond. Measurement of the refractive indices of several crystals, J. Appl. Phys. 36, 1674-1677 (1965)

si_schinke_raw=Matrix(CSV.read(string(@__DIR__,"/si_schinke.txt"),DataFrame,delim=" "))
si_schinke_n=cat(si_schinke_raw[:,1],si_schinke_raw[:,2]+1im*si_schinke_raw[:,3],dims=2)
si_schinke=extrapolate(interpolate((1000si_schinke_raw[:,1],),si_schinke_n[:,2].^2,Gridded(Linear())),Flat())
#C. Schinke, P. C. Peest, J. Schmidt, R. Brendel, K. Bothe, M. R. Vogt, I. Kröger, S. Winter, A. Schirmacher, S. Lim, H. T. Nguyen, D. MacDonald. Uncertainty analysis for the coefficient of band-to-band absorption of crystalline silicon. AIP Advances 5, 67168 (2015)



end # module
