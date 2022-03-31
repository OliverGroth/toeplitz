#20220330

using LinearAlgebra
using GenericLinearAlgebra
using ToeplitzMatrices
using BandedMatrices
using SparseArrays
using BenchmarkTools
using KrylovKit
using Latexify

function bl1(n)
	# LinearAlgebra.diagm skapar full matris
	diagm( 0=>6*ones(n),
		     1=>-4*ones(n-1),
		    -1=>-4*ones(n-1),
		     2=>ones(n-2),
		    -2=>ones(n-2))
end

function bl2(n)
	# LinearAlgebra.kron (Kronecker produkt) skapar full matris
	Symmetric(kron(I(n),6)+
		        kron(diagm(1=>ones(n-1)),-4)+
		        kron(diagm(2=>ones(n-2)),1))
end

function krondemo()
	# Exempel på hur kron funkar
	display(kron(I(4),[1 2;3 4]))
	display(kron([1 2;3 4],I(4)))
	display(kron([1 2;3 4],[5 6;7 8]))
	display(kron([5 6;7 8],[1 2;3 4]))
end

function bl3(n)
	# BandedMatrix sparar endast nollskilda element (och matrisen är bandad)
	BandedMatrix( 0=>6*ones(n),
		            1=>-4*ones(n-1),
		           -1=>-4*ones(n-1),
		            2=>ones(n-2),
		            -2=>ones(n-2))     
end

function bl4(n)
	# ToeplitzMatrices skapar också effektivt lagrad datatyp
	Toeplitz(
		vcat([6,-4,1],zeros(n-3)),
		vcat([6,-4,1],zeros(n-3)))
end

function bl5(n)
	# SparseArrays, lagrar bara nollskillda värden (men vet inget om bandade strukturen)
	I=vcat(collect(1:n),collect(1:n-1),collect(1:n-2),collect(2:n),collect(3:n))
	J=vcat(collect(1:n),collect(2:n),collect(3:n),collect(1:n-1),collect(1:n-2))
	V=vcat(6*ones(n),-4*ones(n-1),ones(n-2),-4*ones(n-1),ones(n-2))
	sparse(I,J,V)
end

function bl6(n)
	# SparseArrays, lagrar bara nollskillda värden (men vet inget om bandade strukturen)
	S=spzeros(n,n)
	S[1:2,1:4]=[ 6 -4 1 0; 
	            -4 6 -4 1]
	for ii=3:n-2
		S[ii,ii-2:ii+2]=[1,-4,6,-4,1]
	end
	S[n-1:n,n-3:n]=[1 -4 6 -4;
	                0  1 -4 6]
	return S
end

function bl7(n)
	# Egen Toeplitz funktion, lagrar full matris
	toeplitz(n, [6,-4,1],[6, -4, 1])
end


# en toeplitz funktion
function toeplitz(n :: Integer, vc :: Array{T,1}, vr :: Array{T,1}) where T
    Tn = zeros(eltype(T),n,n)
    nc = min(n,length(vc))
    nr = min(n,length(vr))
    for ii = 1:nc
      Tn = Tn + kron(diagm(-ii+1=>ones(eltype(T),n-ii+1)),vc[ii])
    end
    for jj = 2:nr
      Tn = Tn + kron(diagm( jj-1=>ones(eltype(T),n-jj+1)),vr[jj])
    end
    return Tn
end


function bl8(n)
	# BigFloat datatyp till skillnad från Float64 (double precision)
	BigFloat.(bl6(n))
end


function bigfloatproblem()
	# var försiktiga hur ni hanterar BigFloats
	display(BigFloat(1/3))
	display(BigFloat(1)/3)
end

function typexamples()
	# Olika datatyper av samma vektor
	display(ones(3))
	display(ones(Int64,3))
	display(ones(BigFloat,3))
	display(ones(Rational{BigInt},3))
	display(ones(Complex{Float64},3))
	display(ones(Complex{BigFloat},3))
	display(typeof(ones(Complex{BigFloat},3)))
end



function bench1()
		# Benchmarka egenvärden av en full matris 1000x1000
# långsamt och mycket minne
	n=1000
	Tfull=bl7(n)
	@benchmark  eigvals($Tfull)
end

function bench2()
	
	# Benchmarka egenvärden av en bandad matris 1000x1000
	# notera att vi måste ha symmetrisk matris med Symmetric()
	# snabbare och mindre minne än förra benchmark
	n=1000
	Tbanded=bl3(n)
	@benchmark  eigvals(Symmetric($Tbanded))
end


function bench3()
	# Benchmarka egenvärden av en sparse BigFloat matris 1000x1000
	# notera den måte göras full! med Matrix()
	# lösaren i GenericLinearAlgebra används
	n=100
	Tbigfloat = bl8(n)
	# @benchmark  eigvals($Tbigfloat)
	@benchmark  eigvals(Matrix($Tbigfloat))
end

function bench4()
	# visa hur många decimaler som räknas ut 
	n=100
	Tfull=bl7(n)
	Tbigfloat = bl8(n)
	display(eigvals(Tfull))
	display(eigvals(Matrix(Tbigfloat)))
	display(precision(BigFloat))
	setprecision(BigFloat,1024)
	Tbigfloat = bl8(n)
	display(eigvals(Matrix(Tbigfloat)))
end

function iterativ()
	# Exempel på att använda en iterativ egenvärdeslösare på enspare matris
	# räknar ut ett egenvärde ochdess egenvektor
	# jämför med egenvärdet från bandad matris med vanliga lösaren
	n=1000
	Tsparse=bl6(n)
	ee=eigsolve(Tsparse)
	display(ee[1])
	display(ee[2])
	display(ee[3])
	display(eigvals(Symmetric(bl3(n)))[end])
end
# eigen för egenvektorer också

# Tips: Latexify.jl
function visamatris()
	latexify(toeplitz(6,[1,2,3],[1,4]))
end

function vcathcat()
	a=[1,2,3]
	b=[4,5,6]
	c=[7 8 9]
	display(a)
	display(b)
	display(c)
	display(a')
	display(vcat(a,b))
	display(hcat(a,b))
	append!(a,b)
	display(a)
end


function algorithm23(n,v,j)
	# n storlek på matris
	# v defineirar symmetris symbol v=[6,-4,1] betyder alltså matricen [1 -4 6 -4 1] med 6 i huvuddiagonalen
	# j vilket ehenvärde ska approximeras
end
# Roots.jl
#Polynomails.jl
#wolframalpha.com