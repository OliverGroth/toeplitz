using LinearAlgebra
using SparseArrays
using ToeplitzMatrices
using BandedMatrices
using Arpack
using Plots
using Roots

function v(A,m) #Requires nxn matrix A where 1 <= m <= n-1
	return A[:,m+1]
end

function f(ff,gg,vv,qq,yy,alpha)
	for j in range(1,alpha):
		F = append(ff,0) - (gg-vv'*ff*)*yy/qq

end

function w(ww,F,H,yy)

end

function y(w)
	return append!(w,-1)
end

function q(a,lmb,vv,ww)
	return a - lmb + vv'*ww
end

function main()
 
end