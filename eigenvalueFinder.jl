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

function f(f,g,v,f,q,y)

end

function g(j,m,G)

end

function w(w,F,H,y)

end

function y(w)

end


function q(A,lmb,v,w)
	
	q_1 = A[1,1] - lmb
	N = size(A)[1]
	Q = zeros(N)
	Q[1] = q_1
	for m in range(2,N)
		q = A[m,m] - lmb - v(A,m-1)
		Q[m] = q
	end
	return Q

end

function main()
 
end