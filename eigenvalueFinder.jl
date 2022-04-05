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

function f(j,m,lmb)

end

function g(j,m,G)

end


function w(A,m,lmb)

end


function q(A,lmb)
	
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