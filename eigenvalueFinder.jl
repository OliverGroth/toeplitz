using LinearAlgebra
using SparseArrays
using ToeplitzMatrices
using BandedMatrices
using Arpack
using Plots
using Roots

function banded(n) #To create banded matrix
    e1 = ones(n-2)
    e4 = -4*ones(n-1)
    e6 = 6*ones(n)
    A = BandedMatrix(-2 => e1, -1 => e4, 0 => e6, 1 => e4, 2 => e1) #Banded matrix with BigFloat
    return sparse(A)
end

function v(A,m) #Requires nxn matrix A where 1 <= m <= n-1
	return A[:,m+1]
end

function Fm(ff,G,vv,qq,yy,alpha)
	F = zeros(size(ff)[1]+1,alpha)
	gg = last(G[:,1])
	F[1,:] = -yy*gg/qq
	for j in range(2,alpha)
		gg = last(G[end,j])
		F[j,:] = F[j-1,:] - yy*(gg-v'*F[j,:])/qq
	end
	return F
end

function w(ww,F,H,yy)
	return append!(0,ww) - F*transpose(H)*yy
end

function y(w)
	return append!(w,-1)
end

function q(a,lmb,vv,ww)
	return a - lmb + vv'*ww
end

function GHFinder(A)
	B = displacements(A)
	S = svd(B)
	X=S.U[:,1:2]
    D=diagm(S.S[1:2])
    Y=S.Vt[1:2,:]

    #display(X*D)
    #display(Y')
    #display(S.S)
    alpha = 0
    for s in S.S
    	if s > 10^(-6)
    		display(s)
    		alpha += 1
    	end
    end
    println("alpa = ", alpha)
    return X*D,Y'
end

function displacements(A)
	n = size(A)[1]
	Z = diagm(-1 => ones(n-1))
	return A*Z - Z*A
end

function abFinder(a,b) #(a,b) is starting guess for intervall

end

function qFinder(A,lmb)

end

function main(n)
# Runs the solver for the nxn bi-Laplace matrix
	A = banded(n)


end