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
    return X*D,Y',alpha #G,H,alpha
end

function displacements(A)
	n = size(A)[1]
	Z = diagm(-1 => ones(n-1))
	return A*Z - Z*A
end

function abFinder(a,b,TOL) #(a,b) is starting guess for intervall
	# abFinder måste använda qFinder som ska returnera lista {q_1(lambda), ..., {q_m(lambda)}
	# Vi är nöjda med a och b när:
	# Neg_n(a) = i - 1 och Neg_n(b) = i
	# vilket är ekvivalent med 
	# q_n(a) > 0 och q_n(b) < 0
	N = 1
	maxiter = 10^6
	# a och b ska utifrån startgissningarna tas fram med hjälp av bisektion vilket är en "root finder"-method
	while N <= maxiter
		c = (a+b)/q2
		q_nc = qFinder(A,c)[end]
		if q_nc == 0 || (b-a)/2 < TOL
			# solution found
			return [a, b] # osäker på om [a, b] eller c är mest användbart som return
		else
			N += 1
			if sign(a) == sign(c)
				a = c
			else
				b = c
	end
end

function qFinder(A,lmb)
	n = size(A)[1]
	GHa = GHFinder(A) #Finds G, H, and alpha
	Gn = GHa[1]
	Hn = GHa[2]
	alpha = GHa[3]

	q1 = A[1,1] - lmb
	wm = A[1,2]/q1
	F = zeros(alpha)
	for j in range(1,alpha)
		F[j] = Gn[1,j]/q1
	end

	Q = zeros(n)
	Q[1] = q1
	q2 = a[2,2] - lmb - A[1,2]'*wm
	Q[2] = q2
	ym = [wm -1]

	for m in range(2,n)
		Am = A[1:m,1:m]
		Gm = Gn[1:m,1:m]
		Hm = Hn[1:m,1:m]
		vm = A[1:m-1,m]
		

		q = A[m,m] - lmb .- vm'*wm
		Q[m] = q
		ym = y(w)
		F = Fm(F,Gm,vm',q,y,alpha)
		wm = w(wm,F,Hm,y)
	end
	return Q


end

function main(n)
# Runs the solver for the nxn bi-Laplace matrix
	A = banded(n)


end