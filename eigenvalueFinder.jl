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

function matrixMaker(n) #Makes non-sparse nxn bi-Laplace
	return diagm(-2 => ones(n-2), -1 => -4*ones(n-1), 0 => 6*ones(n),
		1 => -4*ones(n-1), 2 => ones(n-2))
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
		F[j,:] = F[j-1,:] .- yy*(gg.-vv'*F[j,:])/qq
	end
	return F
end

function w(ww,F,H,yy)
	return [0;ww] .- F*transpose(H)*yy
end

function y(w)
	return [w -1]
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

function abFinder(a,b,i,A) #(a,b) is starting guess for intervall
	# abFinder måste använda qFinder som ska returnera lista {q_1(lambda), ..., {q_m(lambda)}
	# Vi är nöjda med a och b när:
	# Neg_n(a) = i - 1 och Neg_n(b) = i
	# och
	# q_n(a) > 0 och q_n(b) < 0

	# Vill först kontrollera om Neg_n(a) ≤ i - 1 och Neg_n(b) ≥ i vilket är ett krav för
	# att vi ska kunna hitta våra a och b

	Neg_a = count(x->x<-0,qFinder(A,a))
	Neg_b = count(x->x<-0,qFinder(A,b))
	
	if Neg_a <= (i - 1) && Neg_b ≥ i
		# Startgissning på a och b OK! Börja iteration
		N = 1
		maxiter = 10^6
		# antag att för a och b så gäller Neg_n(a) ≤ i - 1 och Neg_n(b) ≥ i
		while N <= maxiter
			Neg_a = count(x->x<=0,qFinder(A,a))
			Neg_b = count(x->x<=0,qFinder(A,b))
			qn_a = qFinder(A,a)[end]
			qn_b = qFinder(A,b)[end]

			if Neg_a == i-1 && Neg_b == i && qn_a > 0 && qn_b < 0
				# Solution found
				return [a, b]
			else
				c = (a+b)/2
				Neg_c = count(x->x<-0,qFinder(A,c))

				if Neg_c ≤ i-1
					a = c
					N += 1
				elseif Neg_c ≥ i
					b = c
					N += 1
				else
					println("Unknown error in abFinder.")
				end
			end
		end
		println("Too many iterations while trying to find alpha and beta.")
	else
		println("Starting guess of alpha and beta does not satisty the conditions.")
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
	println(F)
	Q = zeros(n)
	Q[1] = q1

	for m in range(2,n)
		Am = A[1:m,1:m]
		Gm = Gn[1:m,1:m]
		Hm = Hn[1:m,1:m]
		vm = A[1:m-1,m][1]

		q = A[m,m] - lmb .- vm'*wm
		Q[m] = q
		ym = y(wm)
		F = Fm(F,Gm,vm',q,ym,alpha)
		wm = w(wm,F,Hm,ym)
	end
	return Q


end

function main(n)
# Runs the solver for the nxn bi-Laplace matrix
	A = banded(n)


end