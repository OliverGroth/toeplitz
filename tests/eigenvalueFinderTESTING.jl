using LinearAlgebra
using GenericLinearAlgebra
using GenericSVD
using SparseArrays
using ToeplitzMatrices
using BandedMatrices
using Arpack
using Plots
using Roots
using BenchmarkTools
using Profile
using Polynomials
using JLD2
using DoubleFloats
using Base.Threads
# one big one small indicates vector
# one/two small scalar (or small name (like lmb))
# Function names are one letter to avoid confusion with variables
#include("ML.jl")

function toeplitz(n,vc)
  T=eltype(vc)
  Tn=BandedMatrix(0=>vc[1]*ones(T,n))
  for kk=1:length(vc)-1
    Tn=Tn+BandedMatrix(kk=>vc[kk+1]*ones(T,n-kk))
  end
  Symmetric(Tn)
end

function banded(n,T=Float64) #To create banded matrix
    e1 = ones(T,n-2)
    e4 = -4*ones(T,n-1)
    e6 = 6*ones(T,n)
    A = BandedMatrix(-2 => e1, -1 => e4, 0 => e6, 1 => e4, 2 => e1) #Banded matrix with BigFloat
    return A
end

function matrixMaker(n,T=Float64) #Makes non-sparse nxn bi-Laplace
	return diagm(-2 => ones(T,n-2), -1 => -4*ones(T,n-1), 0 => 6*ones(T,n),
		1 => -4*ones(T,n-1), 2 => ones(T,n-2))
end


function GHFinder(A)
	B = displacements(A)
	S = GenericSVD.svd(B)
	X=S.U[:,1:2]
    D=diagm(S.S[1:2])
    Y=S.Vt[1:2,:]
    alpha = 0
    for s in S.S
    	if s > 10^(-6)
    		alpha += 1
    	end
    end
    return X*D,Y',alpha #G,H,alpha
end

function get_GH(n,vc)
  T=eltype(vc)
  nk=length(vc)
  G=zeros(T,n,2)
  G[1,1]=1
  G[n-nk+2:n,2]=vc[nk:-1:2]
  H=zeros(T,n,2)
  H[1:nk-1,1]=vc[2:nk]
  H[n,2]=-1
  if T == BigFloat
  	return BigFloat.(G),BigFloat.(H)
	end
  return Float64.(G),Float64.(H)
end

function displacements(A)
	T = eltype(A)
	n = size(A)[1]
	Z = diagm(-1 => ones(T,n-1))
	return A*Z - Z*A
end

function q(a,lmb,Vv,Ww)
	# Uses values a(=A[m,m]) and lmb (lambda) aswell as
	# vectors vv and ww (v_{m-1} and w_{m-1}) to create
	# q_m (a value) 
	return a - lmb - Vv'*Ww
end

function y(Ww)
	# Uses Ww (w_{m-1} (a vector) to create y_m (a vector))
	return [Ww;-1]
end

function F(F_old,G,vv,qq,yy)
	# Uses F_old ((m-1) x alpha matrix), G ((m-1) x alpha matrix
	# to extract g_{mj} as 
	# bottom value at column j), vv (vector v_{m-1}), qq (value q_m)
	# and yy (vector y_m)
	T = eltype(G)
	m = (size(F_old)[1]+1)
	alpha = size(F_old)[2]
	Fm = zeros(T,m,alpha)
	gg = G[end,1]
	if m-1 == 1
		Fm[:,1] = [F_old[:,1];0] - yy*(gg-(vv'*F_old[:,1])[1])/qq
	else
		Fm[:,1] = [F_old[:,1];0] - yy*(gg-vv'*F_old[:,1])/qq
	end

	for j in range(2,alpha)
		gg = G[end,j]

		if m-1 == 1
			Fm[:,j] = [F_old[:,j-1];0] - yy*(gg - ((vv'*F_old[:,j-1])[1]))/qq
		else
			Fm[:,j] = [F_old[:,j-1];0] - yy*(gg - (vv'*F_old[:,j-1]))/qq
		end
	end
	return Fm
end

function w(Ww_old,FF,H,Yy)
	# Uses Ww_old (vector w_{m-1}), FF (m x alpha matrix F),
	# H (m x alpha matrix H_m) and Yy (vector y_m)


	return [0;Ww_old] - FF*transpose(H)*Yy
end


function abFinder(a,b,i,A) #(a,b) is starting guess for intervall
	# abFinder måste använda qFinder som ska returnera lista {q_1(lambda), ..., {q_m(lambda)}
	# Vi är nöjda med a och b när:
	# Neg_n(a) = i - 1 och Neg_n(b) = i
	# och
	# q_n(a) > 0 och q_n(b) < 0

	# Vill först kontrollera om Neg_n(a) ≤ i - 1 och Neg_n(b) ≥ i vilket är ett krav för
	# att vi ska kunna hitta våra a och b
	T = eltype(A)
	qFind(lmb) = qFinder(A,lmb)[1]
	Neg_a = count(x->x<-0,qFind(a))
	Neg_b = count(x->x<-0,qFind(b))
	
	if Neg_a <= (i - 1) && Neg_b ≥ i
		# Startgissning på a och b OK! Börja iteration
		N = 1
		maxiter = 10^6
		# antag att för a och b så gäller Neg_n(a) ≤ i - 1 och Neg_n(b) ≥ i
		while N <= maxiter
			#println("Iteration ",N)
			q_a = qFind(a)
			q_b = qFind(b)
			Neg_a = count(x->x<=0,q_a)
			Neg_b = count(x->x<=0,q_b)
			qn_a = q_a[end]
			qn_b = q_b[end]

			if Neg_a == i-1 && Neg_b == i && qn_a > 0 && qn_b < 0
				# Solution found
				return [a, b]
			else
				c = (a+b)/convert(T,2)
				Neg_c = count(x->x<-0,qFind(c))

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
		println("Starting guess of alpha and beta does not satisfy the conditions.")
	end	
end

function qFinder(A,lmb)
	# Finds vector of q_1(lmb),...,q_n(lmb) for a Hermitian 
	# nxn Toeplitz matrix A
	T = eltype(A)
	n = size(A)[1]
	
	if T == Float64
		G,H = get_GH(n,Float64.([6,-4,1]))
	elseif T == BigFloat
		G,H = get_GH(n,BigFloat.([6,-4,1]))
	else
		diplay("Error")
	end
	alpha = size(G)[2]
	q_1 = A[1,1] - lmb
	w_1 = A[1,2] / q_1
	v_1 = A[1,2]

	f_1 = zeros(T,1,alpha)
	for j in range(1,alpha)
		f_1[1,j] = G[1,j] / q_1
	end

	qvector = zeros(T,n)
	qvector[1] = q_1 

	v_prev = v_1
	w_prev = w_1
	f_prev = f_1
	y_m = []
	V_test = []

	for m in range(2,n)
		A_m = A[1:m,1:m] # The A_n matrix cutting off all rows and columns at > m 
		G_m = G[1:m,:] # dropping rows m+1 to n, g_j will be the jth column of G_m
		H_m = H[1:m,:]

		if m != n
			v_m = A[1:m,m+1] # PROBLEM, Will run into problem at m = n since m + 1 will be too large, 
		end
		qvector[m] = q(A[m,m],lmb,v_prev,w_prev) # Vv = V_{m-1}, Ww = W_{m-1}
		y_m = y(w_prev)
		if m == n/2
			V_test = y_m
		end

		f_m = F(f_prev,G_m,v_prev,qvector[m],y_m)
		w_m = w(w_prev, f_m, H_m, y_m)

		w_prev = w_m
		if m != n
			v_prev = v_m
		end
		f_prev = f_m
	end
	return qvector, y_m, V_test
end



function eigFinder(A,I = 0)
	# An eigenvalue finder using algorithm 2.3 combined with bisection to find 
	# interval (alpha,beta) (abFinder). It then uses bisection on this interval 
	# to find the eigenvalue.
	# Takes a matrix and three types of values for I,
	# No value: calculates all eigenvalues
	# I is tuple: calculates values between I[1] and I[2]
	# I is vector: calculates the values specified in vectors
	# The algorithm calculates the ith smallest eigenvalue, thus I specifies
	# which eigenvalue by index that is calculated.

	T = eltype(A)
	
	a = -1 # Specific for Bi-Laplace
	b = 17 # Specific for Bi-Laplace

	qFind(lmb) = qFinder(A,lmb)[1]
	if I == 0
		N = size(A)[1] 	# How many eigenvalues we look for. Here we assume that
		E = zeros(T,N) 	# there are n eigenvalue for an nxn matrix (multiplicites?)
		V = zeros(T,N,N)
		@threads for i in range(1,N)
			x = convert.(T,abFinder(a,b,i,A))	#Interval (alpha,beta) for iterate through
			E[i] = find_zero(lmb->qFind(lmb)[end],(x[1],x[2]))
			V[:,i] = qFinder(A,E[i])[2]
		end
	elseif typeof(I) == Float64 || typeof(I) == Int64
		V = zeros(T,size(A)[1])
		x = convert.(T,abFinder(a,b,I,Float64.(A)))
		E = find_zero(lmb->qFind(lmb)[end],(x[1],x[2]))
		V =  qFinder(A,E)[2] 
	elseif typeof(I) == Tuple{Int64,Int64}
		N = I[2] - I[1] + 1
		E = zeros(T,N)
		V = zeros(T,size(A)[1],N)
		k = I[1] # To keep track of eigenvalue we are on
		for i in range(1,N)
			x = convert.(T,abFinder(a,b,k,A))
			E[i] = find_zero(lmb->qFind(lmb)[end],(x[1],x[2]))
			V[:,i] = qFinder(A,E[i])[2]
			k += 1
		end
	elseif typeof(I) == Vector{Int64} || typeof(I) == UnitRange{Int64} 	# We also allow
		N = length(I)													# a range like a:b
		E = zeros(T,N)
		V = zeros(T,size(A)[1],N)													
		for i in range(1,N)
			x = convert.(T,abFinder(a,b,I[i],A))
			E[i] = find_zero(lmb->qFind(lmb)[end],(x[1],x[2]))
			V[:,i] = qFinder(A,E[i])[2]
		end
	else 
		println("Illegal eigenvalue!")
		return 
	end
	return E,V
end

function MLtest(n1,nf,alpha,T=Float64)
    @. f(t) = 6 - 8*cos(t) + 2*cos(2*t) 
    C = compute_c(n1,alpha,convert.(T,[6.0,-4.0,1.0]),f)
    E_ML = intext(nf,C,f)
    E_true = eigvals(toeplitz(nf,convert.(T,[6.0,-4.0,1.0])))
    return E_ML,E_true
end


function compute_c(n1 :: Integer, alpha :: Integer, vc, f) 
  
  T=eltype(vc)
  if T == Float64
    datatype="Float64"
  elseif T==Double64
    datatype="Double64"
  else
    datatype="BigFloat$(precision(BigFloat))"
  end
  j1 = 1:n1
  E  = zeros(T,alpha,n1)
  hs = zeros(real(T),alpha)
  t1 = LinRange(convert(T,pi)/(n1+1),convert(T,pi)*n1/(n1+1),n1)
  ft1= f(t1)
  for kk = 1:alpha
    nk = 2^(kk-1)*(n1+1)-1
    jk = 2^(kk-1)*j1
    hs[kk] = convert(real(T),1)/(nk+1)
    filename = "eigs/eTn$(nk)$(datatype).jld2"
    #=if isfile(filename)
      @load filename eTn
    else
      eTn=eigvals(toeplitz(nk,vc))
      @save filename eTn
    end=#
    eTn=eigvals(toeplitz(nk,vc))
    E[kk,:] = eTn[jk] - ft1
  end
  V = zeros(real(T),alpha,alpha)
  for ii = 1:alpha, jj = 1:alpha
    V[ii,jj] = hs[ii]^jj
  end
  return C=V\E 
end

function intext(nf,C,f)
  T=eltype(C)
  alpha=size(C,1)
  beta = (alpha+2)*ones(Int64,alpha)
  n1=size(C,2)
  hf=convert(T,1)/(nf+1)
  tf=(1:nf)*hf
  lambdas=zeros(T,nf)
  poly_evals = zeros(T,alpha+1)
  for jj=1:nf
    ell = tf[jj]*(n1+1)
    poly_evals[1]=f(pi*tf[jj])
    for kk = 1:alpha
        m=min(max(floor(Int64,ell-beta[kk]/2),0),n1-beta[kk])
        ii=m+1:m+beta[kk]
        tt=convert.(T,ii)/(n1+1)
        ff=fit(tt,C[kk,ii],beta[kk]-1)
        poly_evals[kk+1]=ff(tf[jj])
    end
    hosymbol=Polynomial(poly_evals)
    lambdas[jj] = hosymbol(hf)
  end
  lambdas
end

function startWithML(A,lmbs)
	k = length(lmbs)
	T = eltype(A)
	V = zeros(T,size(A))
	@threads for i = 1:k 
		V[:,i] = qFinder(A,lmbs[i])[2]
	end
	return V
end

function MLVApprox(n1,nf,alpha,vc,T=Float64)
	vc = convert.(T,vc)
	@. f(t) = 6 - 8*cos(t) + 2*cos(2*t)
	C = compute_c(n1,alpha,vc,f)
	E_ML = intext(nf,C,f)
	A = toeplitz(nf,vc)
	V = startWithML(A,E_ML)
	return V
end

function eigFinderML(n1,nf,alpha,vc,T = Float64)
	# DOES NOT WORK, FAILED TO CONVERGE
	@. f(t) = 6 - 8*cos(t) + 2*cos(2*t)
	vc = convert.(T,vc)
	C = compute_c(n1,alpha,vc,f)
	E_approx = intext(nf,C,f)
	A = toeplitz(nf,vc)
	E = zeros(T,nf)

	qFind(lmb) = qFinder(A,lmb)[1]
	for i=1:nf
		E[i] = find_zero(lmb -> qFind(lmb)[end],E_approx[i])
	end
	return E
end

function eigBench(n,name,T=Float64)
	A = matrixMaker(n,T)
	b = @benchmark eigFinder($A)
	@save name b
end

function eigvecML(n1,nf,alpha,vc,name,T=Float64)
	VE = MLVApprox(n1,nf,alpha,vc,T)
	@save name VE
end

function eigVBench(n1,nf,alpha,vc,name1,name2,T=Float64)
	b = @benchmark eigvecML($n1,$nf,$alpha,$vc,$name2,$T)
	@save name1 b

end

function main(n)
# Runs the solver for the nxn bi-Laplace matrix
	A = matrixMaker(n)
	a = eigmin(A)
	b = eigmax(A)
	B = banded(n)
	X = zeros(2,n)
	E = zeros(n)
	for i in range(1,n)
		x = abFinder(a-1,b+1,i,A)
		E[i] = find_zero(lmb->qFinder(A,lmb)[end],(x[1],x[2]), Bisection())
	end
	return E
end