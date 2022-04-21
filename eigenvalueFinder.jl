using LinearAlgebra
using SparseArrays
using ToeplitzMatrices
using BandedMatrices
using Arpack
using Plots
using Roots

# Only one or two big letters indicate matrix
# one big one small indicates vector
# one/two small scalar (or small name (like lmb))
# Function names are one letter to avoid confusion with variables

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

function q(a,lmb,Vv,Ww)
	# Uses values a(=A[m,m]) and lmb (lambda) aswell as
	# vectors vv and ww (v_{m-1} and w_{m-1}) to create
	# q_m (a value) 
	display(a)
	display(lmb)
	display(Vv')
	display(Ww)
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
	m = (size(F_old)[1]+1)
	alpha = size(F_old)[2]
	Fm = zeros(m,alpha)
	gg = G[end,1]
	display(yy)
	display(m)
	Fm[:,1] = -gg*yy/qq
	for j in range(2,alpha)
		gg = G[end,j]
		display([F_old[:,j];0])
		display(vv'*F_old[:,j])
		display(m)
		if m-1 == 1
			Fm[:,j] = [F_old[:,j];0] - yy*(gg - ((vv'*F_old[:,j])[1]))/qq
		else
			Fm[:,j] = [F_old[:,j];0] - yy*(gg - (vv'*F_old[:,j]))/qq
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
	# Finds vector of q_1(lmb),...,q_n(lmb) for a Hermitian 
	# nxn Toeplitz matrix A

	n = size(A)[1]
	
	G,H,alpha = GHFinder(A) # G and H are matrix n x alpha
	q_1 = A[1,1] - lmb
	w_1 = A[1,2] / q_1
	v_1 = A[1,2]

	f_1 = zeros(1,alpha)
	for j in range(1,alpha)
		f_1[1,j] = G[1,j] / q_1
	end
	display(f_1)

	qvector = zeros(n)
	qvector[1] = q_1

	v_prev = v_1
	w_prev = w_1
	f_prev = f_1

	for m in range(2,n)
		A_m = A[1:m,1:m] # The A_n matrix cutting off all rows and columns at > m 
		G_m = G[1:m,:] # dropping rows m+1 to n, g_j will be the jth column of G_m
		H_m = H[1:m,:]

		if m != n
			v_m = A[1:m,m+1] # PROBLEM, Will run into problem at m = n since m + 1 will be too large, 
		end
		qvector[m] = q(A[m,m],lmb,v_prev,w_prev) # Vv = V_{m-1}, Ww = W_{m-1}
		y_m = y(w_prev)
		f_m = F(f_prev,G_m,v_prev,qvector[m],y_m)
		w_m = w(w_prev, f_m, H_m, y_m)

		w_prev = w_m
		if m != n
			v_prev = v_m
		end
		f_prev = f_m
	end
	return qvector
end
	


function main(n)
# Runs the solver for the nxn bi-Laplace matrix
	A = matrixMaker(n)
	display(qFinder(A,10))
end
