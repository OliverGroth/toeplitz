using LinearAlgebra
using SparseArrays
using ToeplitzMatrices
using BandedMatrices
using Arpack
using Plots

function matrix(a)
    A=[a 0 2; 2 a 0; 0 1 a]
    display(sparse(A))
end

function RHS(t,y)
    return -cos(y)
end

function odeSolver()
    k = 0.01
    T = 10
    N = convert(Int64,floor(T/k))
    f = 1
    y = zeros(N)
    y[1] = f
    t = LinRange(0,T,N)
    for i in range(1,N-1)
        y[i+1] = y[i] + k*RHS(t[i],y[i])
    end
    plot(t,y)

end


function spd(n)
    e1 = ones(n-2)
    e4 = -4*ones(n-1)
    e6 = 6*ones(n)
    A = spdiagm(-2 => e1, -1 => e4, 0 => e6, 1 => e4, 2 => e1) #Spase diagonal matrix
    return A
end

function tp(n) #Broken
    e = [1,-4,6,-4,1]
    z = zeros(((n-5)/2))
    w = hcat(z,e)
    w = hcat(w,z)
    A = Toeplitz(1:n,w)
end

function banded(n)
    e1 = ones(BigFloat,n-2)
    e4 = -4*ones(BigFloat,n-1)
    e6 = 6*ones(BigFloat,n)
    A = BandedMatrix(-2 => e1, -1 => e4, 0 => e6, 1 => e4, 2 => e1) #Banded matrix with BigFloat
    return sparse(A)
end

<<<<<<< HEAD
##function toeplitz(vc,vr)
 #   n = length(vc)
 #   A = zeros(n,n)
 #   for i in 0:(n-1)
 #       for j in 1:(n-i)
 #       A[i+1,i+1] = vc[j]
 #   end
 #   return A
#end

function toeplitz(col,row)
# Size of result.
    m = length(col[:]);  n = length(row[:])
# Locate the nonzero diagonals.
    [ic,sc] = find(col[:])
    row[1] = 0;  # not used
    [ir,sr] = find(row[:])
# Use spdiags for construction.
    d = inverse([ ir-1 1-ic ])
    B = repeat( [ sr sc ]'.', minimum(m,n),1 )
    T = spdiags( B,d,m,n )
    return T
end

=======
function kronTest(n) #Make nxn bi-Laplace with Kronecker product
    A = Symmetric(-6*I(n) + 
        4*diagm(1 => ones(n-1)) + 
        diagm( 2 => ones(n-2)))
    return sparse(A)
end


>>>>>>> 35280c64806d3db12e349c06a8be89923ec544bb
function main()
    A = spd(10)
    B = banded(10)
    display(A)
    display(B)
    display(eigs(A))
end