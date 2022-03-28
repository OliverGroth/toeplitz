using LinearAlgebra
using SparseArrays
using ToeplitzMatrices
using BandedMatrices
using Arpack

function matrix(a)
    A=[a 1 2; 2 a 1; 2 1 a]
    println(A)
end

function first(n)
    e1 = ones(n-2)
    e4 = -4*ones(n-1)
    e6 = 6*ones(n)
    A = spdiagm(-2 => e1, -1 => e4, 0 => e6, 1 => e4, 2 => e1) #Spase diagonal matrix
    return A
end

function second(n)
    e = [1,-4,6,-4,1]
    z = zeros(((n-5)/2))
    w = hcat(z,e)
    w = hcat(w,z)
    A = Toeplitz(1:n,w)
end

function third(n)
    e1 = ones(BigFloat,n-2)
    e4 = -4*ones(BigFloat,n-1)
    e6 = 6*ones(BigFloat,n)
    A = BandedMatrix(-2 => e1, -1 => e4, 0 => e6, 1 => e4, 2 => e1)
    return A
end

function main()
    A = first(10)
    B = third(10)
    display(A)
    display(B)
    display(eigs(A))
end