using LinearAlgebra
using SparseArrays
using ToeplitzMatrices
using BandedMatrices

function matrix(a)
    A=[a 1 2; 2 a 1; 2 1 a]
    println(A)
end

function first(n)
    e1 = ones(n-4)
    e4 = -4*ones(n-3)
    e6 = 6*ones(n)
    A = spdiagm(-2 => e1, -1 => e4, 0 => e6, 1 => e4, 2 => e1)
    println(A)
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
    e1 = ones(n-4)
    e4 = -4*ones(n-3)
    e6 = 6*ones(n)
    A = BandedMatrix(-2 => e1, -1 => e4, 0 => e6, 1 => e4, 2 => e1)
    A = sparse(A)
    println(A)
    return A
end

function main()
    A = first(10)
    B = third(10)
    return A==B
end