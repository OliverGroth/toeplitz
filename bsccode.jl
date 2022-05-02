using LinearAlgebra
using GenericLinearAlgebra
using BandedMatrices
using Plots
using Polynomials


# test_intext(101,3,[6.0,-4.0,1.0],t->6 .-8*cos.(t) .+2*cos.(2*t),100000)
function test_intext(n1,alpha,vc,f,nf)
	C=compute_c(n1, alpha, vc, f)
	aeTnf=intext(nf,C,f)
	eTnf=eigvals(toeplitz(nf,vc))
	p1=plot(legend=:none)
	tf=LinRange(pi/(nf+1),nf*pi/(nf+1),nf)
	plot!(p1,tf,log10.(abs.(eTnf-aeTnf)),lw=2)
	display(p1)
end

function toeplitz(n,vc)
	T=eltype(vc)
	Tn=BandedMatrix(0=>vc[1]*ones(T,n))
	for kk=1:length(vc)-1
		Tn=Tn+BandedMatrix(kk=>vc[kk+1]*ones(T,n-kk))
	end
	Symmetric(Tn)
end

function compute_c(n1 :: Integer, alpha :: Integer, vc, f) 
  T=eltype(vc)
  j1 = 1:n1
  E  = zeros(T,alpha,n1)
  hs = zeros(real(T),alpha)
  t1 = LinRange(convert(T,pi)/(n1+1),convert(T,pi)*n1/(n1+1),n1)
  ft1= f(t1)
  for kk = 1:alpha
    nk = 2^(kk-1)*(n1+1)-1
    jk = 2^(kk-1)*j1
    hs[kk] = convert(T,1)/(nk+1)
    E[kk,:] = eigvals(toeplitz(nk,vc))[jk] - ft1
  end
  V = zeros(T,alpha,alpha)
  for ii = 1:alpha, jj = 1:alpha
    V[ii,jj] = hs[ii]^jj
  end
  return C=V\E 
end

function intext(nf,C,f,T=Float64)
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
        # display(m)
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


