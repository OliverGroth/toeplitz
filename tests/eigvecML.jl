include("eigenvalueFinderTESTING.jl")
eigVBench(64,2048,3,[6,-4,1],"benchtest.jld2","vectest.jld2",BigFloat)
using GenericSchur

function eigvecforBench()
	VT = eigvecs(matrixMaker(2048,BigFloat))
	@save "vectest_true.jld2" VT
end

function eigvecBench()
	b = @benchmark eigvecforBench()
	@save "eigvecbenchtest.jld2" b
end

eigvecBench()