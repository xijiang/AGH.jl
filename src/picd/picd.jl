"""
    function gheader(gt)
Retrieve the header infomation of a packed genotype file.
"""
function gheader(gt)
    nlc, nid, dms, gt_majored = nothing, nothing, nothing, true
    open(gt, "r") do io
        header = mmap(io, Vector{Int64}, 3)
        gt_majored = (header[1] == 1)
        nlc, nid = gt_majored ? header[2:3] : header[[3, 2]]
        dms = Tuple(header[2:3])
    end
    gt_majored, nlc, nid, dms
end

        
"""
    function picd(gt, twop; step = 1000, tol = 1e-4)
Given the binary file of `AGH` genotype format, and a vector of
2 × allele frequencies, this function do the `p`ivoted `i`ncomplete
`C`holesky `d`ecomposition.
Finally it returns an approximate G-inverse.

The program 
"""
function picd(gt::String, twop; step = 1000, tol = 1e-4)
    goffset = 24
    gt_majored, nlc, nid, dms = gheader(gt)
    gt_majored || error("The matrix must be genotype majored for speed reasons")
    blk = zeros(dms[1], step)
    # 1/√(2pq), as product is faster than divide.
    rs2pq = 1 ./ sqrt.(twop .* (1 .- twop .* .5))
    open(gt, "r") do io
        gt = mmap(io, Matrix{Int8}, dms, goffset)
        # tp = mean(gt, dims=2) # if gt_majored, or, 1.
        copyto!(blk, view(gt, :, 1:step))
        blk .-= twop
        blk .*= rs2pq
        G = blk'blk
        G, piv, rank, info = LAPACK.pstrf!('L', G, tol)
    end
end

function tst_picd(gt, twop; step = 1000, tol = 1e-4)
    rs2pq = 1 ./ sqrt.(twop .* (1 .- twop .* .5))
    blk = zeros(size(gt)[1], step)
    copyto!(blk, view(gt, :, 1:step))
    blk .-= twop
    blk .*= rs2pq
    G = blk'blk
    G, piv, rank, info = LAPACK.pstrf!('L', G, tol)
end

"""
    function cholesky_decomposition!(A)
This function is just for fun, and should not be taken seriously.
The lower triangle of matrix `A` will be replace by `L`, such that
the original `A = LL'`.
"""
function cholesky_decomposition!(A)
    issymmetric(A) || error("Not a symmetric matrix")
    m = size(A)[1]
    for i in 1:m
        for k in 1:i-1
            A[i, i] -= A[i, k] * A[i, k]
        end
        A[i, i] = sqrt(A[i, i]) # diagonals
        for j in i+1:m
            for k in 1:i-1
                A[j, i] -= A[i, k] * A[j, k]
            end
            A[j, i] /= A[i, i]
        end                     # off diagonals
    end
end
