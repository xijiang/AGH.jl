"""
    function dpca(mf; k = 100)
Diagonal pivoted cross approximation.
This is algorithm 2 of https://arxiv.org/pdf/1505.06195.pdf.

File `mf` contains an Int64 `n` indicating the following `Matrix{Float64}`
of `n×n`.
It is supposed that the matrix is too large to fit in memory.
`Mmap` is used to deal with it.

The function returns `A` of `n×k`, such that `|M - A×A'|_∞ ≤ ϵ`, 
and `piv`, which is the pivoted order.
"""
function dpca(mf; k = 100)
    dim = begin
        t = zeros(Int64, 1)
        read!(mf, t)
        t[1]
    end
    G = mmap(mf, Matrix{Float64}, (dim, dim), 8)
    A = zeros(dim, k)           # to be returned
    d = diag(G)
    #for l in 1:k                # && on ϵ/tol, to be added
    #    i = argmax(d)
    #    γ = d[i]
    #    s = 
    #    A[:, l] = (G[:, i] - 
    ##for i in 1:k
    #ϵ = argmax(d)
    #γ = d[ϵ]
    #
    #@debug "Debug:" ϵ
end

"""
    function ordinary_cholesky_decomposition!(A)
This function is just for fun, and should not be taken seriously.
The lower triangle of matrix `A` will be replace by `L`, such that
the original `A = LL'`.

This is also algorithm 1 of
https://www.maths.manchester.ac.uk/~higham/papers/high09c.pdf
"""
function ordinary_cholesky_decomposition!(A)
    issymmetric(A) || error("Not a symmetric matrix")
    n = size(A)[1]
    for i in 1:n
        A[i, i] = sqrt(A[i, i] - A[i, 1:i-1]'A[i, 1:i-1])
        for j in i+1:n
            A[j, i] -= A[i, 1:i-1]'A[j, 1:i-1]
            A[j, i] /= A[i, i]
        end
    end
end

"""
    function cholesky_decomposition!(A)
https://www.maths.manchester.ac.uk/~higham/papers/high09c.pdf
algorithm 2.
It is also the same as algorithm 1 in
http://www.mucm.ac.uk/Pages/Downloads/Internal%20Reports/INT2.2.1%20LB%20Pivoting%20Cholesky%20Decomposition.pdf.
"""
function cholesky_decomposition!(A; tol=1e-5)
    issymmetric(A) || error("Not a symmetric matrix")
    n = size(A)[1]
    rank = n
    for i in 1:n
        A[i, i] < tol && (return A, i - 1)
        A[i, i] = sqrt(A[i, i])
        r = i+1:n               # range of rest block
        A[r, i] ./= A[i, i]
        A[r, r] -= A[r, i] * A[r, i]'
    end
    A, rank
end

"""
    function pivoted_cholesky_decomposition!(A; tol = 1e-5)
- Algorithm 1 of http://dfg-spp1324.de/download/preprints/preprint076.pdf.
- Algorithm https://mogp-emulator.readthedocs.io/en/latest/methods/proc/ProcPivotedCholesky.html

Return pivoted A, with factor in 'L', `piv`, and `rank`.
"""
function pivoted_cholesky_decomposition!(A; tol = 1e-5)
    issymmetric(A) || error("Not a symmetric matrix")
    n = size(A)[1]              # dimension
    p = collect(1:n)            # pivote
    #for i in 1:n
    for i in 1:1
        d = diag(A)[p[i:end]]   # diagonals
        l = argmax(d) + i - 1
        p[i], p[l] = p[l], p[i] # swap
        A[p[i], p[i]] < tol && return A[p, p], p, i-1
        A[p[i], p[i]] = sqrt(A[p[i], p[i]])
        r = p[i+1:n]
        A[r, p[i]] ./= A[p[i], p[i]]
        A[r, r] -= A[r, p[i]] * A[r, p[i]]'
    end
    A[p, p], p, n
end
