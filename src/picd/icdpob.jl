"""
    function readInt(file)
This small function read the first 64-bit integer of a file, 
which is usually the dimesion of a sqaure matrix.
"""
function readInt64(file)
    t = zeros(Int64, 1)
    read!(file, t)
    t[1]
end

"""
    function update_z!(z, tol)
This is to avoid ugly many-line function.
`z` is of dimension `r` × `s`, where `r` ≥ `s`.
A non-zero return indicates premature updating.
""" 
function update_z!(z, tol)
    r, s = size(z)
    r < s && @error("Not enough rows in z")
    for i in 1:s
        z[i, i] < tol && return i - 1
        z[i, i] = sqrt(z[i, i])
        z[i+1:r, i] ./= z[i, i]
        Threads.@threads for j in i+1:r
            for k in i+1:s
                z[j, k] -= z[j, i] * z[k, i]
            end
        end
    end
    return 0
end

"""
    function cslrg(n::Int64, mem::Int64)
Determine `c`olumn `s`teps of `l`ower `r`ight `G` to deal with.
"""
function cslrg(n::Int64, mem::Int64)
    s = Int(floor(mem * 1024^3 / 8 / n))
    cs = collect(s:s:n)
    n % s == 0 || push!(cs, n)
    cs
end

"""
    function update_n_write_lrg(io, g, f, k, z)
Update lower right G, `g[:, f:k]` with `z` and append to IOStream `io`.
As the update is on a copy of sub `g`, there is no `!` in the function name.
"""
function update_n_write_lrg(io, g, f, k, z)
    t = view(g, :, f:k)
    v = similar(t)
    copyto!(v, t)
    r, c = size(v)
    w = view(z, f:k)
    Threads.@threads for i in 1:r
        for j in 1:c
            v[i, j] -= z[i] * w[j]
        end
    end
    @debug "testing" v

    write(io, v)
end

"""
    function pivote_pz(Z, i, s, p)
When `i` > 1, the lower right G are pivoted on their diagonals.
The previous calculated `Z`, or `p`artial `z`, 
should pivote correspondingly.
"""
function pivote_pz!(Z, i, s, p)
    r, c = size(Z)
    l = length(p)
    z = view(Z, r - l + 1:r, 1:(i-1)*s)
    z[:, :] = z[p, :]
end

"""
    function apprxG(gf, nc; δ = .01, tol = 1e-5, ram=16)
This function is to find a optimized core set of individuals to approximate
G matrix.
- `gf`: Filename for the G matrix. Its first 8 bytes contains the dimension, then
  elements of Float64 column by column.
- `nc`: Number of core animals to be included.
- `δ`: A small value to be added to the diagonal of G.
- `ram`: Available memory for this calculation.
"""
function apprxG(gf, nc; δ = .01, tol = 1e-5, ram=16)
    n = readInt64(gf)
    G = mmap(gf, Matrix{Float64}, (n, n), 8)
    D = diag(G)
    p = sortperm(D)
    #Z = copy(G)
    #Z += δ * I
    #Z, piv, rank, info = LAPACK.pstrf!('L', Z, tol)
end

function mypotri!(M)
    LAPACK.potri!('L', M)
    n, _ = size(M)
    for i in 1:n
        for j in i+1:n
            M[i, j] = M[j, i]
        end
    end
end

"""
    function apprxGi(gf, nc, foo; δ = 0.01, tol = 1e-5, ram=16)
This is now a dirty but 'theoretically' best solution for approx. G.
"""
function apprxGi(gf, nc, foo; δ = 0.01, tol = 1e-5, ram=16)
    n = readInt64(gf)
    G = mmap(gf, Matrix{Float64}, (n, n), 8)
    Z = copy(G)
    Z += δ * I
    Z, piv, rk, info = LAPACK.pstrf!('L', Z, tol)
    b1 = piv[1:nc]
    b2 = piv[nc+1:end]
    G11i = copy(Z[1:nc, 1:nc])
    mypotr!i(G11i)
    G21 = copy(G[b2, b1])
    G12 = copy(G[b1, b2])
    Di = diagm(1 ./ diag(Z)[b2])
    #inv11 = G11i + G11i * G21 * Di
end

function tease(; tol = 1e-5)
    M = begin
        t = rand(10, 5)
        t't
    end
    L = copy(M)
    @debug "Matrix: " M
    _, piv, rk, info = LAPACK.pstrf!('L', L, tol)
    for i in 1:5
        for j in i+1:5
            L[i, j] = 0
        end
    end
    @debug "Factor: " L
    r = piv[1:2]
    X = copy(M[r, r])
    LAPACK.pstrf!('L', X, tol)
    @debug "Factor: " X
end
    
