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

"""
    function mkGroup(d, ng)
Give the diagonals `d`, this function find the index of the largest element.
It then group the indices of `d` into `ng` groups.
Each group will contains the index of the largest element.
"""
function mkGroup(d, ng)
    m, n = findmax(d)[2], length(d)
    p = randperm(n)
    s = n ÷ ng
    ss = collect(s:s:n)
    ss[end] = n
    b = 1
    grp = []
    for s in ss
        t = view(p, b:s)
        findnext(x -> x == m, t, 1) == nothing ?
            push!(grp, sort([t; m])) :
            push!(grp, sort(t))
        b = s + 1
    end
    grp              
end

function findbsub(gf, nc; ng = 3, δ = .01, tol = 1e-5)
    n = readInt64(gf)
    G = mmap(gf, Matrix{Float64}, (n, n), 8)
    d = diag(G)
    #ng = Int(ceil(n/nc))
    r = Int(ceil(nc / ng * 1.8))
    @debug "ratio" r
    grp = mkGroup(d, ng)
    rst = Set{Int64}[]
    for x in grp
        m = length(x)
        g = zeros(m, m)
        copyto!(g, G[x, x])
        _, piv, rk, info = LAPACK.pstrf!('L', g, tol)
        rst = rst ∪ Set(x[piv[1:r]])
        @debug "collected" length(rst)
    end
    rst = sort(collect(rst))
    m = length(rst)
    g = zeros(m, m)
    copyto!(g, G[rst, rst])
    _, piv, rk, info = LAPACK.pstrf!('L', g, tol)
    #sort(piv)[1:nc]
    piv
end

function teaser(op; ng = 3, tol = 1e-5)
    file = "dat/picd/G.bin"
    n = readInt64(file)
    G = mmap(file, Matrix{Float64}, (n, n), 8)
    δ = .01
    
    if op == 1                  # test 1, a full PCD
        g = zeros(n, n)
        copyto!(g, G)
        g += δ * I
        g, piv, rk, info = LAPACK.pstrf!('L', g, tol)
    elseif op == 2 # test 2: partial PCD, with the largest diagonal included in 3 set
        Set{Int64}
    end
end
