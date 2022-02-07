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
    function icdpob(gf, b; s = 1000, tol = 1e-5, tmp = ".", ram = 8)
This is an approximation method of Cholesky decomposition.
The function name is an acroname of `i`ncomplete `C`holesky `d`ecomposition
`pivoted` `o`n `b`locks.
It starts with a ramdom set of `G`.  Then select the next block which has
the largest diagonals for Cholesky decomposition procedure.
This is finished until the required blocks are done.

## Arguments
- `gf`: filename storing the big `G` matrix.
    Its first 8 bytes stores the binary dimension,
    follows the `G` elements column by column in binary of Float64.
- `b`: number of blocks you want to approximate `G`.
- `s`: block size, 1000 columns by default.
- `tol`: decomposition stops prematurely if the conditional diagonal
    is smaller than `tol`.
- `tmp`: the temporary directory (big enough) to store intermediate results.
    The program may re-write updated sub-`G` several times on disk.
- `ram`: available RAM to deal with lower right big `G`, in gigabytes.

## Notes
- The function returns `Z` of `n × bs`, `ZZ' ≈ G`.
- The function starts from a random set of the `G` columns.  
    Hence, it gives different results on different runs.
"""
function icdpob(gf, b; s = 1000, tol = 1e-5, tmp = ".", ram::Int64 = 8)
    n = readInt64(gf)           # dimension, first 8 bytes of G
    bs = b * s;    bs > n && error("You're asking too many blocks")

    @debug "ToDo: Enough memory (G)?" n*bs*8/1024^3
    @debug "ToDo: starting block not random ?"
    Z, piv, rest, file = zeros(n, bs), Int[], collect(1:n), gf
    for i in 1:b
        G = mmap(file, Matrix{Float64}, (n, n), 8)
        d = diag(G)
        i == 1 ? p = shuffle(1:n) : p = sortperm(d, rev = true)
        #p = sortperm(d, rev = true)
        x, y = view(p, 1:s), view(p, s+1:n);    sort!(x);    sort!(y)
        i > 1 && pivote_pz!(Z, i, s, p)
        z, g = begin
            l, r = (i - 1) * s + 1, size(Z)[1]
            view(Z, l:r, l:i*s),   view(G, y, y)
        end
        @debug "testing" g
        copyto!(z, G[p, x]) # columns to z
        append!(piv, rest[x]);  rest = rest[y]
        info = update_z!(z, tol)
        l = (i - 1) * s + info
        info > 0 && return Z[:, 1:l], piv[1:l]
        # update g by blocks and write to tmp/tmpfile
        tf = tempname(tmp)
        n -= s
        i == b && return Z, piv
        open(tf, "w") do io
            write(io, n)
            cs = cslrg(n, ram)   # number of g columns to deal a time
            f = 1
            for k in cs
                update_n_write_lrg(io, g, f, k, z[s+1:end, s])
                f = k + 1
            end
        end
        i > 2 && rm(file)
        file = tf
    end
end
