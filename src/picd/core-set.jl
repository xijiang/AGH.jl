function core_set(gf, b; s = 1000, tol = 1e-5, tmp = ".", mem::Int64 = 8)
    n = readInt64(gf)           # dimension, first 8 bytes of G
    G = mmap(gf, Matrix{Float64}, (n, n), 8)
    d = diag(G)
    for i in 1:b
    end
    d
end
