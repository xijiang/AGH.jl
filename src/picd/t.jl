
function pack_SNP(pg)
    locus_majored, dms = gheader(pg)
    nlc, nid = dms
    haps = BitArray(undef, nlc, 2nid)
    g = Mmap.mmap(pg, Matrix{Int8}, dms, 24)
    nt = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    Threads.@threads for i in 1:nid
        x, y = zeros(Int8, nlc), zeros(Int8, nlc)
        for j in 1:nlc
            if g[j, i] == 0
                x[j], y[j] = 0, 0
            elseif g[j, i] == 1
                x[j], y[j] = 0, 1
            else
                x[j], y[j] = 1, 1
            end
        end
        haps[:, 2i-1] = x
        haps[:, 2i] = y
    end
    BLAS.set_num_threads(nt)
    haps
end

        
function diag_bits(haps, twop)
    s2pq = 1. / (1 .- .5twop)'twop
    nlc, nhp = size(haps)
    nid = nhp ÷ 2
    d = zeros(nid)
    nt = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    Threads.@threads for i in 1:nid
        c = haps[:, 2i-1] + haps[:, 2i]
        c -= twop
        d[i] = c'c * s2pq
    end
    BLAS.set_num_threads(nt)
    d
end


function part_G_old(haps, ids, twop)
    c = 1. / (1 .- .5twop)'twop
    s = twop'twop
    n = length(ids)
    g = zeros(n, n)
    nt = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    for i in 1:n
        x = view(haps, :, 2ids[i]-1:2ids[i])
        Threads.@threads for j in 1:i
            y = view(haps, :, 2ids[j]-1:2ids[j])
            g[i, j] = sum(x'y)
            g[i, j] -= (sum([x y], dims=2)' * twop)[1]
            g[i, j] += s
            g[i, j] *= c
        end
    end
    BLAS.set_num_threads(nt)
    g
end

function part_G_v2(haps, ids, twop)
    c = 1. / (1 .- .5twop)'twop
    s = twop'twop
    n = length(ids)
    g = zeros(n, n)
    g2p = zeros(n)
    for i in 1:n
        x = view(haps, :, 2ids[i]-1:2ids[i])
        g2p[i] = matmul(sum(x, dims=2)', twop)[1]
    end
    nt = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    for i in 1:n
        x = view(haps, :, 2ids[i]-1:2ids[i])
        a = sum(x, dims = 2)
        #g[i, i] = (s + matmul(a', a)[1] - 2g2p[i]) * c
        Threads.@threads for j in 1:i
            y = view(haps, :, 2ids[j]-1:2ids[j])
            b = sum(y, dims = 2)
            g[i, j] = (s + matmul(a', b)[1] - g2p[i] - g2p[j]) * c
        end
    end
    BLAS.set_num_threads(nt)
    g
end
