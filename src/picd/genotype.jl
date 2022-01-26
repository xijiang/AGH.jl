# This file contains the function that manipulate the genotype files.
#=
Note:
=# 
"""
    function gheader(gt)
Retrieve the header infomation of a packed genotype file.
"""
function gheader(gt; target="GBLUP")
    # The header
    header = zeros(Int64, 3)
    read!(gt, header)
    locus_majored = (header[1] == 1)
    dms = Tuple(header[2:3])

    # Check header
    majored = ["locus", "ID"]
    i = locus_majored ? 1 : 2
    msg = ["Matrix storage:",
           "- This matrix is labeled $(majored[i])-majored.",
           "- That is, the genotypes are stored continuously on $(majored[i]).",
           "- Please use `AGH.PICD.transpose_G` to transpose it."]
    (locus_majored && (target == "GBLUP")) ||
        (!locus_majored && (target != "GBLUP")) ||
        @error join(msg, '\n')
    locus_majored, dms
end

"""
    function pack_G(src, out; g_majored = true, skip = 0)
Pack the genotypes in `src` into the packed genotype format.

## Source file
I presume that the source file `src` are in the following format. 
That is each line starts with an ID name, either of an animal or a SNP,
and then the genotypes, separated with space. And nothing else.
The number of header lines to be skipped is defined by `skip`.
The lines started with `#`s are also ignored. `skip` has priority.

## Packed genotypes
The binary file has a 24-byte header, which contains 3 Int64 number.
The first Int64 defines whether the genotypes are genotype-majored.
If `1`, then the genotypes of one ID are stored continuously
in the binary file `out`.  Otherwise, it is ID majored.  The genotypes
of one locus are stored continuously.  These 4 bytes can store other
information later.

The next 16-byte store the dimensions of the genotype matrix.
If `g_majored`, then it is `Nloci`×`Nid`.  Otherwise it is
`Nid`×`Nloci`.  For GBLUP, it is better to save `g-majored` to
speed up calculation.  For SNP-BLUP, a transposed storage is
more suitable.  The rest of the file store the genotypes of Int8.

Julia matrix storage is like FORTRAN, which is column majored.  That is,
if matrix `m = [1 2 3; 4 5 6]`, `write("t.bin", m)` will result in a
binary file that stores 1 4 2 5 3 6, consecutively of `Int64`.

If to read a line and write a line later, and also take the advantage of
sequential I/O, this must be considered very carefully.
"""
function pack_G(src, out; g_majored = true, skip = 0)
    m, n, gdm = 0, 0, zeros(Int64, 3)
    open(out, "w") do io
        write(io, gdm)
        open(src, "r") do gf
            [readline(gf) for _ in 1:skip];
            for line in eachline(src) # OBS!, this will transpose the matrix
                line[1] == '#' && continue
                s = parse.(Int8, split(line)[2:end])
                n += 1          # number of lines
                m += length(s)  # number of columns
                write(io, s)
            end
        end
    end            # genotypes written
    # A simple check point to see if the file is right.
    (m % n == 0) || @error "Genotypes not rectanglar"
    open(out, "r+") do io
        gdm = mmap(io, Vector{Int64}, 3)
        gdm[1], gdm[2], gdm[3] = g_majored, m/n, n # auto to Int64
        Mmap.sync!(gdm)
    end
end

"""
    function transpose_G(src, out)
---
Transpose **G** matrix in file `src`, and save the results to file `out`.
"""
function transpose_G(src, out)
    open(out, "w") do foo
        open(src, "r") do fin
            header = zeros(Int64, 3)
            read!(fin, header)
            header[1] = header[1] == 1 ? 0 : 1
            write(foo, header[1], header[3], header[2])
            gt = mmap(fin, Matrix{Int8}, Tuple(header[2:3]), 24)
            write(foo, gt')
        end
    end
end

"""
    function mat_W(pg, twop, foo)
Stream genotypes in `pg`, normalize them according `twop`,
and stream them to `foo`.

This funciton also uses `mmap` to deal very large dataset that
can't be fit in memory.
"""
function mat_W(pg, twop, foo)
    locus_majored, dms = gheader(pg)
    if locus_majored
        dms[1] == length(twop) || error("Wrong twop")
    else
        dms[2] == length(twop) || error("wrong twop")
    end
    gt = mmap(pg, Matrix{Int8}, dms, goffset)
    open(foo, "w") do io
        w = zeros(dms[1])
        header = zeros(Int64, 3)
        read!(pg, header)
        write(io, header)
        if locus_majored
            for i in 1:dms[2]
                copyto!(w, gt[:, i]) # always col by col, or cont. on disk.
                w -= twop
                write(io, w)
            end
        else
            for i in 1:dms[2]
                copyto!(w, gt[:, i])
                w .-= twop[i]
                write(io, w)
            end
        end
    end
end

"""
    function vr1_G(W, twop, G)
Given file `W` which contains the `W` matrix, which is the normalized G,
this function stream the G or van Raden method I to file `G`.
"""
function vr1_G(W, twop, G; step = 1000) # note: step is on ID for GBLUP
    s2pq = 1. / (1 .- .5twop)'twop
    locus_majored, dms = gheader(W) # read locus_majored.
    nlc, nid = dms
    w = mmap(W, Matrix{Float64}, dms, goffset)
    blk = zeros(nid, step)      # step rows of GBLUP G a time
    steps = collect(step:step:nid)
    nid % step == 0 || push!(steps, nid)
    open(G, "w") do io
        write(io, nid)          # header of dimension, one Int64 only.
        xf = 1                  # fra ID of wx
        for xt in steps         # end ID of wx
            xw = view(w, :, xf:xt)
            bt = xt - xf + 1    # number of rows in blk
            yf = 1              # fra ID of wy
            for yt in steps     # end ID of wy
                yw = view(w, :, yf:yt)
                xy = view(blk', 1:bt, yf:yt)
                matmul!(xy, xw', yw)
                yf += step
            end
            write(io, blk[:, 1:bt] .* s2pq) # written col by col
            xf += step
        end
    end
end
