# This file contains the function that manipulate the genotype files.
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

Note that Julia matrix is column majored.
"""
function pack_G(src, out; g_majored = true, skip = 0)
    m, n, gdm = 0, 0, zeros(Int64, 3)
    open(out, "w") do io
        write(io, gdm)
        open(src, "r") do gf
            [readline(gf) for _ in 1:skip]
            for line in eachline(src)
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
            write(foo, 1 - header[1], header[3], header[2])
            gt = mmap(fin, Matrix{Int8}, Tuple(header[2:3]), 24)
            write(foo, gt')
        end
    end
end
