"""
    function read_g012(file)
---
This function read a `plink 012` genotype `file`.
It return a named tuple, with an `ID` string vector,
and `genotype` Int8 matrix.

Note, this function doesn't do error check.
The file integrity should be checked before this function.

The genotype matrix is of `nLoci` by `nID`, i.e., ID column majored.
Genotypes of the same ID are in the same column.
"""
function read_g012(file)
    # Determine the dimensions of the genotypes
    nid = countlines(file)
    nlc = length(split(readline(file))) - 1
    id = String[]
    genotype = zeros(Int8, nlc, nid)
    i = 1
    for line in eachline(file)
        # ToDo: speed up this reading using buffer.
        t = split(line)
        push!(id, t[1])
        genotype[:, i] = parse.(Int8, t[2:end])
        i += 1
    end
    return (ID = id, genotype = genotype)
end
