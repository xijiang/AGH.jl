"""
    module PICD
---
The module that does `p`ivoted `i`ncomplete `C`holesky `d`ecomposition(PICD).
"""
module PICD
using Mmap, Statistics, Octavian, LinearAlgebra

include("genotype.jl")
include("picd.jl")

function picd_H()
end

end
