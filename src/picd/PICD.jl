"""
    module PICD
---
The module that does `p`ivoted `i`ncomplete `C`holesky `d`ecomposition(PICD).
"""
module PICD
using Mmap, Statistics, Octavian, LinearAlgebra, Random

include("genotype.jl")
include("cholesky-decomposition.jl")
include("icdpob.jl")
include("core-set.jl")

const goffset = 24

end
