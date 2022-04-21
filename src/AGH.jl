module AGH
using LinearAlgebra, SparseArrays, ProgressMeter, Statistics, LinearAlgebra,
    Random, Plots, Octavian

include("A.jl")
include("G.jl")
include("H.jl")
include("misc.jl")
include("file-io.jl")
include("one-step.jl")
include("picd/PICD.jl")

end # module
