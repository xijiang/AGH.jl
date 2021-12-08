using SparseArrays, LinearAlgebra, ProgressMeter, Serialization

function spcT()
    @info "Reading the data and genotyped ID"
    Ai = deserialize("Ainv.ser")
    jx = begin
        tmp = Int[]
        for i in eachline("genotyped_renumid_sort.dat")
            push!(tmp, parse(Int, i))
        end
        tmp
    end
    N = size(Ai, 1)
    n2 = length(jx)
    ix = sort(collect(setdiff(Set(1:N), Set(jx))))
    A11 = Ai[ix, ix]
    A12 = Ai[ix, jx]
    A22 = Matrix(Ai[jx, jx])
    nr = n2 ÷ 10
    stops = collect(nr:nr:n2)
    stops[end] = n2
    start = 1
    @info "Calculating by blocks"
    @showprogress for stop in stops
        rhs = Matrix(A12[:, start:stop])
        sol = A11 \ rhs
        blk = A12'sol
        A22[:, start:stop] -= blk
        start = stop + 1
    end
    serialize("a22.txt", A22)
end
