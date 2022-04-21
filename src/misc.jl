"""
    function acgt201(v)
This function translate SNP ACGT alleles to 0/1.
It reports an error if more than two alternatives, and returns nothing.
It reports a warning if homozygous.
"""
function acgt201(v)
    d = Dict()
    for x in v
        d[x] = 1
    end
    if length(d) > 2
        @error "More than two alleles"
        return
    end
    length(d) == 1 && @warn "Homozygous"
    i = 0
    for (key, _) in d
        d[key] = i
        i += 1
    end
    z = zeros(Int8, length(v))
    for i in 1:length(v)
        z[i] = d[v[i]]
    end
    @debug "Dictionary" d
    z
end
