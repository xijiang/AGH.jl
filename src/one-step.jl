"""
    function pcdos(file; threshold = 0.02)
---
`pcdos` means PCD method for one-step.

This function:
- construct a `G` matrix with the genotypes in `file`.
- does a pivoted Cholesky decomposition, say, `C`.
- take individuals whose corresponding `C.U` is greater than `threshold`.
- calculate `G` again and return ID used and `G` inverse.
"""
function pcdos(file; threshold = 0.02)
    @info "Reading genotypes"
    dat = read_g012(file)
    @info "Calculating GRM with VanRaden method I"
    G = grm(dat.genotype)
    @info "Pivoted Cholesky decomposition"
    C = cholesky(G, Val(true), check=false)
    @info "Find individuals whose correspoding diagonals in `C.U > $threshold`"
    i = 0
    for t in diag(C.U)
        t < threshold && break
        i += 1
    end
    id = sort(C.p[1:i])
    @info "Calculate reduced G and its inverse"
    G = grm(dat.genotype[:, id])
    return (id = dat.ID[id], gi = inv(G))
end
