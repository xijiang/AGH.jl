"""
    function grm(gt; nz=0.)
---
This is an instance of Van Raden method I.
Genotypes are of value 0, 1 or 2, and of ID column majored.
They should be stored per ID per column.
This function uses vanRaden method I.
I do no genotype and frequency check here,
As it should be done somewhere else.
"""
function grm(g; nz=0.)
    twop = mean(g, dims = 2)
    s2pq = (1 .- .5 .* twop)'twop # one element array
    r2pq = 1. / s2pq[1]
    Z = g .- twop
    G = Z'Z .* r2pq
    !iszero(nz) && (G += nz .* I)
    G
end
