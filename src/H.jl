"""
---
A^{-1}_{22} = A^{22} - A^{21}(A^{11}^{-1})A^{12}
"""
function A_22inv(Ai, n1, nc)
    n2 = size(Ai, 1) - n1
    A22 = Ai[n1+1:end, n1+1:end]
    A12 = Ai[1:n1, n1+1:end]
    A11 = Ai[1:n1, 1:n1]

    v = A11 \ vec(A12[:, 1])
    A12'v
end

