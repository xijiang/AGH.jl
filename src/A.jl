"""
    function AMat(ped)
---
"""
function AMat(ped)
end

function read_ped(file)
    @info "Reading pedigree from $file"
end

"""
    function code_ped()
---
"""
function code_ped()
end

function read_srt_ped(file)
    ped = Int32[]
    for line in eachline(file)
        (lstrip(line)[1] == '#') && continue
        append!(ped, parse.(Int32, split(line)[[1, 2]]))
    end
    reshape(ped, 2, :)'        
end

function AInv(ped)
end

"""
    function A_matrix(ped)
---
Give a matrix of 2 columns:
- Sire
- Dam

With row number as ID number, this function return a **A** matrix.
`Sire`s and `Dam`s are of integers less than their offspring `ID`.
This function cannot deal very large pedigree, i.e.,
with 32G momory, this function can only calculate a pedigree with
~90k ID.
With half precision, this can be ~130k.
With double precision, this function can only handle ~65k ID.
"""
function A_matrix(ped)
    @warn "Not right, needs debug"
    n = size(ped, 1)
    A = zeros(Float32, n, n)
    for i in 1:n
        ip, im = ped[i, :]      # ID i's pa and ma
        A[i, i] = 1.            # below avoid swap
        ip >0 && im >0 && (A[i, i] += .5(A[ip, im] + A[im, ip]))
        for j in i+1:n          # upper triangular
            jp, jm = ped[j, :]
            jp > 0 && (A[i, j] += .5A[i, jp])
            jm > 0 && (A[i, j] += .5A[i, jm])
            A[j, i] = A[i, j]
        end
        i % 100 == 0 && (print("\t$i"))
    end
    A
end

"""
    function kinship(ped, i, j)
---
This function is handy if just to calculate relationship of a few (pairs of) ID.
It can also speed up by adding `Thread.@threads` before your pair loop.
"""
function kinship(ped, i, j)
    (i == 0 || j == 0) && return 0
    ipa, ima = ped[i, :]          # used both for below and the last
    i == j && (return 1 + .5kinship(ped, ipa, ima))
    if i < j
        jpa, jma = ped[j, :]
        return .5(kinship(ped, i, jpa) + kinship(ped, i, jma))
    end
    return .5(kinship(ped, j, ipa) + kinship(ped, j, ima))
end

"""
    function kinship(ped, i, j, A)
---
This function is handy if just to calculate relationship of a few (pairs of) ID.
It can also speed up by adding `Thread.@threads` before your pair loop.
"""
function kinship(ped, i, j, A)
    (i == 0 || j == 0) && return 0
    
    !iszero(A[i, j]) && return A[i, j] # already calculated
    
    ip, im = ped[i, :]          # i's pa and ma
    i == j && begin
        A[i, i] = 1 + .5kinship(ped, ip, im, A)
        return A[i, i]
    end
    
    i < j && begin
        jp, jm = ped[j, :]      # j's pa and ma
        A[i, j] = .5(kinship(ped, i, jp, A) + kinship(ped, i, jm, A))
        return A[i, j]
    end

    A[i, j] = .5(kinship(ped, j, ip, A) + kinship(ped, j, im, A))
    return A[i, j]
end

"""
    function sprtA(ped, idx)
---
Part `A` of pedigree `ped` using sparse arrays.
List `idx` is the relationship ID to be calculated.

## Notes:
1. `idx` is better sorted.
2. `idx`'s biggest number should not be bigger than `size(ped, 1)`.
3. `ped` is of `Array{Int32, 2}`, with `pa` and `ma` column, and row number is `id` number.
"""
function sprtA(ped, idx)
    N = size(ped, 1)
    A = spzeros(N, N)
    n = length(idx)
    for i in n:-1:1
        x = idx[i]
        for j in i:-1:1
            y = idx[j]
            print("\r",
                  lpad("Row progress: $(n-i+1)/$n", 25),
                  lpad("ID: $x", 15),
                  ", ",
                  lpad("Col progress: $(i-j+1)/$i", 25),
                  lpad("ID: $y", 15))
            A[x, y] = A[y, x] = kinship(ped, x, y, A)
        end
    end
    println()
    Matrix(A[idx, idx])
end

#=
  case 'F':
  case 'f':
    for(auto id{1}; id<=nid; ++id)
      cout << Amat(id, id, ped, mid) - 1. <<'\n';
    break;
  case 'D':
  case 'd':
  case 'T':
  case 't':
    pma[0] = -1;

    cout << nid << '\n';

    for(auto id{1}; id<=nid; ++id){
      const auto&[pa, ma] = ped[id];
      if(pma.find(pa) == pma.end()) pma[pa] = Amat(pa, pa, ped, mid) - 1.;
      if(pma.find(ma) == pma.end()) pma[ma] = Amat(ma, ma, ped, mid) - 1.;
      cout << 1./(.5 - .25*(pma[pa] + pma[ma])) << '\n';
    }

    for(auto id{1}; id<=nid; ++id){
      const auto&[pa, ma] = ped[id];
      if(pa<ma) putT(pa, ma, id);
      else      putT(ma, pa, id);
    }
    break;
=#

function assignCSC(i, j, d, k, Ti, Tj, Td)
    k += 1
    Ti[k] = i
    Tj[k] = j
    Td[k] = d
    k
end

function A_inverse(ped)
    @info "Calculating A inversed"
    N = size(ped, 1)
    Ai = spzeros(N, N)
    A = spzeros(N, N)
    
    prt = collect(Set(ped))
    npm = length(prt)
    ele = zeros(npm)
    @showprogress "Inbreeding coefficients of parents" for i in 1:npm
        ele[i] = kinship(ped, prt[i], prt[i], A) - 1
    end

    pm = Dict(prt .=> ele)

    @info "Computing T & D"
    D = zeros(N)
    Ti = zeros(Int, 3N)
    Tj = zeros(Int, 3N)
    Td = zeros(3N)
    k = 0
    @showprogress "Matrix D and T" for i in 1:N
        pa, ma = ped[i, :]
        D[i] = 1. / (.5 - .25(pm[pa] + pm[ma]))
        pa > 0 && (k = assignCSC(i, pa, -.5, k,  Ti, Tj, Td))
        ma > 0 && (k = assignCSC(i, ma, -.5, k,  Ti, Tj, Td))
        k = assignCSC(i, i, 1., k,  Ti, Tj, Td)
    end
    T = sparse(Ti[1:k], Tj[1:k], Td[1:k])
    
    @info "Constructing A^-1"
    T'Diagonal(D)T
end
