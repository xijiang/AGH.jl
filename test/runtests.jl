using AGH
using Octavian, Test

@testset "Octavian validity" begin
    a = rand(2, 3)
    b = rand(3, 2)
    c = matmul(a, b)
    @test c == a * b
end
