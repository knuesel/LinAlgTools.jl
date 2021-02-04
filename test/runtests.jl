using LinAlgTools
using Test

@testset "LinAlgTools.jl" begin
  A = [(i+j)//j for i in 1:3, j in 1:4]
  rref!(A)
  @test A == [1 0 -1//3 -1//2; 0 1 4//3 3//2; 0 0 0 0]
end
