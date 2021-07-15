using LinAlgTools
using Test

@testset "LinAlgTools.jl" begin
  @testset "In-place calls" begin
    A = [(i+j)//j for i in 1:3, j in 1:4]
    rref!(A)
    @test A == [1 0 -1//3 -1//2; 0 1 4//3 3//2; 0 0 0 0]

    newrows = row_swap!(A, 1, 3)
    @test A == [0 0 0 0; 0 1 4//3 3//2; 1 0 -1//3 -1//2]
    @test newrows == [0 0 0 0; 1 0 -1//3 -1//2]
    @test !Base.mightalias(newrows, A)

    mulrow = row_mul!(A, 3 => 10)
    @test A == [0 0 0 0; 0 1 4//3 3//2; 10 0 -10//3 -5]
    @test mulrow == [10, 0, -10//3, -5]
    @test !Base.mightalias(mulrow, A)

    addrow = row_add!(A, 3 => 10 => 2)
    @test A == [0 0 0 0; 100 1 -32 -97//2; 10 0 -10//3 -5]
    @test addrow == [100, 1, -32, -97//2]
    @test !Base.mightalias(addrow, A)
  end

  @testset "Non-mutating calls" begin
    A = [(i+j)//j for i in 1:3, j in 1:4]
    Aref = ref(A)
    @test Aref ==  [2 3//2 4//3 5//4; 0 -1//4 -1//3 -3//8; 0 0 0 0]
    @test !Base.mightalias(Aref, A)

    Arref = rref(A)
    @test Arref == [1 0 -1//3 -1//2; 0 1 4//3 3//2; 0 0 0 0]
    @test !Base.mightalias(Arref, A)

    Aswap = row_swap(A, 1, 3)
    @test Aswap == [4 5//2 2 7//4; 3 2 5//3 3//2; 2 3//2 4//3 5//4]
    @test !Base.mightalias(Aswap, A)

    Amul = row_mul(A, 3 => 10)
    @test Amul == [2 3//2 4//3 5//4; 3 2 5//3 3//2; 40 25 20 35//2]
    @test !Base.mightalias(Amul, A)

    Aadd = row_add(A, 3 => 10 => 2)
    @test Aadd == [2 3//2 4//3 5//4; 43 27 65//3 19; 4 5//2 2 7//4]
    @test !Base.mightalias(Aadd, A)
  end
end
