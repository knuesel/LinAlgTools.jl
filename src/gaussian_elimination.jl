export row_swap!, row_mul!, row_add!, ref!, rref!

"""
    row_swap!(A, i, j)

Swap rows `i` and `j` in matrix `A`.

Returns the two rows in the new order.

# Example

```jldoctest
julia> A = reshape(1:12, 3, 4) |> collect
3×4 Array{Int64,2}:
 1  4  7  10
 2  5  8  11
 3  6  9  12

julia> row_swap!(A, 1, 3)
2×4 Array{Int64,2}:
 3  6  9  12
 1  4  7  10

julia> A
3×4 Array{Int64,2}:
 3  6  9  12
 2  5  8  11
 1  4  7  10
```
"""
row_swap!(A, i, j) = A[[i,j], :] = A[[j,i], :]

"""
    row_mul!(A, i, λ)

Multiply row `i` of matrix `A` by a factor `λ`.

Returns the new row `i`.

# Example

```jldoctest
julia> A = reshape(1:12, 3, 4) |> collect
3×4 Array{Int64,2}:
 1  4  7  10
 2  5  8  11
 3  6  9  12

julia> row_mul!(A, 2, 100)
4-element Array{Int64,1}:
  200
  500
  800
 1100

julia> A
3×4 Array{Int64,2}:
   1    4    7    10
 200  500  800  1100
   3    6    9    12
```
"""
row_mul!(A, i, λ) = A[i, :] *= λ

"""
    row_add!(A, i, λ, i₀)

Add to row `i` of matrix `A` row `i₀` multiplied by `λ`.

Returns the new row `i`.

# Example

```jldoctest
julia> A = [i//j for i in 1:3, j in 1:4]
3×4 Array{Rational{Int64},2}:
 1//1  1//2  1//3  1//4
 2//1  1//1  2//3  1//2
 3//1  3//2  1//1  3//4

julia> row_add!(A, 3, -3, 1)
4-element Array{Rational{Int64},1}:
 0//1
 0//1
 0//1
 0//1

julia> A
3×4 Array{Rational{Int64},2}:
 1//1  1//2  1//3  1//4
 2//1  1//1  2//3  1//2
 0//1  0//1  0//1  0//1
```
"""
row_add!(A, i, λ, i₀) =  A[i, :] += λ * A[i₀, :]

"""
    ref_pass!(A)

Perform one pass of the row-echelon-form algorithm.

Returns the column index of the pivot for this pass, or `nothing` if `A` is all
zeros.
"""
function ref_pass!(A)
  # Find first non-zero column
  col = findfirst(!iszero, collect(eachcol(A)))
  isnothing(col) && return nothing

  # Move pivot to first row if not there already
  if A[1, col] == 0
    row = findfirst(!iszero, A[:, col])
    row_swap!(A, 1, row)
  end

  # Use row operations to bring zeros under the pivot
  for row in 2:size(A, 1)
    if A[row, col] != 0
      row_add!(A, row, -A[row, col]/A[1, col], 1)
    end
  end

  return col
end

"""
    log(A, txt)

Show matrix `A` with title text `txt`.
"""
function log(A, txt)
  println(txt)
  println(repeat('-', length(txt)))
  display(A)
  println()
end

"""
    ref!(A; show_steps=false)

Put A in row echelon form.

Returns the (row, column) indices of the pivots.

# Example

```jldoctest
julia> A = Float64[mod(i+j,3) for i in 1:4, j in 1:6]
4×6 Array{Float64,2}:
 2.0  0.0  1.0  2.0  0.0  1.0
 0.0  1.0  2.0  0.0  1.0  2.0
 1.0  2.0  0.0  1.0  2.0  0.0
 2.0  0.0  1.0  2.0  0.0  1.0

julia> ref!(A)
3-element Array{Tuple{Int64,Int64},1}:
 (1, 1)
 (2, 2)
 (3, 3)

julia> A
4×6 Array{Float64,2}:
 2.0  0.0   1.0  2.0  0.0   1.0
 0.0  1.0   2.0  0.0  1.0   2.0
 0.0  0.0  -4.5  0.0  0.0  -4.5
 0.0  0.0   0.0  0.0  0.0   0.0
```
"""
function ref!(A; show_steps=false)
  show_steps && log(A, "Start")
  pivots = Tuple{Int,Int}[]

  # Iterate on the rows, each time doing one pass and "hiding" another row.
  for startrow = 1:size(A, 1)
    col = ref_pass!(@view A[startrow:end, :])

    # No pivot found -> nothing to do
    col == nothing && break

    # Store the position of the pivot
    push!(pivots, (startrow, col))
    show_steps && log(A, "Pivot $(pivots[end])")
  end
  return pivots
end

"""
    rref!(A; show_steps)

Compute the reduced row echelon form of A.

Returns the (row, column) indices of the pivots.

# Example

```jldoctest
julia> A = Float64[mod(i+j,3) for i in 1:4, j in 1:6]
4×6 Array{Float64,2}:
 2.0  0.0  1.0  2.0  0.0  1.0
 0.0  1.0  2.0  0.0  1.0  2.0
 1.0  2.0  0.0  1.0  2.0  0.0
 2.0  0.0  1.0  2.0  0.0  1.0

julia> rref!(A)
3-element Array{Tuple{Int64,Int64},1}:
 (1, 1)
 (2, 2)
 (3, 3)

julia> A
4×6 Array{Float64,2}:
  1.0   0.0  0.0   1.0   0.0  0.0
  0.0   1.0  0.0   0.0   1.0  0.0
 -0.0  -0.0  1.0  -0.0  -0.0  1.0
  0.0   0.0  0.0   0.0   0.0  0.0
```
"""
function rref!(A; show_steps=false)
  # First calculate non-reduced form
  pivots = ref!(A; show_steps)

  # Backward phase: iterate from last to first pivot
  for (row, col) in reverse(pivots)
    # Multiply row to have pivot = 1
    pivot = A[row, col]
    row_mul!(A, row, inv(pivot))
    show_steps && log(A, "Reduce $((row,col)) $pivot -> 1")

    # Ensure all zeros above the pivot
    for i in 1:row-1
      row_add!(A, i, -A[i, col], row)
    end
    show_steps && log(A, "Reduce above/below $((row,col))")
  end
  return pivots
end
