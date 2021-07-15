export row_swap!, row_mul!, row_add!, ref!, rref!,
       row_swap,  row_mul,  row_add,  ref,  rref

# Call f on a copy of A and return the (modified) copy
function with_copy(f, A, args...; kwargs...)
  B = copy(A)
  f(B, args...; kwargs...)
  return B
end

"""
    row_swap!(A, i, j)

Swap rows `i` and `j` in matrix `A`.

Returns a copy of the two rows in the new order.

# Example

```jldoctest
julia> A = reshape(1:12, 3, 4) |> collect
3×4 Matrix{Int64}:
 1  4  7  10
 2  5  8  11
 3  6  9  12

julia> row_swap!(A, 1, 3)
2×4 Matrix{Int64}:
 3  6  9  12
 1  4  7  10

julia> A
3×4 Matrix{Int64}:
 3  6  9  12
 2  5  8  11
 1  4  7  10
```
"""
row_swap!(A, i, j) = A[[i,j], :] = A[[j,i], :]

"""
    row_swap(A, i, j)

Create a matrix equal to `A` but with rows `i` and `j` swapped.

# Example

```jldoctest
julia> A = reshape(1:12, 3, 4) |> collect
3×4 Matrix{Int64}:
 1  4  7  10
 2  5  8  11
 3  6  9  12

julia> row_swap(A, 1, 3)
3×4 Matrix{Int64}:
 3  6  9  12
 2  5  8  11
 1  4  7  10
```
"""
row_swap(A, i, j) = with_copy(row_swap!, A, i, j)

"""
    row_mul!(A, i => λ)

Multiply row `i` of matrix `A` by a factor `λ`.

Returns a copy of the new row `i`.

# Example

```jldoctest
julia> A = reshape(1:12, 3, 4) |> collect
3×4 Matrix{Int64}:
 1  4  7  10
 2  5  8  11
 3  6  9  12

julia> row_mul!(A, 2 => 100)
4-element Vector{Int64}:
  200
  500
  800
 1100

julia> A
3×4 Matrix{Int64}:
   1    4    7    10
 200  500  800  1100
   3    6    9    12
```
"""
row_mul!(A, (i, factor)) = A[i, :] *= factor

"""
    row_mul(A, i => factor)

Create a matrix equal to `A` but with row `i` multiplied by `λ`.

# Example

```jldoctest
julia> A = reshape(1:12, 3, 4) |> collect
3×4 Matrix{Int64}:
 1  4  7  10
 2  5  8  11
 3  6  9  12

julia> row_mul(A, 2 => 100)
3×4 Matrix{Int64}:
   1    4    7    10
 200  500  800  1100
   3    6    9    12
```
"""
row_mul(A, (i, factor)) = with_copy(row_mul!, A, i => factor)

"""
    row_add!(A, i => λ => j)

Add row `i` multiplied by `λ` to row `j` of matrix `A`.

Returns a copy of the new row `j`.

# Example

```jldoctest
julia> A = [i//j for i in 1:3, j in 1:4]
3×4 Matrix{Rational{Int64}}:
 1//1  1//2  1//3  1//4
 2//1  1//1  2//3  1//2
 3//1  3//2  1//1  3//4

julia> row_add!(A, 1 => -3 => 3)
4-element Vector{Rational{Int64}}:
 0//1
 0//1
 0//1
 0//1

julia> A
3×4 Matrix{Rational{Int64}}:
 1//1  1//2  1//3  1//4
 2//1  1//1  2//3  1//2
 0//1  0//1  0//1  0//1
```
"""
row_add!(A, (i, (λ, j))) = A[j, :] += λ * A[i, :]

"""
    row_add(A, i => λ => j)

Create a matrix equal to `A` but with `λ` times row `i` added to row `j`.

# Example

```jldoctest
julia> A = [i//j for i in 1:3, j in 1:4]
3×4 Matrix{Rational{Int64}}:
 1//1  1//2  1//3  1//4
 2//1  1//1  2//3  1//2
 3//1  3//2  1//1  3//4

julia> row_add(A, 1 => -3 => 3)
3×4 Matrix{Rational{Int64}}:
 1//1  1//2  1//3  1//4
 2//1  1//1  2//3  1//2
 0//1  0//1  0//1  0//1
```
"""
row_add(A, (i, (λ, j))) = with_copy(row_add!, A, i => λ => j)

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
    row = findfirst(!iszero, @view A[:, col])
    row_swap!(A, 1, row)
  end

  # Use row operations to bring zeros under the pivot
  for row in 2:size(A, 1)
    if A[row, col] != 0
      row_add!(A, 1 => -A[row, col]/A[1, col] => row)
    end
  end

  return col
end

"""
    log(A, txt)

Write matrix `A` with title text `txt` to standard output.
"""
function log(A, txt)
  println(txt)
  println(repeat('-', length(txt)))
  show(stdout, MIME("text/plain"), A)
  println()
  println()
end

"""
    ref!(A; show_steps=false)

Put `A` in row echelon form.

Returns the `(row, column)` indices of the pivots.

# Example

```jldoctest
julia> A = Float64[mod(i+j,3) for i in 1:4, j in 1:6]
4×6 Matrix{Float64}:
 2.0  0.0  1.0  2.0  0.0  1.0
 0.0  1.0  2.0  0.0  1.0  2.0
 1.0  2.0  0.0  1.0  2.0  0.0
 2.0  0.0  1.0  2.0  0.0  1.0

julia> ref!(A)
3-element Vector{Tuple{Int64, Int64}}:
 (1, 1)
 (2, 2)
 (3, 3)

julia> A
4×6 Matrix{Float64}:
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
    ref(A; show_steps=false)

Compute a row echelon form of `A`.

Returns the `(row, column)` indices of the pivots.

# Example

```jldoctest
julia> A = Float64[mod(i+j,3) for i in 1:4, j in 1:6]
4×6 Matrix{Float64}:
 2.0  0.0  1.0  2.0  0.0  1.0
 0.0  1.0  2.0  0.0  1.0  2.0
 1.0  2.0  0.0  1.0  2.0  0.0
 2.0  0.0  1.0  2.0  0.0  1.0

julia> ref(A)
4×6 Matrix{Float64}:
 2.0  0.0   1.0  2.0  0.0   1.0
 0.0  1.0   2.0  0.0  1.0   2.0
 0.0  0.0  -4.5  0.0  0.0  -4.5
 0.0  0.0   0.0  0.0  0.0   0.0
```
"""
ref(A; show_steps=false) = with_copy(ref!, A; show_steps)

"""
    rref!(A; show_steps)

Put `A` in the reduced row echelon form.

Returns the `(row, column)` indices of the pivots.

# Example

```jldoctest
julia> A = Float64[mod(i+j,3) for i in 1:4, j in 1:6]
4×6 Matrix{Float64}:
 2.0  0.0  1.0  2.0  0.0  1.0
 0.0  1.0  2.0  0.0  1.0  2.0
 1.0  2.0  0.0  1.0  2.0  0.0
 2.0  0.0  1.0  2.0  0.0  1.0

julia> rref!(A)
3-element Vector{Tuple{Int64, Int64}}:
 (1, 1)
 (2, 2)
 (3, 3)

julia> A
4×6 Matrix{Float64}:
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
    if pivot != 1
        row_mul!(A, row => inv(pivot))
        show_steps && log(A, "Reduce $((row,col)) $pivot -> 1")
    end

    # Ensure all zeros above the pivot
    if !all(iszero, A[1:row-1, col])
        for i in 1:row-1
            row_add!(A, row => -A[i, col] => i)
        end
        show_steps && log(A, "Reduce above $((row,col))")
    end
  end
  return pivots
end

"""
    rref(A; show_steps)

Compute the reduced row echelon form of `A`.

# Example

```jldoctest
julia> A = Float64[mod(i+j,3) for i in 1:4, j in 1:6]
4×6 Matrix{Float64}:
 2.0  0.0  1.0  2.0  0.0  1.0
 0.0  1.0  2.0  0.0  1.0  2.0
 1.0  2.0  0.0  1.0  2.0  0.0
 2.0  0.0  1.0  2.0  0.0  1.0

julia> rref(A)
4×6 Matrix{Float64}:
  1.0   0.0  0.0   1.0   0.0  0.0
  0.0   1.0  0.0   0.0   1.0  0.0
 -0.0  -0.0  1.0  -0.0  -0.0  1.0
  0.0   0.0  0.0   0.0   0.0  0.0
```
"""
rref(A; show_steps=false) = with_copy(rref!, A; show_steps)
