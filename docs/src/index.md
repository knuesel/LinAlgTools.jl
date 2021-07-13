```@meta
CurrentModule = LinAlgTools
```

# LinAlgTools

Tools for teaching linear algebra. For now just a textbook implementation of
Gaussian elimination.

## Gaussian Elimination

The algorithm rarely works with matrices of integer type since they generally
cannot store the results of the multiplications by non-integers. However
rational arrays generally work well for small matrices:

```@repl gauss
A = Rational[ 4  2   0  -5   0   5
              4  3  -5  -2  -2   3
              2  2   2  -3   2   5
              5  4   3  -3  -2  -2 ]
```

Of course as the number of operations increases the fractions can become
unwieldy, and overflow unless e.g. `Rational{BigInt}` is used.

With SymPy.jl it is also possible to use symbolic matrices.

The functions [`ref`](@ref) and [`rref`](@ref) compute the row echelon form and
reduced row echelon form respectively. An optional `show_steps` argument can be
used to see each pass of the algorithm:

```@repl gauss
using LinAlgTools

B = ref(A, show_steps=true)
```

The variants [`ref!`](@ref) and [`rref!`](@ref) modify directly the matrix
given as parameter, and return the list of pivots, each as a `(row, column)`
tuple:

```@repl gauss
rref!(B)
B
```

Elementary row operations are available as [`row_swap`](@ref),
[`row_mul`](@ref) and [`row_add`](@ref):

```@repl gauss
row_swap(B, 1, 2)
row_mul(B, 1, by=10)
row_add(B, 1=>2, -1)
```

The variants [`row_swap!`](@ref), [`row_mul!`](@ref), [`row_add!`](@ref) modify
the given matrix. These functions return copies of the modified rows (a
potentially useful side-effect of keeping the implementation as simple as
possible).

```@repl gauss
row_swap!(B, 1, 2)
B
```
