```@meta
CurrentModule = LinAlgTools
```

# LinAlgTools

Tools for playing with linear algebra. For now just a textbook implementation
of Gaussian elimination.

## Gaussian Elimination

The algorithm rarely works with matrices of integer type since they generally
cannot store the results of the multiplications by non-integers. However
rational arrays work well:

```@repl gauss
A = Rational[ 4  2   0  -5   0   5
              4  3  -5  -2  -2   3
              2  2   2  -3   2   5
              5  4   3  -3  -2  -2 ]
```

Calls such as [`ref`](@ref) and [`rref`](@ref) return a modified copy of the
matrix, and the optional `show_steps` argument can be used to see each pass of
the algorithm:

```@repl gauss
using LinAlgTools

B = ref(A, show_steps=true)
```

Calls such as [`ref!`](@ref) and [`rref!`](@ref) modify the given matrix and
return a list of pivots, each as a `(row, column)` tuple:

```@repl gauss
rref!(B)
B
```

Elementary row operations are available as [`row_swap`](@ref),
[`row_mul`](@ref), [`row_add`](@ref), which return the modified matrix:

```@repl gauss
row_swap(B, 1, 2)
row_mul(B, 1, by=10)
row_add(B, 1=>2, -1)
```

The in-place variants [`row_swap!`](@ref), [`row_mul!`](@ref),
[`row_add!`](@ref) return copies of the modified rows (a potentially useful
side-effect of keeping the implementation as simple as possible).

```@repl gauss
row_swap!(B, 1, 2)
B
```
