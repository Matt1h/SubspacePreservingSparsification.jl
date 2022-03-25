```@meta
EditURL = "<unknown>/docs/literate/example.jl"
```

# example
Lets try to find a sparse representation for a matrix, so that the Frobeniusnorm of the difference
is as small as possible.
Lets consider a symmetric matrix $M$:

````@example example
using SubspacePreservingSparsification
using SparseArrays
using LinearAlgebra
M = [17.05 16.98 0.3 6.99 7; 16.98 0.2 7.1 6.9 0; 0.3 7.1 -12 0.01 17; 6.99 6.9 0.01 -11.97 0; 7 0 17 0 -0.1]
````

First we calculate a sparsity pattern:

````@example example
M_id = sparsity_pattern(M, 0.6, 2)
````

For every column and every row of the Matrix the function solves the optimization problem
```math
max_{Z}\, ||x - x\circ Z(x)||_{0}\,\, s.t. \\
```
```math
\begin{align*}
(Z(x))_{i} &= 0 \, \text{or}\, 1 \\
||Z(x)||_{0} &\ge N \\
||x - x\circ Z(x)||_{p} 	\,&\le\,  (1-q)||x||_{p},
\end{align*}
```
where $Z$ is the sparsity pattern for the column or row. $N$ is a maximum number for non zero entries,
$q$ is a factor that controls how sparse the the row or column should be. All the individual row and column
patter are overlayed, if either row or column or both return a one for an entry the entry is one.
This ensures that the pattern will preserve different Subspaces of the Matrix.

Than we can take the sparsity pattern and modify it to a binning pattern:

````@example example
bin_sparse_matrix!(M, M_id, 200)
````

The function looks at the maximum and minimum values of the matrix and then finds a partition
of the negative and positive area according to the specified maximum number of bins.
For each entry of the matrix it is checked to which bin it belongs to. The binning controls the
number of unknowns and so the computation cost to solve the optimization problem. It also typically
improves the conditioning of the optimization problem.

With the binning pattern we can find the Matrix so that the norm of the difference
is as small as possible and so that the binning constrains are fullfilled:

````@example example
binned_minimization(M_id, zeros(5, 5), Matrix{Float64}(I, 5, 5), M)
````

The function solves the optimization problem
```math
min ||M - X||_{F} \,\, s.t.
```
```math
X\, \text{has the specified sparsity pattern}\, B(M) = M_{id}
```
Of course our example optimization problem is not too interesting, because without the binning
only some entries are set to zero, while the others remain the same. However,
more sophisticated optimization problems can be used here.

We can also use the sparsify function on $M$:

````@example example
sparsify(M, 0.6, 2, 200)
````

This function calculates the sparsity and binning pattern like we did and uses as optimization problem
```math
min\, \frac{1}{2} \sum\limits_{i=1}^r \frac{1}{\sigma_i^2} ||\left( X - M\right) v_i||_2^2
+ \frac{1}{2} \sum\limits_{i=1}^r \frac{1}{\sigma_i^2} ||\left( X^* - M^*\right) u_i||_2^2,
```
where r is the rank, $\sigma_i$ are the r biggest singular values, $v_i$ the corresponding right singular vectors and $u_i$
the corresponding left singular vectors of M. $M^*$ denotes the conjugate transpose of $M$. The function compares the action of
the unknown matrix $X$ with the action of M on the singular vectors of $M$ and penalizes the differences in near null-space with
larger weights. It can be formulatet in the form
```math
X M^+(M^+)^* + (M^+)^*M^+X = MM^+(M^+)^* (M^+)^*M^+M,
```
where $M^+$ is the pseudo inverse of $M$.
We can also set `impose_null_spaces` true:

````@example example
sparsify(M, 0.6, 2, 200, true)
````

Then after the described optimization problem was solved the function also solves an additional optimization problem that
ensures that the left and right null-spaces are preserved exactly. Because our matrix has full rank this makes no difference in our case.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

