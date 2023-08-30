# Levitt-Martinsson-HBS

Levitt-Martinsson-HBS implements the HBS matrix recovery algorithm from matrix-vector products described in the preprint _Linear-Complexity Black-Box Randomized Compression of Rank-Structured Matrices_ (J. Levitt and P.G. Martinsson, 2022): https://arxiv.org/abs/2205.02990 

## Testing

To test the algorithm on a basic HBS matrix, you can run the following few lines of code. First, we generate a semiseparable matrix given by random low-rank factors. The off-diagonal blocks are rank-1 by construction.

```
N = 128;
u = randn(N, 1);
v = randn(N, 1);
A = triu(u*v',0) + tril(v*u',-1);
```

We define the function handle for matrix-vector multiplication with A.

```
Ax = @(x) A*x;
```

Now, we pass this in along with the matrix size and off-diagonable block rank, 1. This gives us the matrix recovered by the Levitt-Martinsson algorithm. 

```
Approx = LevittMartinssonHBS(Ax, N, 1)
```
