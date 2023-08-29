# Levitt-Martinsson-HBS

Levitt-Martinsson-HBS implements the HBS matrix recovery algorithm from matrix-vector products described in the preprint _Linear-Complexity Black-Box Randomized Compression of Rank-Structured Matrices_ (J. Levitt and P.G. Martinsson, 2022): https://arxiv.org/abs/2205.02990 

This code relies on the [hm-toolbox](https://github.com/numpi/hm-toolbox) code, which allows us to create and manipulate HODLR matrix objects. We build on this code by adding assignFactor.m, which allows us to access and modify particular off-diagonal blocks of the matrix. 

## Testing

To test the algorithm on a basic HBS matrix, you can run the following few lines of code. First, we generate a semiseparable matrix given by random low-rank factors.

```
N = 64;
u = randn(N, 1);
v = randn(N, 1);
A = diag(randn(N,1)) + triu(u*v',1) + tril(v*u',-1);
```

We define the function handle for matrix-vector multiplication with A.

```
Ax = @(x) A*x;
```

Now, we pass this in along with the matrix size and off-diagonable block rank, 2. 

```
Approx = LevittMartinssonHBS(Ax, N, 2);
```
