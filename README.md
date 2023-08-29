# Levitt-Martinsson-HBS

Levitt-Martinsson-HBS implements the HBS matrix recovery algorithm from matrix-vector products described in the preprint _Linear-Complexity Black-Box Randomized Compression of Rank-Structured Matrices_ (J. Levitt and P.G. Martinsson, 2022): https://arxiv.org/abs/2205.02990 

This code relies on the [hm-toolbox](https://github.com/numpi/hm-toolbox) code, which allows us to create and manipulate HODLR matrix objects. We build on this code by adding assignFactor.m, which allows us to access and modify particular off-diagonal blocks of the matrix. 
