# Regularization solvers for graph matching problem

<img src = "https://github.com/rongxuan-li/graph-match/blob/main/image/bodymatch_cover.png" height="200" />

This repository contains codes for $L_p$ norm regularization solvers and linear reweighted regularization solvers for graph matching problem.

Any comments and suggestions are welcome. 

==========================================================================================================

Solvers:

* `linear_reweighted_solver.m`: Linear reweighted regularization slover for solving graph matching problem.

* `lp_norm_solver.m`: $L_p$ norm regularization slover for solving graph matching problem. [Implemented based on Algorithm 2 in Jiang, B., Liu, Y.F. and Wen, Z., 2016. Lp-norm regularization algorithms for optimization over permutation matrices. SIAM Journal on Optimization, 26(4), pp.2284-2313.]

Data:

* `matrix_generator.m`: Main code for generating random matrix and random permutation matrix.

* `data_set_A_B.mat`: One possible sets of random generated matrix A, B and random permutation matrix P.

