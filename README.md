# Regularization solvers for graph matching problem

<img src = "https://github.com/rongxuan-li/graph-match/blob/main/image/cover.png" height="300"/>

This repository contains codes for $L_p$ norm regularization solver and linear reweighted regularization solver for graph matching problem.

Any comments and suggestions are welcome. 

=================================

Solvers:

* `linear_reweighted_solver.m`: Linear reweighted regularization slover for solving graph matching problem.

* `lp_norm_solver.m`: $L_p$ norm regularization slover for solving graph matching problem. [Implemented based on Algorithm 2 in Jiang, B., Liu, Y.F. and Wen, Z., 2016. Lp norm regularization algorithms for optimization over permutation matrices. SIAM Journal on Optimization, 26(4), pp.2284-2313.]

Data:

* `matrix_generator.m`: Main code for generating random matrix and random permutation matrix.

* `dataset_A_B.mat`: One possible set of random matrix A, B and random permutation matrix P.

Results:

* `main_result.m` : Main code for observing the experimental results of different solvers.

<img src = "https://github.com/rongxuan-li/graph-match/blob/main/image/result_plot.png" height="300"/>


