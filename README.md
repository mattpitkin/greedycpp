# Overview

A fast, scalable and easy-to-use parallelization of the greedy algorithm for building application-specific basis, empirical interpolants and reduced order quadrature rules. 

In many cases, for example when implemented on a computer and using the Euclidean inner product, the greedy algorithm is QR with column pivoting. And so this code may also be used to compute QR decompositions of a large, dense matrix.

For details on the code's scaling and QR-based model reduction:

[1] Harbir Antil, Dangxing Chen, and Scott Field. "A note on QR-based model reduction". In preparation. 


Code's coordinates: https://sfield83@bitbucket.org/sfield83/greedycpp.git

# Acknowledgements

Priscilla Canizares, Collin Capano, Michael Pürrer, and Rory Smith for careful error reporting, code improvements and testing of early versions of the code. Peter Diener for numerous performance optimizations. See contributors.txt for a complete list of people who have contributed to the code.

# Wiki
View the project's Bitbucket [wiki](https://bitbucket.org/sfield83/greedycpp/wiki/Home) or clone it

    git clone https://bitbucket.org/sfield83/greedycpp.git/wiki