This repository contains a collection of scripts used to explore the construction
of sparse discrete partial differential operators in 1D, on triangles, and
tetrahedra, using orthogonal polynomials. In 1D, we consider the Jacobi polynomials. On simplices, we consider the higher-dimensional analogs of Jacobi polynomials, such as those due to Koornwinder on triangles. The goal here is to construct efficient, compressed, and spectrally accurate solvers for a broad class of problems. 

The scripts in this repo can be used to
- compute Jacobi matrices in 1D and on triangles 
- compute multiplication operators in 1D and on triangles
- compute sparse discrete PDOs in 1D and on triangles
- compute Gaussian quadrature for any Jacobi weight in 1D
- compute Gaussian-like quadrature for any Koornwinder weight on the triangle
- solve most second order ODEs
- solve Poisson and variable-coefficient Helmholtz problems on the triangle
- explore the connection between Gaussian-like quadrature and 
  Joint-Eigenvalue-Decomposition algorithms, as employed in signal 
  processing problems like blind source seperation 

My writeup `jacpoly.pdf` gives a more detailed picture.
