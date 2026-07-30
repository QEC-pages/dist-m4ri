---
title: README
author: 
 - Leonid P. Pryadko
 - Weilei Zeng
date: 2024-08-01
---
# dist-m4ri - distance of a classical or quantum CSS code

## Overview

The program implements two algorithms for calculating the distance of
a classical or quantum CSS binary code:

- Random information set (also, random window, or `RW`), to calculate
  the upper distance bound, and
- Connected cluster (`CC`), to calculate the actual distance or lower
  distance bound of an LDPC code (quantum or classical).

For a classical (binary linear) code, only matrix `H`, the
parity-check matrix, should be specified.

For a quantum CSS code, matrix `H=Hx` and either `G=Hz` or `L=Lx`
matrices are needed.

Alternatively, a detector error model (DEM) file from `stim` can be specified
using `fdem=[str]`. This file contains both detector and logical error
associations, which are used to reconstruct `Hx` and `Lx` matrices.

All matrices with entries in `GF(2)` should have the same number of
columns, `n`, and obey the following orthogonality conditions:
$$H_XH_Z^T=0,\quad H_XL_Z^T=0,\quad L_XH_Z^T=0,\quad L_XL_Z^T=I,$$
where $I$ is an identity matrix.  Notice that the latter identity is
not required; it is sufficient that `Lx` and `Lz` matrices have the
same full row rank `=k`, the dimension of the code, each row of `Lx`
has a non-zero scalar product with a row of `Lz`, and vice versa.

## How it works: RW algorithm (`method=1`)

Given the error model, i.e., the matrices $H=H_x$, $L=L_x$ ($L$ is empty
for a classical code), the program searches for smallest-weight binary
`codewords` $c$ such that $Hc=0$, $Lc\neq0$.

It repeatedly calculates reduced row echelon form of `H`, with columns
taken in random order, which uniquely fixes the information set
(non-pivot columns).  Generally, column permutation and row reduction
gives $H=U\,(I|A)\,P$, where $U$ is invertible, $P$ is a permutation
matrix, $I$ is the identity matrix of size given by the rank of $H$,
and columns of $A$ form the information set of the corresponding
binary code.  The corresponding codewords can be drawn as the rows of
the matrix $(A^T|I')\,P$.  Because of the identity matrix $I'$, the
distribution of such vectors is tilted toward smaller weight, which
qualitatively explains why it works.

To speed up the distance calculation, you can use the parameter `wmin` and `wmax`
(by default, `wmin=1` and `wmax=-1`).  When non-zero, if a code word of weight `w`
$\le$ `wmin` is found, the distance calculation is terminated
immediately, and the result `-w` with a negative sign is returned.
This is useful, e.g., if we need to construct a code with big enough
distance.

Additional command-line parameters relevant for this method: 

- `steps` the number of RW decoding steps (the number of information
  sets to be constructed).

## How it works: CC algorithm (`method=2`).

The program tries to construct a codeword recursively, by starting
with a non-zero bit in a position `i` in the range from $0$ to $n-1$,
where $n$ is the number of columns, and then recursively adding the
additional bits in the support of unsatisfied checks starting from the
top.  The complexity to enumerate all codewords of weight up to $w$
can be estimated as $n\,(\Delta-1)^{w-1}$, where $\Delta$ is the
maximum row weight.

Additional command-line parameters relevant for this method: 

- `wmax` the maximum size of the connected cluster.

- `start` the position to start the cluster.  In this case only one
  starting position `i=start` will be used.  This is useful, e.g., if
  the code is symmetric (as, e.g., for cyclic codes).

- `noscan` if set to 1, start the recursion directly with `w=wmax`,
  without scanning smaller weights.  Only works with `method=2`.
  
### Calculating confinement profile 

With `method=2`, the program also calculates the confinement profile, minimum
syndrome weight for a given error weight.  Notice that it should properly be
calculated for a given `irreducible` error weight, thus the confinement profile
may only be trusted up to the half of the distance.

```sh
$ ./dist_m4ri method=2 finH= ../examples/surf_d5_H.mmx finL= ../examples/surf_d5_L.mmx wmax=5 debug=0
# confinement: 1,1,1,1,1
5

```

```sh
$ ./dist_m4ri method=2 finH= ../examples/QX150.mtx finG= ../examples/QZ150.mtx wmax=8
# read H <- file '../examples/QX150.mtx'
# read G <- file '../examples/QZ150.mtx'
# recursively searching for w=1 codewords wmax=8 beg=0 end=149
# recursively searching for w=2 codewords wmax=8 beg=0 end=148
# recursively searching for w=3 codewords wmax=8 beg=0 end=147
# recursively searching for w=4 codewords wmax=8 beg=0 end=146
# recursively searching for w=5 codewords wmax=8 beg=0 end=145
# recursively searching for w=6 codewords wmax=8 beg=0 end=144
# w=1 min non-zero syndrome weight 2
# w=2 min non-zero syndrome weight 2
# w=3 min non-zero syndrome weight 2
# w=4 min non-zero syndrome weight 2
# w=5 min non-zero syndrome weight 2
# w=6 min non-zero syndrome weight 2
### Cluster (actual min-weight codeword found): d=6

```


## How to run it

For help, just run `./dist_m4ri -h` or `./dist_m4ri --help`.  This
shows the following 
```sh
$ ./dist_m4ri --help 
src/dist_m4ri: distance of a classical or quantum CSS code
	usage: src/dist_m4ri parameter=value [...]

   Required parameter:
	method=[int]: bitmap for method used (no default): 

		1: random window (RW) algorithm. Options:
		   steps=[int]: how many information sets to use (1)
		   wmin=[int]:  minimum distance of interest (1)
			 immediately stop and return '-w' on a cw of weight w<=wmin
			 use this option to quickly scan over a large number of codes
		   dmax=[int]:  if non-zero, ignore vectors of this and larger wgt (0)
			 this option accelerates the search somewhat

		2: connected cluster (CC) algorithm.  Options:
		   wmax=[int]:  maximum cluster weight to construct, inclusive (0)
			 must be non-zero for CC only, otherwise use upper bound from RW
		   smax=[int]:  maximum syndrome weight of interest, inclusive (20)
			 must be non-zero to calculate confinement profile
		   start=[int]: use only this position to start (-1)
		   noscan=[int]: start CC directly with wmax (0)

   General parameters:
	fdem=[str]: detector error model (DEM) file from stim (NULL)
	pmin=[float]: minimum error probability to keep for DEM (0.0)
	finH=[str]: parity check matrix Hx (NULL)
	finG=[str]: matrix Hz (quantum CSS code only) (NULL)
	finL=[str]: matrix Lx (quantum CSS code only) (NULL)
		 Either L=Lx or G=Hz matrix is required for a quantum CSS code
	fin=[str]:  base name for input files ("try")
		 set finH->"${fin}X.mtx"  finG->"${fin}Z.mtx"
	css=[int]:  reserved for future use (1)
	seed=[int]: rng seed [use 0 for time(NULL)] (0)
	debug=[int]:	 bitmap for aux information to output (3)
		0: clear the entire debug bitmap to 0.
		1: output misc general info (on by default)
		2: output more general info (on by default)
		4: debug command line arguments parsing
		8: output progress reports every 1000 steps
		16: output new min-weight codewords found (cut large vectors)
		32: output matrices (unless n is large)
		64: reserved
		128: reserved
		256: print out neighbor lists
		512: print out vectors/syndrome weights during recursion
		1024: print piv/skip_pivs/reserved
		2048: allow big matrix / large vector output
		   see the source code for more options
	  Multiple 'debug' parameters are XOR combined except for 0.
	  Use debug=0 as the 1st argument to suppress all debug messages.
   -h gives this help (also '--help')
```

## Compilation

The program is intended for use with recent `gcc` compilers under
linux.  Download the distribution from `github` then run from the
`dist-m4ri/src` directory ```sh make -j all ``` This should compile
the executable `dist_m4ri`.

The program uses `m4ri` library for binary linear algebra.  To install
under Ubuntu, run
```
sudo apt-get update -y
sudo apt-get install -y libm4ri-dev
```

## References

If you use this program, please cite:

*   A. Dumer, A. A. Kovalev, and L. P. Pryadko, "Distance verification for classical and quantum LDPC codes," *IEEE Transactions on Information Theory*, vol. 63, no. 7, pp. 4675-4690, 2017. [doi:10.1109/TIT.2017.2690381](https://doi.org/10.1109/TIT.2017.2690381).

Other related papers and software:

*   **vecdec Repository** (implements the Random Information Set (RW) algorithm and can handle different error weights/probabilities):
    [QEC-pages/vecdec](https://github.com/QEC-pages/vecdec).

*   **QDistRnd GAP Package** (describing the Random Information Set algorithm for quantum codes over arbitrary finite fields):
    L. P. Pryadko, V. A. Shabashov, and V. K. Kozin, "QDistRnd: A GAP package for computing the distance of quantum error-correcting codes," *Journal of Open Source Software*, vol. 7, no. 71, p. 4120, 2022. [doi:10.21105/joss.04120](https://doi.org/10.21105/joss.04120).

*   **Performance Comparison** (comparing the performance of this program with other available distance-finding programs):
    M. Webster, A. Jacob, and O. Higgott, "Distance-Finding Algorithms for Quantum Codes and Circuits," arXiv:2603.22532 [quant-ph], 2026. [arXiv:2603.22532](https://arxiv.org/abs/2603.22532).
