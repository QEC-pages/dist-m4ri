Author: Leonid P. Pryadko & Weilei Zeng

Date: 2024-08-01

## Overview

The program implements two algorithms for calculating the distance of
a classical linear binary code or a quantum CSS qubit code:

- Random information set (also, random window, or `RW`), to calculate
  the upper distance bound, and
- Connected cluster (`CC`), to calculate the actual distance or lower
  distance bound of an LDPC code (quantum or classical).

For a classical (binary linear) code, only matrix `H`, the
parity-check matrix, should be specified.

For a quantum CSS code, matrix `H=Hx` and either `G=Hz` or `L=Lx`
matrices are needed.

All matrices with entries in `GF(2)` should have the same number of columns,
`n`, and obey the following orthogonality conditions: $$H_XH_Z^T=0,\ \ \ \
H_XL_Z^T=0,\ \ \ \ L_XH_Z^T=0,\ \ \ \ L_XL_Z^T=I,$$ where \f$I\f$ is an identity
matrix.  Notice that the latter identity is not required; it is sufficient that
matrices `Lx`, `Lz`, and the product \f$L_XL_Z^T\f$ have the same full row rank
\f$k\f$, which is also the dimension of the code (number of encoded qubits).


## How it works: RW algorithm (method=1)

Given the error model, i.e., the matrices \f$H=H_x\f$, \f$L=L_x\f$ (\f$L\f$ is empty
for a classical code), the program searches for smallest-weight binary
`codewords` \f$c\f$ such that \f$Hc=0\f$, \f$Lc\neq0\f$.


It repeatedly calculates reduced row echelon form of `H`, with columns
taken in random order, which uniquely fixes the information set
(non-pivot columns).  Generally, column permutation and row reduction
gives \f$H=U\,(I|A)\,P\f$, where \f$U\f$ is invertible, \f$P\f$ is a permutation
matrix, \f$I\f$ is the identity matrix of size given by the rank of \f$H\f$,
and columns of \f$A\f$ form the information set of the corresponding
binary code.  The corresponding codewords can be drawn as the rows of
the matrix \f$(A^T|I')\,P\f$.  Because of the identity matrix \f$I'\f$, the
distribution of such vectors is tilted toward smaller weight, which
qualitatively explains why it works.

To speed up the distance calculation, you can use the parameter `wmin`
(by default, `wmin=1`).  When non-zero, if a code word of weight `w`
\f$\le\f$ `wmin` is found, the distance calculation is terminated
immediately, and the result `-w` with a negative sign is returned.
This is useful, e.g., if we need to construct a code with big enough
distance.

Additional command-line parameters relevant for this method: 

- `steps` the number of RW decoding steps (the number of information
  sets to be constructed).

## How it works: CC algorithm (method=2).

The program tries to construct a codeword recursively, by starting
with a non-zero bit in a position `i` in the range from \f$0\f$ to \f$n-1\f$,
where \f$n\f$ is the number of columns, and then recursively adding the
additional bits in the support of unsatisfied checks starting from the
top.  The complexity to enumerate all codewords of weight up to \f$w\f$
can be estimated as \f$n\,(\Delta-1)^{w-1}\f$, where \f$\Delta\f$ is the
maximum row weight.

Additional command-line parameters relevant for this method: 

- `wmax` the maximum size of the connected cluster.

- `start` the position to start the cluster.  In this case only one
  starting position `i=start` will be used.  This is useful, e.g., if
  the code is symmetric (as, e.g., for cyclic codes).
  
With `debug&2` non-zero, the program in this mode also displays the minimum
weight of the syndrome found for each error weight `w`.  Example (run from
`dist-m4ri/src/`)

```sh
./dist_m4ri debug=3 finH= ../examples/c204H.mmx method=2 wmax=10
# read H <- file '../examples/c204H.mmx'
# recursively searching for w=1 codewords wmax=10 beg=0 end=194
# recursively searching for w=2 codewords wmax=10 beg=0 end=194
# recursively searching for w=3 codewords wmax=10 beg=0 end=194
# recursively searching for w=4 codewords wmax=10 beg=0 end=194
# recursively searching for w=5 codewords wmax=10 beg=0 end=194
# recursively searching for w=6 codewords wmax=10 beg=0 end=194
# recursively searching for w=7 codewords wmax=10 beg=0 end=194
# recursively searching for w=8 codewords wmax=10 beg=0 end=194
# w=1 min syndrome weight 2
# w=2 min syndrome weight 3
# w=3 min syndrome weight 3
# w=4 min syndrome weight 2
# w=5 min syndrome weight 1
# w=6 min syndrome weight 1
# w=7 min syndrome weight 1
# w=8 min syndrome weight 0
### Cluster (actual min-weight codeword found): dmin=8
success  (found min-weight codeword) d=8

```

## How to run it

For help, just run `./dist_m4ri -h` or `./dist_m4ri --help`.  This
shows the following 
```sh
$ ./dist_m4ri --help 
./dist_m4ri: distance of a classical or quantum CSS code
	usage: ./dist_m4ri parameter=value [...]

   Required parameter:
	method=[int]: bitmap for method used (required, default 0: none): 
		1: random window (RW) algorithm. Options:
		   steps=[int]: how many information sets to use (1)
		   wmin=[int]:  minimum distance of interest (1)
		2: connected cluster (CC) algorithm.  Options:
		   wmax=[int]:  maximum cluster weight (5) 
		   start=[int]: use only this position to start (-1)

   General parameters:
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

The program is intended for use with recent `gcc` compilers under linux.
Download the distribution from `github` then run from the `dist-m4ri/src`
directory ```sh make -j all ``` This should compile the executable `dist_m4ri`.

The program uses `m4ri` library for binary linear algebra.  To install
under Ubuntu, run
```
sudo apt-get update -y
sudo apt-get install -y libm4ri-dev
```


## Documentation 
The document for this package is available at
[https://qec-pages.github.io/dist-m4ri/](https://qec-pages.github.io/dist-m4ri/). The
package is hosted on [github](https://github.com/QEC-pages/dist-m4ri)

This file (`introduction.md`) is used as the main page for the doxygen documentation.

