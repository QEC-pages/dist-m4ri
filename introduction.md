Author: Leonid P. Pryadko & Weilei Zeng

Date: 2024-08-01

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

All matrices with entries in `GF(2)` should have the same number of
columns, `n`, and obey the following orthogonality conditions:
$$H_XH_Z^T=0,\quad H_XL_Z^T=0,\quad L_XH_Z^T=0,\quad L_XL_Z^T=I,$$
where \f$I\f$ is an identity matrix.  Notice that the latter identity is
not required; it is sufficient that `Lx` and `Lz` matrices have the
same full row rank `=k`, the dimension of the code, each row of `Lx`
has a non-zero scalar product with a row of `Lz`, and vice versa.



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

## How to run it

For help, just run `./dist_m4ri -h` or `./dist_m4ri --help`.  This
shows the following 
```sh
$ ./dist_m4ri --help 
./dist_m4ri: calculate the minumum distance of a q-LDPC code
        usage: ./dist_m4ri [arguments [...]]
Supported parameters:
        debug=[int]:     bitmap for aux information (3)
        fin=[string]: base name for input files ("try")
                 finH->"${try}X.mtx"  finG->"${try}X.mtx"
        finH=[str]: parity check matrix Hx (NULL)
        finG=[str]: matrix Hz or NULL for classical code (NULL)
        finL=[str]: matrix Lx or NULL for classical code (NULL)
                 Either L=Lx or G=Hz matrix is required for a quantum CSS code
        css=1: this is a CSS code (the only supported one) (1)
        seed=[int]: rng seed  [0 for time(NULL)]
        method=[int]: bitmap for method used:
                1: random window (RW) algorithm
                2: connected cluster (CC) algorithm
        steps=[int]: how many RW decoding cycles to use (1)
        wmax=[int]: max cluster weight in CC (5)
        wmin=[int]: min distance of interest in RW (1)
        -h or --help gives this help
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


## Documentation 
The document for this package is available at [https://qec-pages.github.io/dist-m4ri/](https://qec-pages.github.io/dist-m4ri/). The package is hosted on [github](https://github.com/QEC-pages/dist-m4ri)

This file (`introduction.md`) is used as the main page for the doxygen documentation.

### Formula syntax

Inline formula: here \f$ a \f$ is random variable. \f$J \in {-1,1} \f$, 
and (display mode) \f[ x \in (0,\infty) \f]

Then one can compute \f$y\f$ as (display mode)
$$ y=a \exp(-Jx) $$

$single dollar sign$ is not supported for inline formulas

### Code blocks

A `Python` function can be defined as
```
def add(a,b):
    return a+b
```

### Inline code block `y=\exp(Jx)` is not supported in section titles, how about just use (method = 1)
