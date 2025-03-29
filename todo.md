# todo notes for `dist_m4ri` program

- [x] Enable `wmax` for `method=1` to speed things up:
  - For quantum codes, only check for degeneracy the codewords of
    weight smaller than `wmax` (suppose we already have a codeword of
    weight `wmax`).
  - Second, when recording indices, stop at `wmax` (as it currently
    done with `minW`).
- [ ] Maintain the count of valid min-weight codewords encountered
      (???) --- will it be compatible with `minW`?
- [ ] Maintain the count of trivial codewords of different weights
      (below `wmax`) encountered (???)
- [ ] Make sure any debugging information is sent to `stderr`, while
      only the distance is printed to `stdout` (positive or negative),
      with perhaps additional data in the same line after a space.
- [ ] Use a faster sorting routine.
- [ ] Use a hash to record constructed codewords and the number of
      times they have been encountered.  Perhaps, bound the hash size
      by a relatively small limit, like it is done in `QDistRnd`, to
      be able to set up an alternative stopping condition.
- [x] To calculate the `confinement` correctly, we need to set up a
      hash of syndrome vectors and corresponding min-weight errors.
      **Do we also need to separate by logical sectors?** ~~Right now
      the degeneracy is not taken into account, i.e., zero syndrome
      weight is calculated for a trivial error of weight 4 (since it
      is equivalent to a zero-weight error this should not be counted).~~
- [ ] make sure `debug` and `method` are read first (as a separate items)
- [ ] given `method`, check parameters for relevancy
- [ ] see if parameters need disambiguation when `method=3` (both) 
- [ ] see what can be done regarding early recursion termination --
      **currently it is removed**.  Perhaps, separate `confinement` to a
      different method, so that full enumeration can be done.
- [ ] see if the constructed hash can be used to find minimum code
      distance better (using half-weight vectors).
## operation mode: distance verification.

We believe we know the distance (in reality, an upper bound), and we
just would like to verify this to be the case.  Say, we believe the
distance is `d=10`, and we are not interested in codes with smaller
distances. Run with `method=1 wmax=10 wmin=9`; the program will not
bother to verify any codewords of weight 10 and above, and will
terminate immediately after finding a valid (non-degenerate in the
case of a quantum code) codeword of weight `9` or below.  For example,
if a codeword of weight `9` has been found, it will print out the
result as `-9`.  If no codewords have been found below `wmax`, the
program would output `0`.

## operation mode: find upper bound on the distance.

Use hash to count the codewords found in the range `wmin < w < wmax`
(do we also need `dw`?)
