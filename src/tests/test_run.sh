#!/bin/bash

# Get directories
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
SRC_DIR="$SCRIPT_DIR/.."
WS_ROOT="$SCRIPT_DIR/../.."
BIN="$SRC_DIR/dist_m4ri_old"
BIN_FORK="$SRC_DIR/dist_m4ri"
EXAMPLES_DIR="$WS_ROOT/examples"

if [ ! -f "$BIN" ]; then
    echo "Error: dist_m4ri_old binary not found at $BIN"
    exit 1
fi

if [ ! -f "$BIN_FORK" ]; then
    echo "Error: dist_m4ri binary not found at $BIN_FORK"
    exit 1
fi

FAILED=0

assert_output() {
    local cmd="$1"
    local expected_exit="$2"
    local expected_out_regex="$3"
    local expected_err_regex="$4"
    
    echo "Running: $cmd"
    local stdout_file=$(mktemp)
    local stderr_file=$(mktemp)
    
    eval "$cmd" > "$stdout_file" 2> "$stderr_file"
    local exit_code=$?
    
    if [ $exit_code -ne $expected_exit ]; then
        echo "  [FAIL] Expected exit code $expected_exit, got $exit_code"
        FAILED=1
    else
        if [ -n "$expected_out_regex" ]; then
            if ! grep -q -E "$expected_out_regex" "$stdout_file"; then
                echo "  [FAIL] Stdout did not match regex: $expected_out_regex"
                echo "         Got: $(cat $stdout_file)"
                FAILED=1
            fi
        fi
        if [ -n "$expected_err_regex" ]; then
            if ! grep -q -E "$expected_err_regex" "$stderr_file"; then
                echo "  [FAIL] Stderr did not match regex: $expected_err_regex"
                echo "         Got: $(cat $stderr_file)"
                FAILED=1
            fi
        fi
    fi
    
    rm -f "$stdout_file" "$stderr_file"
}

# Test 1: CC baseline
assert_output "$BIN method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx wmax=5 debug=0" \
    0 "^5$" ""

# Test 2: CC noscan
assert_output \
    "$BIN method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx wmax=5 noscan=1 debug=0" \
    0 "^5$" ""

# Test 3: CC fdem baseline
assert_output "$BIN method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=3 debug=0" 0 "^3$" ""

# Test 4: CC fdem with pmin
assert_output "$BIN method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=3 pmin=0.01 debug=0" 0 "^3$" ""

# Test 5: CC fdem with high pmin
assert_output "$BIN method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=3 pmin=0.02 debug=0" 0 "^-3$" ""

# Test 6: RW fdem baseline
assert_output "$BIN method=1 fdem=$EXAMPLES_DIR/surf_d3.dem steps=100 debug=0" 0 "^3$" ""

# Test 7: Validation: noscan=1 with method=1
assert_output \
    "$BIN method=1 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx wmax=5 noscan=1 debug=0" \
    255 "" "noscan=1 only works with method=2"

# Test 8: Validation: noscan=1 with method=3
assert_output \
    "$BIN method=3 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx wmax=5 noscan=1 debug=0" \
    255 "" "noscan=1 only works with method=2"

# Test 9: Validation: fdem and finH together
assert_output "$BIN method=2 fdem=$EXAMPLES_DIR/surf_d3.dem finH=$EXAMPLES_DIR/surf_d5_H.mmx wmax=3 debug=0" \
    255 "" "Cannot specify matrix files.*along with fdem"

# Test 10: Validation: pmin without fdem
assert_output "$BIN method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx wmax=3 pmin=0.01 debug=0" \
    255 "" "pmin can only be used when fdem is specified"

# Create temp DEM file with repeat blocks
TEMP_DEM=$(mktemp --suffix=.dem)
cat << 'EOF' > "$TEMP_DEM"
# Detector Error Model
error(0.1) D0 D1
repeat 2 {
  error(0.2) D0 D1
  shift_detectors 2
}
error(0.3) D0 ^ L0
error(0.4) D0
EOF

# Test 11: Repeat and shift_detectors parser verification
assert_output "$BIN method=2 fdem=$TEMP_DEM wmax=2 debug=0" 0 "^2$" ""

rm -f "$TEMP_DEM"

# Test 12: Regression test for nextelement out-of-bounds read / segfault
SEGFAULT_H=$(mktemp --suffix=.mmx)
python3 "$SCRIPT_DIR/gen_segfault_matrix.py" > "$SEGFAULT_H"
assert_output "$BIN method=1 finH=$SEGFAULT_H steps=1 debug=0" 0 "^65$" ""
rm -f "$SEGFAULT_H"

# Test 13: CC codeword saving and loading
TEMP_CWS1=$(mktemp --suffix=.nz)
TEMP_CWS2=$(mktemp --suffix=.nz)

# Step 1: Run and save
assert_output "$BIN method=2 fdem=$EXAMPLES_DIR/surf_d3.dem drop=0 wmax=3 outC=$TEMP_CWS1 debug=0" 0 "^3$" ""

# Step 2: Verify file 1 is not empty
if [ ! -s "$TEMP_CWS1" ]; then
    echo "  [FAIL] Codewords file is empty"
    FAILED=1
fi

# Step 3: Run again loading file 1 and saving to file 2, check for read message
assert_output "$BIN debug=33 method=2 fdem=$EXAMPLES_DIR/surf_d3.dem drop=0 wmax=3 finC=$TEMP_CWS1 outC=$TEMP_CWS2" \
    0 "" "read 128 codewords from"

# Step 4: Verify files are identical
if ! diff -q "$TEMP_CWS1" "$TEMP_CWS2" >/dev/null; then
    echo "  [FAIL] Codewords files differ"
    FAILED=1
fi

rm -f "$TEMP_CWS1" "$TEMP_CWS2"

# Test 14: Codeword verification (orthogonality checks)
INVALID_CWS=$(mktemp --suffix=.nz)
cat << 'EOF' > "$INVALID_CWS"
%% NZLIST
% valid
3 1 3 23
% invalid (not orthogonal to H)
1 1
% trivial (orthogonal to H and L)
5 1 2 7 17 28
EOF

echo "Running Test 14: Codeword verification"
STDOUT_FILE=$(mktemp)
STDERR_FILE=$(mktemp)
$BIN method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=3 finC=$INVALID_CWS debug=0 > "$STDOUT_FILE" 2> "$STDERR_FILE"
EXIT_CODE=$?

if [ $EXIT_CODE -ne 0 ]; then
    echo "  [FAIL] Expected exit code 0, got $EXIT_CODE"
    FAILED=1
fi
if ! grep -q "skipped 2 invalid codewords" "$STDERR_FILE"; then
    echo "  [FAIL] Warning message missing"
    FAILED=1
fi
if ! grep -q -E "^3$" "$STDOUT_FILE"; then
    echo "  [FAIL] Expected distance 3 output missing"
    FAILED=1
fi

rm -f "$INVALID_CWS" "$STDOUT_FILE" "$STDERR_FILE"

# Test 15: Classical mode auto-detection (H only)
assert_output "$BIN method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx wmax=3 debug=0" 0 "^2$" ""

# Test 16: Conflict detection (classical=0 with only H)
assert_output "$BIN method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx wmax=3 classical=0 debug=0" \
    255 "" "L matrix.*is required for quantum code"

# Test 17: Conflict detection (classical=1 with finL)
assert_output \
    "$BIN method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx wmax=3 classical=1 debug=0" \
    255 "" "Conflict: classical=1 specified"

# Test 18: Discarding L with fdem and classical=1
assert_output "$BIN debug=1 method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=3 classical=1" 0 "" "discarding L matrix"

# Test 19: Coordinate general integer matrix format
assert_output "$BIN method=2 finH=$SCRIPT_DIR/test_crd_gen_int.mmx wmax=2 debug=0" 0 "^2$" ""
assert_output "$BIN method=2 finH=$SCRIPT_DIR/test_crd_gen_int.mmx wmax=1 debug=0" 0 "^-1$" ""

# Test 20: Coordinate symmetric integer matrix format (doubling off-diagonal)
assert_output "$BIN method=2 finH=$SCRIPT_DIR/test_crd_sym_int.mmx wmax=2 debug=0" 0 "^2$" ""
assert_output "$BIN method=2 finH=$SCRIPT_DIR/test_crd_sym_int.mmx wmax=1 debug=0" 0 "^-1$" ""

# Test 21: Coordinate general pattern matrix format (no values)
assert_output "$BIN method=2 finH=$SCRIPT_DIR/test_crd_gen_pat.mmx wmax=2 debug=0" 0 "^2$" ""
assert_output "$BIN method=2 finH=$SCRIPT_DIR/test_crd_gen_pat.mmx wmax=1 debug=0" 0 "^-1$" ""

# Test 22: Array general integer matrix format (dense general)
assert_output "$BIN method=2 finH=$SCRIPT_DIR/test_arr_gen_int.mmx wmax=2 debug=0" 0 "^2$" ""
assert_output "$BIN method=2 finH=$SCRIPT_DIR/test_arr_gen_int.mmx wmax=1 debug=0" 0 "^-1$" ""

# Test 23: Array symmetric integer matrix format (dense symmetric, doubling off-diagonal)
assert_output "$BIN method=2 finH=$SCRIPT_DIR/test_arr_sym_int.mmx wmax=2 debug=0" 0 "^2$" ""
assert_output "$BIN method=2 finH=$SCRIPT_DIR/test_arr_sym_int.mmx wmax=1 debug=0" 0 "^-1$" ""

# =========================================================================
# dist_m4ri multithreading tests (methods 1, 2, 3, dexp, timeout, outC)
# =========================================================================

# Test 24: dist_m4ri method=2 (multithreaded CC exact distance, rw_steps=0)
assert_output \
    "$BIN_FORK method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx wmax=5 debug=0 threads=4" \
    0 "^5 5 0$" ""

# Test 25: dist_m4ri method=2 (CC lower bound when distance > wmax, rw_steps=0)
assert_output \
    "$BIN_FORK method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx wmax=3 debug=0 threads=4" \
    0 "^4 0 0$" ""

# Test 26: dist_m4ri method=1 (multithreaded RW, reported rw_steps > 0)
assert_output "$BIN_FORK method=1 fdem=$EXAMPLES_DIR/surf_d3.dem steps=100 debug=0 threads=4" 0 "^1 3 [0-9]+$" ""

# Test 27: dist_m4ri method=3 (bracketing mode with dexp)
assert_output "$BIN_FORK method=3 fdem=$EXAMPLES_DIR/surf_d3.dem dexp=3 timeout=10 debug=0 threads=4" \
    0 "^3 3 [0-9]+$" ""

# Test 28: dist_m4ri method=3 (bracketing mode with dest alias)
assert_output "$BIN_FORK method=3 fdem=$EXAMPLES_DIR/surf_d3.dem dest=3 timeout=10 debug=0 threads=4" \
    0 "^3 3 [0-9]+$" ""

# Test 29: dist_m4ri method=3 on surf_d5
assert_output \
    "$BIN_FORK method=3 finH=$EXAMPLES_DIR/surf_d5_H.mmx \
finL=$EXAMPLES_DIR/surf_d5_L.mmx dexp=5 timeout=10 debug=0 threads=4 steps=500" \
    0 "^5 5 [0-9]+$" ""

# Test 30: dist_m4ri codeword saving and loading
TEMP_FORK_CWS1=$(mktemp --suffix=.nz)
TEMP_FORK_CWS2=$(mktemp --suffix=.nz)
assert_output "$BIN_FORK method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=3 outC=$TEMP_FORK_CWS1 debug=0 threads=4" \
    0 "^3 3 0$" ""
assert_output \
    "$BIN_FORK method=3 fdem=$EXAMPLES_DIR/surf_d3.dem \
wmax=3 finC=$TEMP_FORK_CWS1 outC=$TEMP_FORK_CWS2 debug=0 threads=4" \
    0 "^3 3 [0-9]+$" ""
if ! diff -q "$TEMP_FORK_CWS1" "$TEMP_FORK_CWS2" >/dev/null; then
    echo "  [FAIL] dist_m4ri codeword files differ"
    FAILED=1
fi
rm -f "$TEMP_FORK_CWS1" "$TEMP_FORK_CWS2"

# Test 31: dist_m4ri timeout graceful termination
assert_output \
    "$BIN_FORK method=1 finH=$EXAMPLES_DIR/surf_d5_H.mmx \
finL=$EXAMPLES_DIR/surf_d5_L.mmx steps=10000000 timeout=0.5 debug=0 threads=4" \
    0 "^1 5 [0-9]+$" ""

# Test 32: dist_m4ri classical mode
assert_output "$BIN_FORK method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx wmax=3 debug=0 threads=4" 0 "^2 2 0$" ""

# Test 33: dist_m4ri method=2 with dW=1 and outC
TEMP_DW_CWS=$(mktemp --suffix=.nz)
assert_output "$BIN_FORK debug=15 method=2 fdem=$EXAMPLES_DIR/surf_d3.dem dW=1 wmax=4 outC=$TEMP_DW_CWS threads=4" \
    0 "^3 3 0$" "continuing up to w=4 for dW=1"
rm -f "$TEMP_DW_CWS"

# Test 34: dist_m4ri_old method=2 with dW=1 reporting
TEMP_M4RI_CWS=$(mktemp --suffix=.nz)
assert_output "$BIN debug=1 method=2 fdem=$EXAMPLES_DIR/surf_d3.dem dW=1 wmax=4 outC=$TEMP_M4RI_CWS" \
    0 "^3$" "CC round w=.*searched with dW=1"
rm -f "$TEMP_M4RI_CWS"

# Test 35: dist_m4ri method=2 timeout lower bound correctness
assert_output \
    "$BIN_FORK method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx \
finL=$EXAMPLES_DIR/surf_d5_L.mmx wmax=5 timeout=0.0001 debug=0 threads=4" \
    0 "^[1-5] 0 0$" ""

# Test 36: dist_m4ri method=3 timeout bound correctness
assert_output \
    "$BIN_FORK method=3 finH=$EXAMPLES_DIR/surf_d5_H.mmx \
finL=$EXAMPLES_DIR/surf_d5_L.mmx wmax=5 timeout=0.0001 debug=0 threads=4" \
    0 "^[1-5] [0-5] [0-9]+$" ""

# Test 37: dmin parameter starting CC search directly at dmin
assert_output \
    "$BIN_FORK method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx \
finL=$EXAMPLES_DIR/surf_d5_L.mmx dmin=5 wmax=5 debug=0 threads=4" \
    0 "^5 5 0$" ""

# Test 38: dmax parameter in method 1
assert_output \
    "$BIN_FORK method=1 finH=$EXAMPLES_DIR/surf_d5_H.mmx \
finL=$EXAMPLES_DIR/surf_d5_L.mmx dmax=5 steps=100 debug=0 threads=4" \
    0 "^1 5 [0-9]+$" ""

# Test 39: dmin and dmax in method 3
assert_output \
    "$BIN_FORK method=3 finH=$EXAMPLES_DIR/surf_d5_H.mmx \
finL=$EXAMPLES_DIR/surf_d5_L.mmx dmin=4 dmax=5 timeout=5 debug=0 threads=4" \
    0 "^5 5 [0-9]+$" ""

# Test 40: multiple debug arguments are OR-combined (decimal or hexadecimal), the first one replaces the default;
# debug=0 alone is silent; invalid values are rejected
assert_output "$BIN_FORK debug=1 debug=0x40 method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=2 threads=2" \
    0 "^3 0 0$" "^# debug=65 \(0x41\)$"
assert_output "$BIN debug=1 debug=2 method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=2" 0 "^-2$" ""
assert_output "$BIN_FORK debug=-1 method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=2" 255 "" "invalid 'debug=-1'"
echo "Running Test 40: debug=0 is silent"
ERR40=$($BIN_FORK debug=0 method=3 fdem=$EXAMPLES_DIR/surf_d3.dem threads=2 2>&1 >/dev/null)
if [ -n "$ERR40" ]; then
    echo "  [FAIL] Test 40: debug=0 gives stderr '$ERR40'"
    FAILED=1
fi

# Test 41: dist_m4ri error on no arguments (short help)
assert_output "$BIN_FORK" 255 "" "Allowed parameters:"

# Test 42: dist_m4ri error on unrecognized parameter (short help)
assert_output "$BIN_FORK unrecognized_param=123" 255 "" \
    "unrecognized parameter \"unrecognized_param=123\" at position 1"

# Test 43: dist_m4ri --help (returns 0, points to --morehelp)
assert_output "$BIN_FORK --help" 0 "morehelp" ""

# Test 44: dist_m4ri --morehelp (returns 0, lists all parameters including classical, finC, outC)
assert_output "$BIN_FORK --morehelp" 0 "classical=\[0\|1\]" ""

# Test 45: dist_m4ri_old --help
assert_output "$BIN --help" 0 "morehelp" ""

# Test 46: dist_m4ri_old --morehelp
assert_output "$BIN --morehelp" 0 "classical=\[0\|1\]" ""

# Test 47: dist_m4ri --version
assert_output "$BIN_FORK --version" 0 "dist_m4ri version 0.15.0" ""

# Test 48: dist_m4ri_old --version
assert_output "$BIN --version" 0 "dist_m4ri version 0.15.0" ""

# Test 49: dist_m4ri RW with ksub subspace sketching
assert_output "$BIN_FORK method=1 fdem=$EXAMPLES_DIR/surf_d3.dem steps=200 ksub=32 debug=0 threads=4" \
    0 "^1 3 [0-9]+$" ""

# Test 50: dist_m4ri RW with localized window permutation (kwin / win_mode)
assert_output \
    "$BIN_FORK method=1 fdem=$EXAMPLES_DIR/surf_d3.dem steps=200 ksub=32 kwin=48 win_mode=0 debug=0 threads=4" \
    0 "^1 3 [0-9]+$" ""

# Test 51: dist_m4ri RW with min_hits and cov_cws early convergence
assert_output \
    "$BIN_FORK debug=1 method=1 fdem=$EXAMPLES_DIR/surf_d3.dem steps=5000 ksub=32 min_hits=3 cov_cws=2 threads=4" \
    0 "^1 3 [0-9]+$" "RW convergence reached"

# Test 52: dist_m4ri RW with periodic adaptive basis refresh
assert_output \
    "$BIN_FORK method=1 finH=$EXAMPLES_DIR/surf_d5_H.mmx \
finL=$EXAMPLES_DIR/surf_d5_L.mmx steps=300 ksub=32 refresh=50 debug=0 threads=4" \
    0 "^1 5 [0-9]+$" ""

# Test 53: dist_m4ri_old RW with ksub, kwin, min_hits, and refresh
assert_output "$BIN method=1 fdem=$EXAMPLES_DIR/surf_d3.dem steps=200 ksub=32 win=48 min_hits=3 refresh=50 debug=0" \
    0 "^3$" ""

# Test 54: dist_m4ri ksub automatic fallback warning when m < nu (2m < n)
assert_output "$BIN_FORK debug=1 method=1 fdem=$EXAMPLES_DIR/surf_d3.dem drop=0 steps=200 ksub=32 threads=4" \
    0 "^1 3 [0-9]+$" "Warning: ksub=32 requested, but m=24 < n-m=197 \(<= nu\); falling back to full-matrix RW"

# Test 55: dist_m4ri method=3 with min_hits stops only RW workers and lets CC finish exact proof
assert_output "$BIN_FORK debug=3 method=3 fdem=$EXAMPLES_DIR/surf_d3.dem steps=50000 min_hits=2 cov_cws=2 threads=4" \
    0 "^3 3 [0-9]+$" "RW convergence reached"

# Test 56: default method=3 and smax=0 warning in .mtx mode under debug&1==1
assert_output \
    "$BIN_FORK debug=1 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx wmax=5 threads=4 steps=200" \
    0 "^5 5 [0-9]+$" "Warning: smax=0, confinement profile is not computed"

# Test 57: no smax=0 warning in DEM mode under debug&1==1
echo "Running Test 57: no smax=0 warning in DEM mode"
STDERR_DEM=$(mktemp)
$BIN_FORK debug=1 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=3 threads=4 > /dev/null 2> "$STDERR_DEM"
if grep -q "confinement profile is not computed" "$STDERR_DEM"; then
    echo "  [FAIL] Unexpected smax=0 warning in DEM mode"
    FAILED=1
fi
rm -f "$STDERR_DEM"

# Expert CC options: noscan, start (list, unlimited clusters), cbeg/cend (split runs)
S5="finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx"

# Test 58: noscan=1 without dmin: codeword at w=wmax only gives dmax (dmin not certified), warning
assert_output "$BIN_FORK method=2 $S5 wmax=5 noscan=1 debug=0 threads=4" \
    0 "^1 5 0$" "WARNING: noscan=1 \(expert option\)"

# Test 59: noscan=1 with dmin=wmax: certified, exact
assert_output "$BIN_FORK method=2 $S5 wmax=5 dmin=5 noscan=1 debug=0 threads=4" \
    0 "^5 5 0$" "the supplied dmin=5 being a certified lower bound"

# Test 60: noscan=1 without codewords at w=wmax: dmin not raised
assert_output "$BIN_FORK method=2 $S5 wmax=4 noscan=1 debug=0 threads=4" 0 "^1 0 0$" "WARNING: noscan=1"

# Test 61: dist_m4ri_old noscan=1 without codewords at w=wmax prints 0 (no lower bound certified)
assert_output "$BIN method=2 $S5 wmax=4 noscan=1 debug=0" 0 "^0$" "WARNING: noscan=1"

# Test 62: start list (unlimited clusters) with a warning
assert_output "$BIN_FORK method=2 $S5 wmax=7 start=0,41 debug=0 threads=4" \
    0 "^5 5 0$" "WARNING: start=0,41 \(expert option\)"

# Test 63: start list with all columns reproduces the full CC result
S3="finH=$EXAMPLES_DIR/surf_d3_H.mmx finL=$EXAMPLES_DIR/surf_d3_L.mmx"
START_ALL=$(seq -s, 0 290)
assert_output "$BIN_FORK method=2 $S3 wmax=4 start=$START_ALL debug=0 threads=4" \
    0 "^3 3 0$" "WARNING: start=0,1,2,3,4,5,6,7,\.\.\."

# Test 64: dist_m4ri_old start list (clusters from column 12 are not limited to columns > 12)
assert_output "$BIN method=2 $S5 wmax=7 start=12 debug=0" 0 "^6$" "WARNING: start=12"

# Test 65: cbeg/cend split run with a warning
assert_output "$BIN_FORK method=2 $S5 wmax=7 cbeg=30 debug=0 threads=4" \
    0 "^5 5 0$" "WARNING: cbeg/cend \(expert option for split runs\)"

# Test 66: invalid start list
assert_output "$BIN_FORK method=2 $S5 wmax=7 start=1,,2 debug=0" 255 "" "invalid start='1,,2'"

# Test 67: start together with cbeg
assert_output "$BIN_FORK method=2 $S5 wmax=7 start=0 cbeg=3 debug=0" 255 "" \
    "Cannot specify start along with cbeg or cend"

# Test 68: start column out of range
assert_output "$BIN_FORK method=2 $S3 wmax=3 start=0,291 debug=0" 255 "" \
    "start column 291 cannot be larger than nvar-1=290"

# Test 69: legacy start=-1 (all columns) is still accepted, without a warning
echo "Running Test 69: legacy start=-1 means all columns"
OUT69=$(mktemp)
ERR69=$(mktemp)
$BIN_FORK method=2 $S5 wmax=7 start=-1 debug=0 threads=4 > "$OUT69" 2> "$ERR69"
if ! grep -q -E "^5 5 0$" "$OUT69" || grep -q "WARNING: start" "$ERR69"; then
    echo "  [FAIL] start=-1: stdout '$(cat "$OUT69")', stderr '$(cat "$ERR69")'"
    FAILED=1
fi
rm -f "$OUT69" "$ERR69"

# Test 70: method=2 collecting codewords (outC or maxC, no wmax): no CC rounds after the distance is found
echo "Running Test 70: method=2 with outC / maxC stops after the round w=d"
CWS70=$(mktemp --suffix=.nz)
for EXTRA in "outC=$CWS70" "maxC=1000"; do
    OUT70=$(mktemp)
    ERR70=$(mktemp)
    $BIN_FORK method=2 $S3 $EXTRA timeout=10 debug=3 threads=2 > "$OUT70" 2> "$ERR70"
    if ! grep -q -E "^3 3 0$" "$OUT70" || grep -q "extra dW round" "$ERR70" || \
       ! grep -q "CC found min-weight codeword: d=3" "$ERR70"; then
        echo "  [FAIL] $EXTRA: stdout '$(cat "$OUT70")', stderr '$(cat "$ERR70")'"
        FAILED=1
    fi
    rm -f "$OUT70" "$ERR70"
done
if [ ! -s "$CWS70" ]; then
    echo "  [FAIL] Test 70: codewords file is empty"
    FAILED=1
fi

# Test 71: method=2 with a known upper bound dmax: no CC round w=dmax once dmin=dmax
echo "Running Test 71: method=2 skips the round w=dmax once dmin=dmax"
OUT71=$(mktemp)
ERR71=$(mktemp)
$BIN_FORK method=2 $S5 dmax=5 wmax=5 debug=11 threads=4 > "$OUT71" 2> "$ERR71"
if ! grep -q -E "^5 5 0$" "$OUT71" || grep -q "searching w=5" "$ERR71" || \
   ! grep -q "bounds coincide: dmin = dmax = 5" "$ERR71"; then
    echo "  [FAIL] dmax=5: stdout '$(cat "$OUT71")', stderr '$(cat "$ERR71")'"
    FAILED=1
fi
rm -f "$OUT71" "$ERR71"

# Test 72: method=2 with the upper bound dmax from finC codewords (no wmax): same as Test 71
CWS72=$(mktemp --suffix=.nz)
$BIN_FORK method=2 $S5 wmax=5 outC=$CWS72 debug=0 threads=4 > /dev/null 2>&1
assert_output "$BIN_FORK method=2 $S5 finC=$CWS72 timeout=10 debug=1 threads=4" \
    0 "^5 5 0$" "bounds coincide: dmin = dmax = 5"
rm -f "$CWS70" "$CWS72"

# Test 73: method=3 with few RW steps limits only the RW threads, CC rounds use all threads (debug=11: timing)
assert_output "$BIN_FORK method=3 $S5 steps=10 threads=8 debug=11" \
    0 "^5 5 [0-9]+$" "RW limited to 1 of 8 threads"
assert_output "$BIN_FORK method=3 $S5 steps=10 threads=8 debug=11" \
    0 "^5 5 [0-9]+$" "CC round w=4 started: 8 CC threads"

# Test 74: method=3 with steps=0 runs pure CC on all threads without RW threads
assert_output "$BIN_FORK method=3 $S5 steps=0 threads=4 debug=11" 0 "^5 5 0$" "steps=0, no RW"

# Test 75: method=3 with a low dexp: CC rounds w > dexp are paused (not ended) until RW finds a codeword
assert_output "$BIN_FORK method=3 finH=$EXAMPLES_DIR/c1920H.mmx dexp=2 timeout=2 threads=4 debug=0" \
    0 "^([4-9]|[1-9][0-9]) [0-9]+ [0-9]+$" ""

# Test 76: method=3: the round w=dmax-1, which certifies dmin=dmax, runs on all threads
assert_output "$BIN_FORK method=3 $S5 dmax=5 threads=8 nothrottle=1 debug=11" \
    0 "^5 5 [0-9]+$" "CC round w=4 started: 8 CC threads"

# Test 77: method=3 with a supplied dmin: the first CC round (work not known yet) runs on half of the threads
assert_output "$BIN_FORK method=3 $S5 dmin=3 min_hits=0 threads=8 nothrottle=1 debug=11" \
    0 "^5 5 [0-9]+$" "CC round w=3 started: 4 CC threads"

# Test 78: min_hits with a single codeword (classical [5,1] repetition code, k=1): RW converges after a few steps
REP78=$(mktemp --suffix=.mtx)
cat << 'EOF' > "$REP78"
%%MatrixMarket matrix coordinate integer general
4 5 8
1 1 1
1 2 1
2 2 1
2 3 1
3 3 1
3 4 1
4 4 1
4 5 1
EOF
echo "Running Test 78: min_hits convergence with a single codeword (k=1)"
OUT78=$(mktemp)
ERR78=$(mktemp)
$BIN_FORK method=1 finH=$REP78 steps=100000 threads=2 debug=1 > "$OUT78" 2> "$ERR78"
if ! grep -q -E "^1 5 ([0-9]|[1-9][0-9])$" "$OUT78" || ! grep -q "RW convergence reached" "$ERR78"; then
    echo "  [FAIL] stdout '$(cat "$OUT78")', stderr '$(cat "$ERR78")'"
    FAILED=1
fi
rm -f "$REP78" "$OUT78" "$ERR78"

# Test 79: method=3 with outC (dW=0) exports all minimum-weight codewords, as method=2 (CC round w=d after dmin=dmax)
CWS79=$(mktemp --suffix=.nz)
assert_output "$BIN_FORK method=3 fdem=$EXAMPLES_DIR/surf_d3.dem drop=0 outC=$CWS79 debug=1 threads=4" \
    0 "^3 3 [0-9]+$" "exported 128 codewords"
if [ "$(grep -c '^3 ' "$CWS79")" != "128" ]; then
    echo "  [FAIL] Test 79: expected 128 codewords of weight 3 in the exported file"
    FAILED=1
fi
rm -f "$CWS79"

# Test 80: codewords of weight 1 (H with zero columns): with outC, all of them are collected (method=2 and method=3)
ZC80=$(mktemp --suffix=.mtx)
cat << 'EOF' > "$ZC80"
%%MatrixMarket matrix coordinate integer general
2 5 4
1 1 1
1 2 1
2 2 1
2 3 1
EOF
for M in 2 3; do
    CWS80=$(mktemp --suffix=.nz)
    assert_output "$BIN_FORK method=$M finH=$ZC80 wmax=3 outC=$CWS80 debug=1 threads=2" 0 "^1 1 0$" \
        "exported 2 codewords"
    if [ "$(grep -c '^1 ' "$CWS80")" != "2" ]; then
        echo "  [FAIL] Test 80 (method=$M): expected 2 codewords of weight 1 in the exported file"
        FAILED=1
    fi
    rm -f "$CWS80"
done
rm -f "$ZC80"

# Test 81: outC exports only the collection window (here: the minimum weight found), not the heavier codewords kept
# for the min_hits statistic
echo "Running Test 81: outC exports only the codewords of the minimum weight found (dW=0)"
CWS81=$(mktemp --suffix=.nz)
OUT81=$($BIN_FORK method=1 finH=$EXAMPLES_DIR/c96H.mmx steps=50000 outC=$CWS81 debug=0 threads=4)
W81=$(echo "$OUT81" | awk '{print $2}')
if [ -z "$W81" ] || [ ! -s "$CWS81" ] || grep -v '^%' "$CWS81" | awk -v w="$W81" '$1 != w {bad=1} END {exit !bad}'; then
    echo "  [FAIL] Test 81: stdout '$OUT81', weights: $(grep -v '^%' "$CWS81" | awk '{print $1}' | sort | uniq -c)"
    FAILED=1
fi
rm -f "$CWS81"

# Test 82: debug bits do not change the algorithm: the same RW result with debug=1 and debug=33 (one thread, fixed seed)
echo "Running Test 82: debug does not change the RW result"
O82A=$($BIN_FORK method=1 finH=$EXAMPLES_DIR/QX150.mtx finG=$EXAMPLES_DIR/QZ150.mtx steps=50000 seed=7 threads=1 \
    debug=1 2>/dev/null)
O82B=$($BIN_FORK method=1 finH=$EXAMPLES_DIR/QX150.mtx finG=$EXAMPLES_DIR/QZ150.mtx steps=50000 seed=7 threads=1 \
    debug=33 2>/dev/null)
if [ -z "$O82A" ] || [ "$O82A" != "$O82B" ]; then
    echo "  [FAIL] Test 82: debug=1 gives '$O82A', debug=33 gives '$O82B'"
    FAILED=1
fi

# Test 83: warning for strongly non-uniform RW hit counts (surface code: some minimum-weight logicals are found much
# more often than others by the single columns, pairs=0), with the information-set estimate; no warning for the code
# QX150 (one thread, fixed seed)
echo "Running Test 83: hit-count warning with the information-set estimate"
OUT83=$(mktemp)
ERR83=$(mktemp)
$BIN_FORK method=1 $S5 steps=10000 seed=1 threads=1 pairs=0 debug=1 > "$OUT83" 2> "$ERR83"
if ! grep -q -E "^1 5 10000$" "$OUT83" || \
   ! grep -q "Warning: non-uniform RW hit counts of the 100 codewords of weight 5" "$ERR83" || \
   ! grep -q "information-set estimate (uniform random information sets, n=1958, rank(H)=120)" "$ERR83" || \
   ! grep -q -E "^# RW information sets: n=1958, rank\(H\)=120, steps=10000 \(uniform permutations: 5000\)" "$ERR83"
then
    echo "  [FAIL] Test 83: stdout '$(cat "$OUT83")', stderr '$(cat "$ERR83")'"
    FAILED=1
fi
rm -f "$OUT83" "$ERR83"
echo "Running Test 83: no hit-count warning for QX150"
if $BIN_FORK method=1 finH=$EXAMPLES_DIR/QX150.mtx finG=$EXAMPLES_DIR/QZ150.mtx steps=50000 seed=7 threads=1 \
    debug=1 2>&1 | grep -q "non-uniform"; then
    echo "  [FAIL] Test 83: unexpected hit-count warning for QX150"
    FAILED=1
fi

# Debug bits (see --morehelp): 4 status, 8 timing, 16 code parameters, 32 codewords, 64 arguments, 128 matrices,
# 256 codeword dump; with debug=1, the reason why the run ended

# Test 84: periodic status (debug=4) after 1 s, without the summary lines
echo "Running Test 84: periodic status line (debug=4)"
OUT84=$(mktemp)
ERR84=$(mktemp)
$BIN_FORK method=1 $S5 steps=100000000 min_hits=0 timeout=1.5 threads=2 debug=4 > "$OUT84" 2> "$ERR84"
if ! grep -q -E "^1 5 [0-9]+$" "$OUT84" || \
   ! grep -q -E "^# status 1\.[0-9]s: bounds \[1, 5\]; RW: [0-9]+ of 100000000 steps" "$ERR84" || \
   grep -q -E "^# (input|stopped)" "$ERR84"; then
    echo "  [FAIL] Test 84: stdout '$(cat "$OUT84")', stderr '$(cat "$ERR84")'"
    FAILED=1
fi
rm -f "$OUT84" "$ERR84"

# Test 85: code parameters (debug=16): ranks of H and L (or G), and k
assert_output "$BIN_FORK method=2 fdem=$EXAMPLES_DIR/surf_d3.dem drop=0 wmax=3 threads=2 debug=16" 0 "^3 3 0$" \
    "^# code parameters: n=221, rank\(H\)=24, dim ker\(H\)=197, rank\(L\)=1, k=rank\(\[H;L\]\)-rank\(H\)=1 "
assert_output "$BIN_FORK method=2 finH=$EXAMPLES_DIR/QX150.mtx finG=$EXAMPLES_DIR/QZ150.mtx wmax=2 threads=2 debug=16" \
    0 "^3 0 0$" "rank\(H\)=59, dim ker\(H\)=91, rank\(G\)=59, k=n-rank\(H\)-rank\(G\)=32 "

# Test 86: codeword supports (debug=32), matrices (debug=128), and the codeword dump (debug=256) for the classical
# [5,1] repetition code (one thread, fixed seed)
REP86=$(mktemp --suffix=.mtx)
cat << 'EOF' > "$REP86"
%%MatrixMarket matrix coordinate integer general
4 5 8
1 1 1
1 2 1
2 2 1
2 3 1
3 3 1
3 4 1
4 4 1
4 5 1
EOF
echo "Running Test 86: codeword supports, matrices, and codeword dump (debug=0x1a0)"
OUT86=$(mktemp)
ERR86=$(mktemp)
$BIN_FORK method=1 finH=$REP86 steps=100 seed=3 threads=1 debug=0x1a0 > "$OUT86" 2> "$ERR86"
if ! grep -q -E "^1 5 [0-9]+$" "$OUT86" || \
   ! grep -q "RW codeword of weight 5 (1-based columns): 1 2 3 4 5$" "$ERR86" || \
   ! grep -q "^# matrix H: 4 x 5, 8 nonzeros" "$ERR86" || ! grep -q "^# \.\.11\.$" "$ERR86" || \
   ! grep -q -E "^# cw: \[ 1 2 3 4 5 \] cnt=([0-9]+) \(u=\1, w=0\)$" "$ERR86"; then
    echo "  [FAIL] Test 86: stdout '$(cat "$OUT86")', stderr '$(cat "$ERR86")'"
    FAILED=1
fi
rm -f "$REP86" "$OUT86" "$ERR86"

# Test 87: command-line arguments (debug=64) are echoed, also those before the debug argument
assert_output "$BIN_FORK wmax=3 method=2 fdem=$EXAMPLES_DIR/surf_d3.dem seed=5 threads=2 debug=64" 0 "^3 3 0$" \
    "^# read wmax=3, wmax=3$"
assert_output "$BIN_FORK wmax=3 method=2 fdem=$EXAMPLES_DIR/surf_d3.dem seed=5 threads=2 debug=64" 0 "^3 3 0$" \
    "^# initializing rng from seed=5$"

# Test 88: the reason why the run ended (debug=1)
assert_output "$BIN_FORK method=1 $S5 steps=100000000 min_hits=0 timeout=0.3 threads=2 debug=1" 0 "^1 5 [0-9]+$" \
    "^# stopped after [0-9.]+s: timeout=0\.3s reached, [0-9]+ of 100000000 RW steps done$"
assert_output "$BIN_FORK method=1 $S5 steps=100 min_hits=0 threads=2 debug=1" 0 "^1 [0-9]+ 100$" \
    "^# stopped after [0-9.]+s: all 100 RW steps done$"
assert_output "$BIN_FORK method=2 $S5 wmax=3 threads=2 debug=1" 0 "^4 0 0$" \
    "^# stopped after [0-9.]+s: CC done up to w=3$"
assert_output "$BIN_FORK method=1 $S5 wmin=5 steps=100000 threads=2 debug=1" 0 "^1 [1-5] [0-9]+$" \
    "^# stopped after [0-9.]+s: found a codeword of weight [1-5] <= wmin=5$"
assert_output "$BIN_FORK method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=3 threads=2 debug=1" 0 "^3 3 0$" \
    "^# stopped after [0-9.]+s: CC found min-weight codeword: d=3$"
CWS88=$(mktemp --suffix=.nz)
assert_output "$BIN_FORK method=2 fdem=$EXAMPLES_DIR/surf_d3.dem dW=1 wmax=4 outC=$CWS88 threads=2 debug=1" \
    0 "^3 3 0$" "^# stopped after [0-9.]+s: CC enumerated all codewords of weight 3\.\.4 \(for outC\)$"
rm -f "$CWS88"
assert_output "$BIN_FORK method=3 $S5 dmin=5 dmax=5 threads=2 debug=1" 0 "^5 5 0$" \
    "^# stopped: bounds coincide: dmin = dmax = 5 \(supplied, no search\)$"

# Test 89: predicted vs measured CC work (debug=8)
assert_output "$BIN_FORK method=2 $S5 wmax=4 threads=2 debug=8" 0 "^5 0 0$" \
    "^# CC round w=4 work: [0-9.e+-]+ thread-s measured, [0-9.e+-]+ thread-s predicted"

# Test 90: method=2 without wmax and timeout: a known upper bound dmax (or a stop target dstop) suffices (CC ends once
# dmin=dmax, or dmin=dstop), otherwise an error
assert_output "$BIN_FORK method=2 $S5 dmax=5 timeout=0 threads=2 debug=0" 0 "^5 5 0$" ""
assert_output "$BIN_FORK method=2 $S5 timeout=0 threads=2 debug=0" 255 "" \
    "either parameter wmax>0, dmax>0, dstop>0, or timeout>0 should be specified for CC method=2"

# Test 91: logical check L c != 0 with the column bit masks of L for more than 64 logical operators ([[900,182,8]]:
# three words per column), in RW and in CC
Q900="finH=$EXAMPLES_DIR/QX900.mtx finG=$EXAMPLES_DIR/QZ900.mtx"
assert_output "$BIN_FORK method=1 $Q900 steps=2000 min_hits=0 seed=11 threads=1 debug=0" 0 "^1 8 2000$" ""
assert_output "$BIN_FORK method=2 $Q900 wmax=8 threads=2 debug=0" 0 "^8 8 0$" ""

# Test 92: the experimental ksub prints a warning regardless of debug, without and with min_hits (c1920H: ksub is
# used, m >= nu), and also when it falls back to full-matrix RW (surf_d3.dem: m < nu); no warning without RW
KSUB_WARN="^# WARNING: ksub=[0-9]+ is experimental and should not be used"
assert_output "$BIN_FORK method=1 finH=$EXAMPLES_DIR/c1920H.mmx ksub=32 steps=100 min_hits=0 threads=2 debug=0" \
    0 "^1 [0-9]+ 100$" "$KSUB_WARN"
assert_output "$BIN_FORK method=1 finH=$EXAMPLES_DIR/c1920H.mmx ksub=16 steps=100 threads=2 debug=0" \
    0 "^1 [0-9]+ 100$" "$KSUB_WARN"
assert_output "$BIN_FORK method=1 fdem=$EXAMPLES_DIR/surf_d3.dem ksub=32 steps=200 threads=2 debug=0" \
    0 "^1 3 [0-9]+$" "$KSUB_WARN"
echo "Running Test 92: no ksub warning with method=2"
if $BIN_FORK method=2 fdem=$EXAMPLES_DIR/surf_d3.dem ksub=32 wmax=3 threads=2 debug=0 2>&1 | \
    grep -q "is experimental"; then
    echo "  [FAIL] Test 92: unexpected ksub warning with method=2"
    FAILED=1
fi

# Test 93: space form of the file arguments ("fin= name", "finH= name", ...); a missing file name is an error
assert_output "$BIN_FORK method=2 fin= $EXAMPLES_DIR/try wmax=8 threads=2 debug=0" 0 "^4 4 0$" ""
assert_output \
    "$BIN_FORK method=2 finH= $EXAMPLES_DIR/tryX.mtx finG= $EXAMPLES_DIR/tryZ.mtx wmax=8 threads=2 debug=0" \
    0 "^4 4 0$" ""
assert_output "$BIN_FORK method=2 wmax=8 threads=2 debug=0 fin=" 255 "" "'fin=': missing file name"

# Test 94: hybrid RW (n >= 500): uniform and window permutations alternate with the global step index, i.e., half of
# the steps are uniform also for batches of one step (chunk_size=1)
assert_output "$BIN_FORK method=1 $S5 steps=1000 chunk_size=1 min_hits=0 threads=2 debug=1" 0 "^1 [0-9]+ 1000$" \
    "^# RW information sets: n=1958, rank\(H\)=120, steps=1000 \(uniform permutations: 500\)"

# Test 95: stop target dstop (not an upper bound): method=3 ends once CC certifies dmin=dstop, without waiting for RW
# (which would run until the timeout here); dmax is the lightest codeword found (0: none), never dstop; method=2 needs
# no wmax; a lighter codeword gives the exact distance; a supplied dmin >= dstop needs no search; dstop >= 0
assert_output "$BIN_FORK method=3 $S5 dstop=4 steps=100000000 min_hits=0 timeout=60 threads=4 debug=1" \
    0 "^4 (0|5) [0-9]+$" "^# stopped after [0-9.]+s: lower bound dmin=4 reached dstop=4$"
assert_output "$BIN_FORK method=2 $S5 dstop=4 threads=2 debug=1" 0 "^4 0 0$" \
    "^# stopped after [0-9.]+s: lower bound dmin=4 reached dstop=4$"
assert_output "$BIN_FORK method=3 $S5 dstop=7 threads=2 debug=0" 0 "^5 5 [0-9]+$" ""
assert_output "$BIN_FORK method=1 $S5 dmin=4 dstop=4 debug=1" 0 "^4 0 0$" \
    "^# stopped: lower bound dmin=4 >= dstop=4 \(supplied, no search\)$"
assert_output "$BIN_FORK method=3 $S5 dstop=-1 debug=0" 255 "" "invalid 'dstop=-1': dstop must be non-negative"
# with outC, the codewords of weight up to dstop+dW are collected (here d=5: the round w=6 is not needed)
CWS95=$(mktemp --suffix=.nz)
for m in 2 3; do
    assert_output "$BIN_FORK method=$m $S5 dstop=4 dW=1 outC=$CWS95 threads=4 debug=1" 0 "^5 5 0$" \
        "^# stopped after [0-9.]+s: CC enumerated all codewords of weight 5 \(for outC\)$"
done
rm -f "$CWS95"

# Test 96: conflicting inputs: fin followed by finH/finG/finL is replaced by them (with a warning); fdem together with
# matrix files, and pmin without fdem, are errors
assert_output "$BIN_FORK method=2 fin=$EXAMPLES_DIR/try $S5 wmax=5 threads=2 debug=0" 0 "^5 5 0$" \
    "^# Warning: finH='.*surf_d5_H.mmx' given: ignoring fin='.*try'$"
assert_output "$BIN_FORK method=2 fdem=$EXAMPLES_DIR/surf_d3.dem finH=$EXAMPLES_DIR/surf_d5_H.mmx wmax=3 threads=2 \
debug=0" 255 "" "Cannot specify matrix files \(fin, finH, finG, finL\) along with fdem"
assert_output "$BIN_FORK method=2 $S5 pmin=0.01 wmax=5 threads=2 debug=0" 255 "" \
    "pmin can only be used when fdem is specified"

# Test 97: CC work units: with few start columns per thread, the start columns of a CC round are split into work units
# which also fix the first cluster extensions (here: the round at a supplied dmin, of unknown work).  The units find
# all codewords of weight 6 of the [[150,32,6]] code once each (500 codewords, 500 hits), as whole start columns do;
# also with the expert start list (59 codewords, 59 hits; with cc_prune=0, one of them, which contains two listed
# columns, is found from both: 60 hits) and in method=3
QX97="finH=$EXAMPLES_DIR/QX150.mtx finG=$EXAMPLES_DIR/QZ150.mtx"
CWS97=$(mktemp --suffix=.nz)
assert_output "$BIN_FORK method=2 $QX97 dmin=6 outC=$CWS97 threads=8 nothrottle=1 debug=8" 0 "^6 6 0$" \
    "^# searching w=6 with 8 CC threads, columns \[0, 144\], [0-9]+ work units per start column$"
for thr in "threads=1" "threads=8 nothrottle=1"; do
    assert_output "$BIN_FORK method=2 $QX97 dmin=6 outC=$CWS97 $thr debug=1" 0 "^6 6 0$" \
        "^# codewords accumulated: total=500, min_w=6: cws=500, total_hits=500,"
    assert_output "$BIN_FORK method=2 $QX97 dmin=6 start=0,7,77 outC=$CWS97 $thr debug=1" 0 "^6 6 0$" \
        "^# codewords accumulated: total=59, min_w=6: cws=59, total_hits=59,"
    assert_output "$BIN_FORK method=2 $QX97 dmin=6 start=0,7,77 outC=$CWS97 $thr cc_prune=0 debug=1" 0 "^6 6 0$" \
        "^# codewords accumulated: total=59, min_w=6: cws=59, total_hits=60,"
done
assert_output "$BIN_FORK method=3 $QX97 dmin=6 steps=0 outC=$CWS97 threads=8 nothrottle=1 debug=8" 0 "^6 6 0$" \
    "^# CC round w=6 started: 8 CC threads, .*, [0-9]+ work units per start column\)$"
assert_output "$BIN_FORK method=3 $QX97 dmin=6 steps=0 outC=$CWS97 threads=8 nothrottle=1 debug=1" 0 "^6 6 0$" \
    "^# codewords accumulated: total=500, min_w=6: cws=500, total_hits=500,"
rm -f "$CWS97"

# Test 98: thread allocation: method=2 starts 4 threads, and the others for a CC round predicted to take longer (all at
# once with nothrottle=1); the default number of threads is the number of CPUs available (CPU affinity mask)
assert_output "$BIN_FORK method=2 $S5 wmax=5 threads=8 debug=2" 0 "^5 5 0$" \
    "^# run plan: method=2 \(CC\), 8 threads \(4 at first\), CC from w=1 to wmax=5 \(cc_prune=1\),"
assert_output "$BIN_FORK method=2 $S5 wmax=5 threads=8 debug=8" 0 "^5 5 0$" \
    "^# CC round w=[345]: started 4 more threads, 8 in total \("
assert_output "$BIN_FORK method=2 $S5 wmax=5 threads=8 nothrottle=1 debug=2" 0 "^5 5 0$" \
    "^# run plan: method=2 \(CC\), 8 threads, CC from w=1 to wmax=5 \(cc_prune=1\),"
if command -v taskset > /dev/null && [ "$(nproc --all)" -ge 2 ]; then
    assert_output "taskset -c 0,1 $BIN_FORK method=1 $S5 steps=100 debug=2" 0 "^1 [0-9]+ [0-9]+$" \
        "^# run plan: method=1 \(RW\), 2 threads, steps=100,"
fi

# Test 99: RW elimination: the blocked elimination (Four Russians, from 128 rows; c1920H: 1280 rows), and the choice
# of each thread between it and the column-by-column elimination (128 to 512 rows; surf_d7: 336 rows) find the same
# information sets as the column-by-column elimination: the same codewords and hit counts with a fixed seed and one
# thread (reference values of the column-by-column elimination, with the single columns only: pairs=0)
RW99="steps=60 seed=5 threads=1 min_hits=1000000 cov_cws=5000 pairs=0"
assert_output "$BIN_FORK method=1 finH=$EXAMPLES_DIR/c1920H.mmx $RW99 debug=1" 0 "^1 344 60$" \
    "^# codewords accumulated: total=5151, min_w=344: cws=1, total_hits=1,"
S7="finH=$EXAMPLES_DIR/surf_d7_H.mmx finL=$EXAMPLES_DIR/surf_d7_L.mmx"
RW99="steps=400 seed=6 threads=1 min_hits=1000000 cov_cws=5000 pairs=0"
assert_output "$BIN_FORK method=1 $S7 $RW99 debug=1" 0 "^1 7 400$" \
    "^# codewords accumulated: total=5835, min_w=7: cws=70, total_hits=71,"
assert_output "$BIN_FORK method=1 $S7 $RW99 debug=8" 0 "^1 7 400$" \
    "^# \[thread 0\] RW elimination: (blocked|column by column) \(blocked: [0-9.]+x the time of column by column\)$"

# Test 100: CC with bitset syndromes and a hash lookup for the last column of a cluster (smax=0, the default; here the
# walk of the first unsatisfied row, cc_prune=0) finds the same codewords with the same hit counts as the recursion on
# sorted syndrome lists (smax>0), also with duplicate columns of H (surf_d5 circuit DEM: 1958 columns, 1679 distinct; a
# codeword of weight 5 is found up to 5 times), with the rows of H in the input order (sort_rows=0) and sorted by weight
# (the default for matrices: other hit counts)
CWS100=$(mktemp --suffix=.nz)
for sm in 0 3; do
    assert_output "$BIN_FORK method=2 $S5 wmax=5 smax=$sm outC=$CWS100 threads=8 nothrottle=1 sort_rows=0 cc_prune=0 \
debug=1" 0 "^5 5 0$" "^# codewords accumulated: total=18043, min_w=5: cws=18043, total_hits=18147,"
    assert_output "$BIN_FORK method=2 $S5 wmax=5 smax=$sm outC=$CWS100 threads=8 nothrottle=1 cc_prune=0 debug=1" 0 \
        "^5 5 0$" "^# codewords accumulated: total=18043, min_w=5: cws=18043, total_hits=18063,"
done
assert_output "$BIN_FORK method=2 $S5 wmax=3 threads=2 debug=8" 0 "^4 0 0$" \
    "^# CC: branch and bound \(cc_prune=1\), bitset syndromes, hash of 1958 nonzero columns \(1679 distinct\) for the"
rm -f "$CWS100"

# Test 101: the components of a decomposed DEM error (separated by `^`) are summed mod 2: a detector (here D1) or an
# observable (L0) listed twice in one error line cancels (before 0.12.4, D1 was kept: no codeword of weight 2 here)
DEM101=$(mktemp --suffix=.dem)
printf 'error(0.1) D0 D1 ^ D1 D2\nerror(0.1) D0 D2 L0\n' > "$DEM101"
assert_output "$BIN_FORK method=2 fdem=$DEM101 wmax=2 threads=1 debug=1" 0 "^2 2 0$" \
    "^#   H: 3 x 2 \(detectors\), 4 nonzeros, max row weight 2, max column weight 2$"
printf 'error(0.1) D0 L0 ^ D1 L0\nerror(0.1) D0 D1 L0\n' > "$DEM101"
assert_output "$BIN_FORK method=2 fdem=$DEM101 wmax=2 threads=1 debug=1" 0 "^2 2 0$" \
    "^#   L: 1 x 2 \(observables\), 1 nonzeros, max row weight 1, max column weight 1$"
rm -f "$DEM101"

# Test 102: the equal columns of a DEM (same detectors and observables) are merged (merge=1, the default; exact): the
# column indices (outC, finC, and the printed codewords) are those of the DEM, the codewords found use the first one of
# each set of equal columns (here 1-based columns 1=5, 2=4, and 6=7), and two equal columns of a finC codeword cancel;
# merge=0 keeps all columns (4 codewords of weight 3)
DEM102=$(mktemp --suffix=.dem)
CWS102=$(mktemp --suffix=.nz)
printf 'error(0.1) D0 L0\nerror(0.1) D0 D1\nerror(0.1) D1\nerror(0.1) D0 D1\n' > "$DEM102"
printf 'error(0.1) D0 L0\nerror(0.1) D2\nerror(0.1) D2\n' >> "$DEM102"
for bin in "$BIN_FORK" "$BIN"; do
    assert_output "$bin method=2 fdem=$DEM102 wmax=3 outC=$CWS102 threads=1 debug=1" 0 "^3( 3 0)?$" \
        "^#   merged 3 equal columns of the DEM: 4 of 7 columns kept \(merge=0 to keep all\)$"
    if [ "$(grep -v '^%' "$CWS102" | tr -s ' ')" != "3 1 2 3" ]; then
        echo "  [FAIL] Test 102: expected the single codeword '3 1 2 3' in outC, got: $(grep -v '^%' "$CWS102")"
        FAILED=1
    fi
done
assert_output "$BIN_FORK method=2 fdem=$DEM102 wmax=3 merge=0 drop=0 outC=$CWS102 threads=1 debug=1" 0 "^3 3 0$" \
    "^#   H: 3 x 7 \(detectors\), 9 nonzeros"
if [ "$(grep -c -v '^%' "$CWS102")" != 4 ]; then
    echo "  [FAIL] Test 102: expected 4 codewords in outC with merge=0, got: $(grep -v '^%' "$CWS102")"
    FAILED=1
fi
printf '%%%% NZLIST\n3 3 4 5\n5 1 2 3 6 7\n' > "$CWS102"
assert_output "$BIN_FORK method=2 fdem=$DEM102 wmax=3 finC=$CWS102 threads=1 debug=33" 0 "^3 3 0$" \
    "^# read 1 codewords from .*, total 1, min weight 3$"
assert_output "$BIN_FORK method=2 fdem=$DEM102 wmax=3 finC=$CWS102 threads=1 debug=33" 0 "^3 3 0$" \
    "^# finC: lightest codeword of weight 3 \(1-based columns\): 1 2 3$"
assert_output "$BIN_FORK method=2 fdem=$DEM102 wmax=3 merge=2" 255 "" "invalid 'merge=2': merge should be 0 or 1"
rm -f "$DEM102" "$CWS102"

# Test 103: drop=FILE drops the listed detectors (here D2, of the other basis) if this is exact: every column touching
# them has a twin with the same other detectors and no dropped one (1-based columns 4 -> 1 and 6 -> 3; column 5 has no
# other detector), with the same observables (first DEM) or L + A H (second DEM: A maps D2 to L0). The codewords found
# (outC) are those of the full DEM, and finC codewords are mapped to the twins. Otherwise, a note is printed and the
# full DEM is used (third DEM: column 7 has no twin; fourth DEM: the observables are not L + A H)
DEM103=$(mktemp --suffix=.dem)
DROP103=$(mktemp --suffix=.txt)
CWS103=$(mktemp --suffix=.nz)
printf '# the other basis\n2\n' > "$DROP103"
printf 'error(0.1) D0 L0\nerror(0.1) D0 D1\nerror(0.1) D1\n' > "$DEM103"
printf 'error(0.1) D0 D2 L0\nerror(0.1) D2\nerror(0.1) D1 D2\n' >> "$DEM103"
for bin in "$BIN_FORK" "$BIN"; do
    assert_output "$bin method=2 fdem=$DEM103 drop=$DROP103 wmax=3 outC=$CWS103 threads=1 debug=1" 0 "^3( 3 0)?$" \
        "^#   dropped 1 of 3 detectors \(drop=.*\), exact: the 3 columns touching them have twins \(the same observ"
    if [ "$(grep -v '^%' "$CWS103" | tr -s ' ')" != "3 1 2 3" ]; then
        echo "  [FAIL] Test 103: expected the single codeword '3 1 2 3' in outC, got: $(grep -v '^%' "$CWS103")"
        FAILED=1
    fi
done
printf '%%%% NZLIST\n3 2 4 6\n4 2 3 4 5\n' > "$CWS103"
assert_output "$BIN_FORK method=2 fdem=$DEM103 drop=$DROP103 wmax=3 finC=$CWS103 threads=1 debug=33" 0 "^3 3 0$" \
    "^# read 1 codewords from .*, total 1, min weight 3$"
assert_output "$BIN_FORK method=2 fdem=$DEM103 drop=$DROP103 wmax=3 start=4 threads=1" 255 "" "no start column left"
assert_output "$BIN_FORK method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx drop=$DROP103 wmax=3" 255 "" \
    "drop can only be used when fdem is specified"
printf 'error(0.1) D0 L0\nerror(0.1) D0 D1\nerror(0.1) D1\n' > "$DEM103"
printf 'error(0.1) D0 D2\nerror(0.1) D2 L0\nerror(0.1) D1 D2 L0\n' >> "$DEM103"
assert_output "$BIN_FORK method=2 fdem=$DEM103 drop=$DROP103 wmax=3 threads=1 debug=1" 0 "^3 3 0$" \
    "twins \(L \+ A H with A != 0\), 3 columns kept$"
printf 'error(0.1) D0 L0\nerror(0.1) D0 D1\nerror(0.1) D1\n' > "$DEM103"
printf 'error(0.1) D0 D2 L0\nerror(0.1) D2\nerror(0.1) D1 D2\nerror(0.1) D2 D3\n' >> "$DEM103"
assert_output "$BIN_FORK method=2 fdem=$DEM103 drop=$DROP103 wmax=3 threads=1 debug=0" 0 "^3 3 0$" \
    "^# Note: drop=.* not used: column 7 \(1-based\) touches a dropped detector, but no column with the same other"
printf 'error(0.1) D0 L0\nerror(0.1) D0 D1\nerror(0.1) D1\n' > "$DEM103"
printf 'error(0.1) D0 D2\nerror(0.1) D2\nerror(0.1) D1 D2\n' >> "$DEM103"
assert_output "$BIN_FORK method=2 fdem=$DEM103 drop=$DROP103 wmax=3 threads=1 debug=0" 0 "^3 3 0$" \
    "^# Note: drop=.* not used: the observables of the twins are not L \+ A H for any matrix A; using the full DEM$"
rm -f "$DEM103" "$DROP103" "$CWS103"

# Test 104: --help TOPIC and --morehelp TOPIC print only the entries of the full help for TOPIC (names and aliases,
# without case, '-', '_', and '=VALUE'), in the order of the full help under their section headers (each once); a
# group of the section "Option groups" (also as singular or plural) gives all its entries, and a parameter of a group
# also the other entries of its group and the group (see also); an unknown topic gives the closest names on stderr
# and exit code 1; an argument with '=' is not a topic
assert_output "$BIN_FORK --help maxC" 0 "^  maxC=\[int\] " ""
assert_output "$BIN_FORK --morehelp=dest" 0 "^  dexp=\[int\] " ""
assert_output "$BIN_FORK --help no-throttle" 0 "^  nothrottle=\[0\|1\] " ""
assert_output "$BIN_FORK -h Batch" 0 "^  chunk_size=\[int\] " ""
assert_output "$BIN_FORK --morehelp debug" 0 "^    2048 : No size cutoff" ""
assert_output "$BIN_FORK --help mthod" 1 "" "no help entry for 'mthod'; did you mean method\?"
assert_output "$BIN_FORK --help windos" 1 "" "no help entry for 'windos'; did you mean windows"
assert_output "$BIN_FORK --help maxC,xyz" 1 "^  maxC=\[int\] " "no help entry for 'xyz' \(--morehelp lists all"
assert_output "$BIN_FORK --help finH=x.mtx" 0 "^Extra parameters" ""
assert_output "$BIN --help maxC" 0 "^  maxC=\[int\] " ""
help_outline() {  # the section headers, the first words of the entries, and '-' for a blank line
    awk '/^$/ {print "-"; next} /^[^ ]/ {print; next} /^  [^ ]/ {print $1}'
}
GROUPS104="Option groups (--help GROUP prints all entries of GROUP):"
GOT104=$("$BIN_FORK" --help method,SEED | help_outline)
EXP104=$(printf '%s\n' "Calculation method:" "method=[int]" "-" "Execution, multithreading, and timing:" "seed=[int]")
echo "Running: $BIN_FORK --help method,SEED (entries in the order of the full help)"
if [ "$GOT104" != "$EXP104" ]; then
    echo "  [FAIL] Expected: $EXP104"
    echo "         Got: $GOT104"
    FAILED=1
fi
GOT104=$("$BIN_FORK" --help maxC | help_outline)
EXP104=$(printf '%s\n' "Codeword collection and export:" "outC=[file]" "finC=[file]" "maxC=[int]" "dW=[int]" \
    "irred=[0|1]" "-" "$GROUPS104" "codewords")
echo "Running: $BIN_FORK --help maxC (the codewords group)"
if [ "$GOT104" != "$EXP104" ]; then
    echo "  [FAIL] Got: $GOT104"
    FAILED=1
fi
EXP104=$(printf '%s\n' "Execution, multithreading, and timing:" "kwin=[int]" "win_mode=[0|1]" "win_hyb=[int]" "-" \
    "$GROUPS104" "windows")
for topic in windows window WIN_HYB kwin; do
    echo "Running: $BIN_FORK -h $topic (the windows group)"
    if [ "$("$BIN_FORK" -h $topic | help_outline)" != "$EXP104" ]; then
        echo "  [FAIL] Got: $("$BIN_FORK" -h $topic | help_outline)"
        FAILED=1
    fi
done
GOT104=$("$BIN_FORK" --help experimental | help_outline)
echo "Running: $BIN_FORK --help experimental"
if [ "$GOT104" != "$(printf '%s\n' "Execution, multithreading, and timing:" "ksub=[int]" "refresh=[int]" "-" \
    "$GROUPS104" "experimental")" ]; then
    echo "  [FAIL] Got: $GOT104"
    FAILED=1
fi
GOT104=$("$BIN_FORK" --morehelp debug)
echo "Running: $BIN_FORK --morehelp debug (the whole debug section)"
if [ "$(printf '%s\n' "$GOT104" | head -1)" != "Debug output bitmap (debug=[int], default: 3):" ] || \
   [ "$(printf '%s\n' "$GOT104" | tail -1)" != "  the thread allocation and timing." ] || \
   [ "$(printf '%s\n' "$GOT104" | grep -c '^[^ ]')" != "1" ]; then
    echo "  [FAIL] Got: $GOT104"
    FAILED=1
fi

# Test 105: drop=auto (the default for a DEM; also drop=1) finds the detectors which can be dropped exactly from the
# DEM alone: the 8 X-type (minority-basis) detectors of the surface code (outD=FILE writes them); drop=0 keeps all
# detectors; a group which cannot be dropped is kept, and another one is dropped (here D4); nothing to drop; errors
OUT105=$(mktemp --suffix=.txt)
DEM105=$(mktemp --suffix=.dem)
for drop in "" "drop=auto" "drop=1"; do
    assert_output "$BIN_FORK method=2 fdem=$EXAMPLES_DIR/surf_d3.dem $drop wmax=3 outD=$OUT105 threads=2 debug=1" 0 \
        "^3 3 0$" "^#   dropped 8 of 24 detectors \(drop=auto\), exact: the 166 columns touching them have twins"
    if [ "$(grep -v '^#' "$OUT105")" != "5 7 8 10 13 15 16 18" ] || \
       [ "$(head -1 "$OUT105")" != "# 8 of 24 detectors dropped (drop=auto), 0-based" ]; then
        echo "  [FAIL] Test 105: outD: $(cat "$OUT105")"
        FAILED=1
    fi
done
assert_output "$BIN method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=3 debug=1" 0 "^3$" "dropped 8 of 24 detectors"
assert_output "$BIN_FORK method=2 fdem=$EXAMPLES_DIR/surf_d3.dem drop=0 wmax=3 outD=$OUT105 threads=2 debug=1" 0 \
    "^3 3 0$" "^#   H: 24 x 221 \(detectors\)"
if [ "$(cat "$OUT105")" != "# 0 of 24 detectors dropped (drop=0), 0-based" ]; then
    echo "  [FAIL] Test 105: outD with drop=0: $(cat "$OUT105")"
    FAILED=1
fi
printf 'error(0.1) D0 D1\nerror(0.1) D0 D2 D3\nerror(0.1) D2\nerror(0.1) D3 L0\nerror(0.1) D3 D5\nerror(0.1) D5\n' \
    > "$DEM105"
printf 'error(0.1) D4\nerror(0.1) D3 D4 L0\n' >> "$DEM105"
assert_output "$BIN_FORK method=2 fdem=$DEM105 wmax=4 outD=$OUT105 threads=1 debug=41" 0 "^3 3 0$" \
    "^# drop=auto: 1 of 6 detectors found \(4 groups to try: 1 dropped, 3 kept, 0 not verified;"
assert_output "$BIN_FORK method=2 fdem=$DEM105 wmax=4 threads=1 debug=33" 0 "^3 3 0$" \
    "CC codeword of weight 3 \(1-based columns\): 4 5 6$"
if [ "$(grep -v '^#' "$OUT105")" != "4" ]; then
    echo "  [FAIL] Test 105: outD: $(cat "$OUT105")"
    FAILED=1
fi
printf 'error(0.1) D0 L0\nerror(0.1) D0 D1\nerror(0.1) D1\n' > "$DEM105"
assert_output "$BIN_FORK method=2 fdem=$DEM105 wmax=3 threads=1 debug=1" 0 "^3 3 0$" \
    "^#   drop=auto: no detector can be dropped exactly \(the full DEM is used\)$"
assert_output "$BIN_FORK method=2 fdem=$EXAMPLES_DIR/surf_d3.dem classical=1 wmax=3 debug=1" 0 "^3 3 0$" \
    "^#   H: 24 x 221 \(detectors\)"
assert_output "$BIN_FORK method=2 fdem=$EXAMPLES_DIR/surf_d3.dem drop=$DEM105.none wmax=3" 255 "" \
    "cannot read this file .*; use drop=0, drop=1 or drop=auto, or a file"
assert_output "$BIN_FORK method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx drop=auto wmax=3" 255 "" \
    "drop can only be used when fdem is specified"
assert_output "$BIN_FORK method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx outD=$OUT105 wmax=3" 255 "" \
    "outD can only be used when fdem is specified"
assert_output "$BIN_FORK method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx drop=0 wmax=5" 0 \
    "^5 5 0$" ""
rm -f "$OUT105" "$DEM105"

# Test 106: sort_rows (default: 1 for matrices, 0 for a DEM) sorts the rows of H by increasing weight (stable) for CC;
# the input summary notes it.  The distance, the codewords enumerated (outC), and RW with the same seed do not depend
# on the order of the rows; rows of equal weight are left as they are (QX150); sort_rows=1 also sorts the detectors of
# a DEM; invalid values are an error
S5="finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx"
SORTED106="^#   rows of H sorted by increasing weight for CC \(sort_rows=1; sort_rows=0 keeps their order\)$"
CWS106A=$(mktemp --suffix=.nz)
CWS106B=$(mktemp --suffix=.nz)
for bin in "$BIN_FORK" "$BIN"; do
    assert_output "$bin method=2 $S5 wmax=5 threads=2 debug=1" 0 "^5( 5 0)?$" "$SORTED106"
    if $bin method=2 $S5 wmax=5 sort_rows=0 threads=2 debug=1 2>&1 | grep -q "rows of H sorted"; then
        echo "  [FAIL] Test 106: rows sorted with sort_rows=0 ($bin)"
        FAILED=1
    fi
done
assert_output "$BIN_FORK method=2 $S5 wmax=5 outC=$CWS106A threads=2 debug=0" 0 "^5 5 0$" ""
assert_output "$BIN_FORK method=2 $S5 wmax=5 outC=$CWS106B sort_rows=0 threads=2 debug=0" 0 "^5 5 0$" ""
if [ "$(grep -v '^%' "$CWS106A" | sort)" != "$(grep -v '^%' "$CWS106B" | sort)" ] || \
   [ "$(grep -c -v '^%' "$CWS106A")" -lt 2 ]; then
    echo "  [FAIL] Test 106: different CC codewords with sort_rows=0 and 1:"
    echo "         $(grep -v '^%' "$CWS106A" | sort | tr '\n' ';')"
    echo "         $(grep -v '^%' "$CWS106B" | sort | tr '\n' ';')"
    FAILED=1
fi
RW106A=$($BIN_FORK method=1 $S5 steps=3000 min_hits=0 seed=7 threads=1 outC=$CWS106A debug=0 2>/dev/null)
RW106B=$($BIN_FORK method=1 $S5 steps=3000 min_hits=0 seed=7 threads=1 outC=$CWS106B sort_rows=0 debug=0 2>/dev/null)
if [ "$RW106A" != "$RW106B" ] || [ "$(cat "$CWS106A")" != "$(cat "$CWS106B")" ]; then
    echo "  [FAIL] Test 106: RW results differ with sort_rows=0 and 1: '$RW106A' vs '$RW106B'"
    FAILED=1
fi
if $BIN_FORK method=2 finH=$EXAMPLES_DIR/QX150.mtx finG=$EXAMPLES_DIR/QZ150.mtx wmax=3 debug=1 2>&1 | \
   grep -q "rows of H sorted"; then
    echo "  [FAIL] Test 106: rows of equal weight reported as sorted (QX150)"
    FAILED=1
fi
if $BIN_FORK method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=3 threads=2 debug=1 2>&1 | grep -q "rows of H sorted"; then
    echo "  [FAIL] Test 106: the detectors of a DEM sorted by default"
    FAILED=1
fi
assert_output "$BIN_FORK method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=3 sort_rows=1 threads=2 debug=1" 0 "^3 3 0$" \
    "$SORTED106"
assert_output "$BIN_FORK method=2 $S5 wmax=3 sort_rows=2" 255 "" "invalid 'sort_rows=2': sort_rows should be 0 or 1"
rm -f "$CWS106A" "$CWS106B"

# Test 107: the hybrid RW windows (kwin=0, n >= 500: a window on every other step, uniform permutations on the others):
# win_hyb=0 (the default) gives 256 columns for a DEM and 512 for matrices, win_hyb=N the size (at most 3n/4); with
# debug&1, the input summary gives the rounds of a DEM (here a repetition code with 50 checks in 12 rounds: measurement
# errors flip two detectors D=50 apart, and a shift by D maps errors to errors; n=1162 columns, L = n D / m = 97
# columns per round; no round structure for the same DEM with randomly permuted detectors; none for matrices), and the
# end of RW the codewords found by the uniform and by the window steps; win_hyb is ignored with kwin > 0 and for
# n < 500 (warnings); invalid values are an error
DEM107=$(mktemp --suffix=.dem)
DEM107P=$(mktemp --suffix=.dem)
awk 'BEGIN { k = 50; R = 12
  for (t = 0; t < R; t++) {
    for (i = 0; i <= k; i++) { # data errors: qubit i flips checks i-1 and i (and L0 for i=0)
      s = "error(0.01)"
      if (i > 0) s = s " D" (t * k + i - 1)
      if (i < k) s = s " D" (t * k + i)
      if (i == 0) s = s " L0"
      print s
    }
    if (t < R - 1)
      for (i = 0; i < k; i++) # measurement errors: check i in rounds t and t+1
        print "error(0.01) D" (t * k + i) " D" ((t + 1) * k + i)
  }
}' > "$DEM107"
awk 'BEGIN { srand(1); m = 600; for (d = 0; d < m; d++) p[d] = d
  for (d = m - 1; d > 0; d--) { j = int(rand() * (d + 1)); x = p[d]; p[d] = p[j]; p[j] = x } }
  { s = $1; for (f = 2; f <= NF; f++) s = s " " (($f ~ /^D/) ? "D" p[substr($f, 2) + 0] : $f); print s }' \
    "$DEM107" > "$DEM107P"
RW107="method=1 steps=40 threads=2 debug=1"
WIN107="columns on every other step, uniform permutations on the others"
AUTO107="\(win_hyb=0: 256 for a DEM, 512 for matrices\)$"
STEP107="^# RW step types, w=[0-9]+: uniform steps: cws=[0-9]+, hits=[0-9]+; window steps: cws=[0-9]+, hits=[0-9]+;"
assert_output "$BIN_FORK fdem=$DEM107 $RW107" 0 "^1 [0-9]+ [0-9]+$" \
    "^#   rounds of the DEM: D=50 detectors and L=97 columns per round \(estimated: score 1\.000\)$"
assert_output "$BIN_FORK fdem=$DEM107 $RW107" 0 "" "^#   RW windows: 256 $WIN107 $AUTO107"
assert_output "$BIN_FORK fdem=$DEM107 $RW107" 0 "" "$STEP107 both: cws=[0-9]+ \(of [0-9]+ cws\)$"
assert_output "$BIN_FORK fdem=$DEM107P $RW107" 0 "" \
    "^#   rounds of the DEM: no round structure found \(best score 0\.[0-4][0-9]* < 0\.5\)$"
assert_output "$BIN_FORK fdem=$DEM107 $RW107 win_hyb=2000" 0 "" "^#   RW windows: 871 $WIN107 \(win_hyb=2000\)$"
assert_output "$BIN_FORK $S5 $RW107" 0 "" "^#   RW windows: 512 $WIN107 $AUTO107"
ERR107=$($BIN_FORK fdem=$DEM107 method=2 wmax=1 threads=2 debug=1 2>&1)
if printf '%s\n' "$ERR107" | grep -q -E "RW windows|RW step types" || \
   ! printf '%s\n' "$ERR107" | grep -q "rounds of the DEM: D=50 detectors"; then
    echo "  [FAIL] Test 107: CC alone: $ERR107"
    FAILED=1
fi
if $BIN_FORK $S5 $RW107 2>&1 | grep -q "rounds of the DEM"; then
    echo "  [FAIL] Test 107: rounds reported for a matrix"
    FAILED=1
fi
# every codeword of RW was found by the uniform steps, the window steps, or both
STEPS107=$($BIN_FORK fdem=$DEM107 method=1 steps=200 cov_cws=50 threads=2 debug=1 2>&1 | grep "^# RW step types")
if [ -z "$STEPS107" ] || ! printf '%s\n' "$STEPS107" | awk '{
      u = $0; sub(/.*uniform steps: cws=/, "", u); w = $0; sub(/.*window steps: cws=/, "", w)
      b = $0; sub(/.*both: cws=/, "", b); t = $0; sub(/.*\(of /, "", t)
      if ((u + 0) + (w + 0) - (b + 0) != (t + 0)) bad = 1 } END { exit bad }'; then
    echo "  [FAIL] Test 107: RW step types do not add up: $STEPS107"
    FAILED=1
fi
assert_output "$BIN_FORK fdem=$DEM107 $RW107 kwin=64 win_hyb=100" 0 "" \
    "^# WARNING: win_hyb=100 is ignored with kwin=64 \(a window of kwin columns on every RW step\)$"
assert_output "$BIN_FORK fdem=$EXAMPLES_DIR/surf_d3.dem method=1 steps=20 win_hyb=100 threads=2 debug=0" 0 "" \
    "^# WARNING: win_hyb=100 is ignored for n=[0-9]+ < 500 \(uniform permutations in all RW steps\)$"
assert_output "$BIN_FORK fdem=$DEM107 $RW107 win_hyb=-1" 255 "" \
    "invalid 'win_hyb=-1': win_hyb should be 0 \(256 columns for a DEM, 512 for matrices\) or positive"
rm -f "$DEM107" "$DEM107P"

# Test 108: RW pairs of non-pivot columns (p=2; pairs=l: the pairs of columns which agree on l random pivot rows; -1:
# auto, the default; 0: off): the summary line gives l (auto: l=8 for the surf_d5 circuit DEM with n=1958, rank 120,
# and d=5; pairs=0 without pairs, also with the experimental ksub, with a warning), and the pairs checked per step and
# the codewords found by pairs and by single columns; the pairs find many more minimum-weight codewords (here, 15125 of
# the 18043 of weight 5 in 2000 steps, 1029 without pairs), all valid (read back with finC); the information-set
# estimate includes the pairs; invalid values are an error
RW108="method=1 steps=400 seed=1 threads=1 debug=1"
assert_output "$BIN_FORK $S5 $RW108" 0 "^1 5 400$" \
    "^# RW information sets: n=1958, rank\(H\)=120, steps=400 \(uniform permutations: 200\), ksub=0, pairs=8, "
assert_output "$BIN_FORK $S5 $RW108" 0 "" "^# RW pairs \(p=2, buckets of l=8 random pivot rows\): [0-9.e+]+ pairs \
checked per step; codewords recorded with weight up to the upper bound: [0-9]+ from pairs, [0-9]+ from single columns$"
assert_output "$BIN_FORK $S5 $RW108 pairs=12" 0 "^1 5 400$" "ksub=0, pairs=12, "
ERR108=$($BIN_FORK $S5 $RW108 pairs=0 2>&1)
if ! printf '%s\n' "$ERR108" | grep -q "ksub=0, pairs=0, " || printf '%s\n' "$ERR108" | grep -q "^# RW pairs"; then
    echo "  [FAIL] Test 108: pairs=0: $ERR108"
    FAILED=1
fi
assert_output "$BIN_FORK finH=$EXAMPLES_DIR/c1920H.mmx method=1 steps=20 ksub=32 pairs=8 threads=1 debug=1" 0 "" \
    "^# WARNING: pairs=8 is ignored with ksub=32 \(no RW pairs with the experimental ksub\)$"
CWS108=$(mktemp --suffix=.nz)
N108=""
for pr in 0 -1; do
    $BIN_FORK $S5 method=1 steps=2000 seed=1 threads=1 pairs=$pr outC=$CWS108 debug=0 > /dev/null 2>&1
    N108="$N108 $(grep -vc '^%' "$CWS108")"
    assert_output "$BIN_FORK $S5 method=1 steps=1 finC=$CWS108 cov_cws=0 pairs=0 threads=1 debug=1" 0 "^1 5 1$" \
        "^# read [0-9]+ codewords from .*, total [0-9]+, min weight 5$"
    if $BIN_FORK $S5 method=1 steps=1 finC=$CWS108 cov_cws=0 threads=1 debug=1 2>&1 | grep -q "invalid codewords"; then
        echo "  [FAIL] Test 108: invalid codewords exported (pairs=$pr)"
        FAILED=1
    fi
done
echo "Running Test 108: minimum-weight codewords exported without and with pairs:$N108"
if ! echo "$N108" | awk '{exit !($1 > 0 && $2 >= 5 * $1)}'; then
    echo "  [FAIL] Test 108: the pairs found too few codewords:$N108"
    FAILED=1
fi
rm -f "$CWS108"
assert_output "$BIN_FORK $S5 method=1 steps=3000 seed=1 threads=1 min_hits=1000 cov_cws=0 debug=1" 0 "" \
    "information-set estimate \(uniform random information sets with pairs \(l=8\), n=1958, rank\(H\)=120\)"
for pr in -2 65; do
    assert_output "$BIN_FORK $S5 $RW108 pairs=$pr" 255 "" \
        "invalid 'pairs=$pr': pairs should be -1 \(auto\), 0 \(off\), or 1 to 64 \(random pivot rows of the buckets\)"
done

# Test 109: irred=1 (the default) records only irreducible codewords, not the union of a codeword and a disjoint nonzero
# vector with zero syndrome; the surface-code DEM has trivial codewords of weight 3, so that reducible codewords appear
# from the weight d+3 = 6 on.  The CC rounds w = 3..6 export 21, 100, 356, and 1198 codewords (with irred=0, also the
# reducible ones found: 1300 of weight 6 with cc_prune=0, 1264 with the branch and bound); RW with the pairs (one
# thread, a fixed seed) finds the same 1675 codewords (with irred=0, also reducible sums of two disjoint circuits); a
# reducible finC codeword is replaced by an irreducible part (c96H: the union of two disjoint codewords of weight 6),
# also by dist_m4ri_old; irred is 0 or 1
CWS109=$(mktemp --suffix=.nz)
FIN109=$(mktemp --suffix=.nz)
cnt109() {  # the numbers of codewords of each weight in an .nz file, e.g., "3:21 4:100 "
    grep -v '^%' "$1" | awk '{n[$1]++} END {for (w in n) printf "%s:%d\n", w, n[w]}' | sort -n | tr '\n' ' '
}
assert_output "$BIN_FORK method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=6 dW=3 outC=$CWS109 threads=4 debug=1" 0 \
    "^3 3 0$" "^# reducible codewords skipped \(irred=1, .*\): [0-9]+ times in CC, 0 times in RW$"
echo "Running Test 109: the irreducible codewords of weight 3..6 found by CC: $(cnt109 "$CWS109")"
if [ "$(cnt109 "$CWS109")" != "3:21 4:100 5:356 6:1198 " ]; then
    echo "  [FAIL] Test 109: expected 3:21 4:100 5:356 6:1198"
    FAILED=1
fi
for pr in 0 1; do
    assert_output "$BIN_FORK method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=6 dW=3 outC=$CWS109 threads=4 irred=0 \
cc_prune=$pr debug=1" 0 "^3 3 0$" "irred=0 \(expert option\): reducible codewords are also recorded"
    if [ "$(cnt109 "$CWS109")" != "3:21 4:100 5:356 6:$([ $pr = 0 ] && echo 1300 || echo 1264) " ]; then
        echo "  [FAIL] Test 109: irred=0 cc_prune=$pr: $(cnt109 "$CWS109")"
        FAILED=1
    fi
done
for ir in 1 0; do
    assert_output "$BIN_FORK method=1 fdem=$EXAMPLES_DIR/surf_d3.dem steps=3000 min_hits=0 dW=3 threads=1 seed=5 \
irred=$ir outC=$CWS109 debug=1" 0 "^1 3 3000$" ""
    N109=$(grep -vc '^%' "$CWS109")
    echo "Running Test 109: RW with the pairs, irred=$ir: $N109 codewords"
    if { [ $ir = 1 ] && [ "$N109" != 1675 ]; } || { [ $ir = 0 ] && [ "$N109" -le 1675 ]; }; then
        echo "  [FAIL] Test 109: RW with irred=$ir exported $N109 codewords"
        FAILED=1
    fi
done
printf '%%%%NZLIST\n%% the union of two disjoint codewords of weight 6 of c96H\n%s\n' \
    "12 18 21 25 46 50 59 60 63 64 76 91 93" > "$FIN109"
assert_output "$BIN_FORK method=2 finH=$EXAMPLES_DIR/c96H.mmx classical=1 finC=$FIN109 threads=1 debug=33" 0 \
    "^6 6 0$" "min weight 6; 1 reducible codewords replaced by an irreducible part \(irred=1\)"
assert_output "$BIN_FORK method=2 finH=$EXAMPLES_DIR/c96H.mmx classical=1 finC=$FIN109 threads=1 irred=0 debug=33" 0 \
    "^6 6 0$" "^# read 1 codewords from .*, total 1, min weight 12$"
assert_output "$BIN method=2 finH=$EXAMPLES_DIR/c96H.mmx classical=1 finC=$FIN109 wmax=6 debug=1" 0 "" \
    "min weight 6; 1 reducible codewords replaced by an irreducible part \(irred=1\)"
assert_output "$BIN_FORK method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=3 irred=2" 255 "" \
    "invalid 'irred=2': irred should be 0 or 1"
rm -f "$CWS109" "$FIN109"

# Test 110: the CC search: cc_prune=0 (the walk of the first unsatisfied row), 1 (the default: the branch and bound),
# and 2 (also the cover bounds) find the same codewords and distances, and with cc_prune=1 or 2 each codeword exactly
# once (one CC hit): the surf_d5 circuit DEM as matrices (equal columns: with cc_prune=0, a codeword is found up to 5
# times), the CC rounds w = 3..6 of the surface-code DEM, c96H, and the work units (fixed levels) on 8 threads, also
# with the expert start list (see Test 97); noscan=1 dmin=wmax; cc_prune is 0, 1, or 2
CWS110=$(mktemp --suffix=.nz)
CWS110B=$(mktemp --suffix=.nz)
assert_output "$BIN_FORK method=2 $S5 wmax=5 outC=$CWS110 threads=8 nothrottle=1 cc_prune=0 debug=1" 0 "^5 5 0$" \
    "^# codewords accumulated: total=18043, min_w=5: cws=18043, total_hits=18063,"
grep -v '^%' "$CWS110" | sort > "$CWS110B"
for pr in 1 2; do
    assert_output "$BIN_FORK method=2 $S5 wmax=5 outC=$CWS110 threads=8 nothrottle=1 cc_prune=$pr debug=1" 0 \
        "^5 5 0$" "^# codewords accumulated: total=18043, min_w=5: cws=18043, total_hits=18043,"
    if ! cmp -s <(grep -v '^%' "$CWS110" | sort) "$CWS110B"; then
        echo "  [FAIL] Test 110: other codewords with cc_prune=$pr than with cc_prune=0 (surf_d5)"
        FAILED=1
    fi
    assert_output "$BIN_FORK method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=6 dW=3 outC=$CWS110 threads=4 cc_prune=$pr \
debug=1" 0 "^3 3 0$" ""
    if [ "$(cnt109 "$CWS110")" != "3:21 4:100 5:356 6:1198 " ]; then
        echo "  [FAIL] Test 110: surf_d3 with cc_prune=$pr: $(cnt109 "$CWS110")"
        FAILED=1
    fi
done
C96="finH=$EXAMPLES_DIR/c96H.mmx classical=1"
assert_output "$BIN_FORK method=2 $C96 wmax=8 dW=2 outC=$CWS110 threads=4 cc_prune=0 debug=0" 0 "^6 6 0$" ""
grep -v '^%' "$CWS110" | sort > "$CWS110B"
for pr in 1 2; do
    assert_output "$BIN_FORK method=2 $C96 wmax=8 dW=2 outC=$CWS110 threads=4 cc_prune=$pr debug=0" 0 "^6 6 0$" ""
    if ! cmp -s <(grep -v '^%' "$CWS110" | sort) "$CWS110B" || [ "$(wc -l < "$CWS110B")" != 27 ]; then
        echo "  [FAIL] Test 110: other codewords with cc_prune=$pr than with cc_prune=0 (c96H, 27 codewords)"
        FAILED=1
    fi
    assert_output "$BIN_FORK method=2 $QX97 dmin=6 outC=$CWS110 threads=8 nothrottle=1 cc_prune=$pr debug=9" 0 \
        "^6 6 0$" "^# codewords accumulated: total=500, min_w=6: cws=500, total_hits=500,"
    assert_output "$BIN_FORK method=2 $QX97 dmin=6 outC=$CWS110 threads=8 nothrottle=1 cc_prune=$pr debug=8" 0 \
        "^6 6 0$" "^# searching w=6 with 8 CC threads, columns \[0, 144\], [0-9]+ work units per start column$"
    assert_output "$BIN_FORK method=2 $QX97 dmin=6 start=0,7,77 outC=$CWS110 threads=8 nothrottle=1 cc_prune=$pr \
debug=1" 0 "^6 6 0$" "^# codewords accumulated: total=59, min_w=6: cws=59, total_hits=59,"
done
assert_output "$BIN_FORK method=2 $QX97 noscan=1 dmin=6 wmax=6 threads=4 debug=0" 0 "^6 6 0$" ""
assert_output "$BIN_FORK method=2 $QX97 dmin=6 threads=4 cc_prune=3" 255 "" \
    "invalid 'cc_prune=3': cc_prune should be 0, 1, or 2"
rm -f "$CWS110" "$CWS110B"

# Test 111: the number of nodes of a completed CC round (debug&8: the calls of the CC recursion) does not depend on the
# threads: a node at a fixed level of the work units, visited by all the units with the same earlier branches, is
# counted once (QX150 at the supplied dmin=6, completed with outC; 1 thread, 5 or 25 work units per start column, and
# method=3 without RW steps)
CWS111=$(mktemp --suffix=.nz)
for pr in 0 1; do
    N111=""
    for run in "method=2 threads=1" "method=2 threads=8 nothrottle=1" "method=2 threads=32 nothrottle=1" \
               "method=3 steps=0 threads=8 nothrottle=1"; do
        N111="$N111 $($BIN_FORK $run $QX97 dmin=6 outC=$CWS111 cc_prune=$pr debug=8 2>&1 | \
            sed -n 's/^# CC round w=6 work: .*, \([0-9][0-9]*\) nodes$/\1/p')"
    done
    echo "Running Test 111: the nodes of the CC round w=6 with cc_prune=$pr:$N111"
    if ! echo "$N111" | grep -q '^ \([0-9][0-9]*\) \1 \1 \1$'; then
        echo "  [FAIL] Test 111: the number of nodes depends on the threads (cc_prune=$pr):$N111"
        FAILED=1
    fi
done
rm -f "$CWS111"

if [ $FAILED -ne 0 ]; then
    echo "Some tests failed!"
    exit 1
else
    echo "All tests passed!"
    exit 0
fi



