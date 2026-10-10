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
assert_output "$BIN method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=3 outC=$TEMP_CWS1 debug=0" 0 "^3$" ""

# Step 2: Verify file 1 is not empty
if [ ! -s "$TEMP_CWS1" ]; then
    echo "  [FAIL] Codewords file is empty"
    FAILED=1
fi

# Step 3: Run again loading file 1 and saving to file 2, check for read message
assert_output "$BIN debug=33 method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=3 finC=$TEMP_CWS1 outC=$TEMP_CWS2" \
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
assert_output "$BIN_FORK --version" 0 "dist_m4ri version 0.11.0" ""

# Test 48: dist_m4ri_old --version
assert_output "$BIN --version" 0 "dist_m4ri version 0.11.0" ""

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
assert_output "$BIN_FORK debug=1 method=1 fdem=$EXAMPLES_DIR/surf_d3.dem steps=200 ksub=32 threads=4" \
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
assert_output "$BIN_FORK method=3 fdem=$EXAMPLES_DIR/surf_d3.dem outC=$CWS79 debug=1 threads=4" \
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
# more often than others), with the information-set estimate; no warning for the code QX150 (one thread, fixed seed)
echo "Running Test 83: hit-count warning with the information-set estimate"
OUT83=$(mktemp)
ERR83=$(mktemp)
$BIN_FORK method=1 $S5 steps=10000 seed=1 threads=1 debug=1 > "$OUT83" 2> "$ERR83"
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
assert_output "$BIN_FORK method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=3 threads=2 debug=16" 0 "^3 3 0$" \
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
   ! grep -q -E "^# cw: \[ 1 2 3 4 5 \] cnt=[0-9]+$" "$ERR86"; then
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

# Test 90: method=2 without wmax and timeout: a known upper bound dmax suffices (CC ends once dmin=dmax), otherwise
# an error
assert_output "$BIN_FORK method=2 $S5 dmax=5 timeout=0 threads=2 debug=0" 0 "^5 5 0$" ""
assert_output "$BIN_FORK method=2 $S5 timeout=0 threads=2 debug=0" 255 "" \
    "either parameter wmax>0, dmax>0, or timeout>0 should be specified for CC method=2"

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

if [ $FAILED -ne 0 ]; then
    echo "Some tests failed!"
    exit 1
else
    echo "All tests passed!"
    exit 0
fi



