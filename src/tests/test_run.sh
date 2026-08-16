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
assert_output "$BIN method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx wmax=5 debug=0" 0 "^5$" ""

# Test 2: CC noscan
assert_output "$BIN method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx wmax=5 noscan=1 debug=0" 0 "^5$" ""

# Test 3: CC fdem baseline
assert_output "$BIN method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=3 debug=0" 0 "^3$" ""

# Test 4: CC fdem with pmin
assert_output "$BIN method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=3 pmin=0.01 debug=0" 0 "^3$" ""

# Test 5: CC fdem with high pmin
assert_output "$BIN method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=3 pmin=0.02 debug=0" 0 "^-3$" ""

# Test 6: RW fdem baseline
assert_output "$BIN method=1 fdem=$EXAMPLES_DIR/surf_d3.dem steps=100 debug=0" 0 "^3$" ""

# Test 7: Validation: noscan=1 with method=1
assert_output "$BIN method=1 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx wmax=5 noscan=1 debug=0" 255 "" "noscan=1 only works with method=2"

# Test 8: Validation: noscan=1 with method=3
assert_output "$BIN method=3 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx wmax=5 noscan=1 debug=0" 255 "" "noscan=1 only works with method=2"

# Test 9: Validation: fdem and finH together
assert_output "$BIN method=2 fdem=$EXAMPLES_DIR/surf_d3.dem finH=$EXAMPLES_DIR/surf_d5_H.mmx wmax=3 debug=0" 255 "" "Cannot specify matrix files.*along with fdem"

# Test 10: Validation: pmin without fdem
assert_output "$BIN method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx wmax=3 pmin=0.01 debug=0" 255 "" "pmin can only be used when fdem is specified"

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
assert_output "$BIN debug=33 method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=3 finC=$TEMP_CWS1 outC=$TEMP_CWS2" 0 "" "read 128 codewords from"

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
assert_output "$BIN method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx wmax=3 classical=0 debug=0" 255 "" "L matrix.*is required for quantum code"

# Test 17: Conflict detection (classical=1 with finL)
assert_output "$BIN method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx wmax=3 classical=1 debug=0" 255 "" "Conflict: classical=1 specified"

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
assert_output "$BIN_FORK method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx wmax=5 debug=0 threads=4" 0 "^5 5 0$" ""

# Test 25: dist_m4ri method=2 (CC lower bound when distance > wmax, rw_steps=0)
assert_output "$BIN_FORK method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx wmax=3 debug=0 threads=4" 0 "^4 0 0$" ""

# Test 26: dist_m4ri method=1 (multithreaded RW, reported rw_steps > 0)
assert_output "$BIN_FORK method=1 fdem=$EXAMPLES_DIR/surf_d3.dem steps=100 debug=0 threads=4" 0 "^1 3 [0-9]+$" ""

# Test 27: dist_m4ri method=3 (bracketing mode with dexp)
assert_output "$BIN_FORK method=3 fdem=$EXAMPLES_DIR/surf_d3.dem dexp=3 timeout=10 debug=0 threads=4" 0 "^3 3 [0-9]+$" ""

# Test 28: dist_m4ri method=3 (bracketing mode with dest alias)
assert_output "$BIN_FORK method=3 fdem=$EXAMPLES_DIR/surf_d3.dem dest=3 timeout=10 debug=0 threads=4" 0 "^3 3 [0-9]+$" ""

# Test 29: dist_m4ri method=3 on surf_d5
assert_output "$BIN_FORK method=3 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx dexp=5 timeout=10 debug=0 threads=4 steps=500" 0 "^5 5 [0-9]+$" ""

# Test 30: dist_m4ri codeword saving and loading
TEMP_FORK_CWS1=$(mktemp --suffix=.nz)
TEMP_FORK_CWS2=$(mktemp --suffix=.nz)
assert_output "$BIN_FORK method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=3 outC=$TEMP_FORK_CWS1 debug=0 threads=4" 0 "^3 3 0$" ""
assert_output "$BIN_FORK method=3 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=3 finC=$TEMP_FORK_CWS1 outC=$TEMP_FORK_CWS2 debug=0 threads=4" 0 "^3 3 [0-9]+$" ""
if ! diff -q "$TEMP_FORK_CWS1" "$TEMP_FORK_CWS2" >/dev/null; then
    echo "  [FAIL] dist_m4ri codeword files differ"
    FAILED=1
fi
rm -f "$TEMP_FORK_CWS1" "$TEMP_FORK_CWS2"

# Test 31: dist_m4ri timeout graceful termination
assert_output "$BIN_FORK method=1 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx steps=10000000 timeout=0.5 debug=0 threads=4" 0 "^1 5 [0-9]+$" ""

# Test 32: dist_m4ri classical mode
assert_output "$BIN_FORK method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx wmax=3 debug=0 threads=4" 0 "^2 2 0$" ""

# Test 33: dist_m4ri method=2 with dW=1 and outC
TEMP_DW_CWS=$(mktemp --suffix=.nz)
assert_output "$BIN_FORK debug=15 method=2 fdem=$EXAMPLES_DIR/surf_d3.dem dW=1 wmax=4 outC=$TEMP_DW_CWS threads=4" 0 "^3 3 0$" "continuing up to w=4 for dW=1"
rm -f "$TEMP_DW_CWS"

# Test 34: dist_m4ri_old method=2 with dW=1 reporting
TEMP_M4RI_CWS=$(mktemp --suffix=.nz)
assert_output "$BIN debug=1 method=2 fdem=$EXAMPLES_DIR/surf_d3.dem dW=1 wmax=4 outC=$TEMP_M4RI_CWS" 0 "^3$" "CC round w=.*searched with dW=1"
rm -f "$TEMP_M4RI_CWS"

# Test 35: dist_m4ri method=2 timeout lower bound correctness
assert_output "$BIN_FORK method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx wmax=5 timeout=0.0001 debug=0 threads=4" 0 "^[1-5] 0 0$" ""

# Test 36: dist_m4ri method=3 timeout bound correctness
assert_output "$BIN_FORK method=3 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx wmax=5 timeout=0.0001 debug=0 threads=4" 0 "^[1-5] [0-5] [0-9]+$" ""

# Test 37: dmin parameter starting CC search directly at dmin
assert_output "$BIN_FORK method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx dmin=5 wmax=5 debug=0 threads=4" 0 "^5 5 0$" ""

# Test 38: dmax parameter in method 1
assert_output "$BIN_FORK method=1 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx dmax=5 steps=100 debug=0 threads=4" 0 "^1 5 [0-9]+$" ""

# Test 39: dmin and dmax in method 3
assert_output "$BIN_FORK method=3 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx dmin=4 dmax=5 timeout=5 debug=0 threads=4" 0 "^5 5 [0-9]+$" ""

# Test 40: conflicting debug parameters error
assert_output "$BIN debug=1 debug=2 method=2 fdem=$EXAMPLES_DIR/surf_d3.dem wmax=2" 255 "" "debug parameter specified multiple times with conflicting values"

if [ $FAILED -ne 0 ]; then
    echo "Some tests failed!"
    exit 1
else
    echo "All tests passed!"
    exit 0
fi


