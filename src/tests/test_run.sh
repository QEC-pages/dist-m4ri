#!/bin/bash

# Get directories
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
SRC_DIR="$SCRIPT_DIR/.."
WS_ROOT="$SCRIPT_DIR/../.."
BIN="$SRC_DIR/dist_m4ri"
EXAMPLES_DIR="$WS_ROOT/examples"

if [ ! -f "$BIN" ]; then
    echo "Error: dist_m4ri binary not found at $BIN"
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
assert_output "$BIN method=1 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx wmax=5 noscan=1 debug=0" 255 "noscan=1 only works with method=2" "ERROR in function"

# Test 8: Validation: noscan=1 with method=3
assert_output "$BIN method=3 finH=$EXAMPLES_DIR/surf_d5_H.mmx finL=$EXAMPLES_DIR/surf_d5_L.mmx wmax=5 noscan=1 debug=0" 255 "noscan=1 only works with method=2" "ERROR in function"

# Test 9: Validation: fdem and finH together
assert_output "$BIN method=2 fdem=$EXAMPLES_DIR/surf_d3.dem finH=$EXAMPLES_DIR/surf_d5_H.mmx wmax=3 debug=0" 255 "Cannot specify matrix files.*along with fdem" "ERROR in function"

# Test 10: Validation: pmin without fdem
assert_output "$BIN method=2 finH=$EXAMPLES_DIR/surf_d5_H.mmx wmax=3 pmin=0.01 debug=0" 255 "pmin can only be used when fdem is specified" "ERROR in function"

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

if [ $FAILED -ne 0 ]; then
    echo "Some tests failed!"
    exit 1
else
    echo "All tests passed!"
    exit 0
fi
