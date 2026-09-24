#!/bin/bash

COMPILERS=("gfortran" "ifx" "nvfortran")
# COMPILERS=("gfortran" "flang-new-18" "ifx" "lfortran")
FAILED=()

# Activate Intel oneAPI
source /opt/intel/oneapi/setvars.sh > /dev/null 2>&1

# Activate conda for lfortran
source ~/miniforge3/etc/profile.d/conda.sh > /dev/null 2>&1
conda activate base > /dev/null 2>&1

echo "=============================="
echo " FPM CI - Local"
echo "=============================="

for COMPILER in "${COMPILERS[@]}"; do
    echo ""
    echo "--- Building with $COMPILER ---"
    fpm build --compiler "$COMPILER" --profile debug
    if [ $? -ne 0 ]; then
        echo "❌ BUILD FAILED with $COMPILER"
        FAILED+=("$COMPILER (build)")
        continue
    fi

    echo "--- Testing with $COMPILER ---"
    fpm test --compiler "$COMPILER" --profile release
    if [ $? -ne 0 ]; then
        echo "❌ TESTS FAILED with $COMPILER"
        FAILED+=("$COMPILER (test)")
    else
        echo "✅ $COMPILER passed"
    fi
done

echo ""
echo "=============================="
if [ ${#FAILED[@]} -eq 0 ]; then
    echo "✅ All compilers passed!"
else
    echo "❌ Failures:"
    for F in "${FAILED[@]}"; do
        echo "   - $F"
    done
    exit 1
fi