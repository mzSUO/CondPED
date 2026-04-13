#!/bin/bash
echo "=== CondPED Package Health Check ==="

# 1. File structure
echo -e "\n1. File structure:"
ls -1 R/ tests/testthat/

# 2. Documentation
echo -e "\n2. Documentation generated:"
ls -1 man/

# 3. R CMD check
echo -e "\n3. Running R CMD check..."
R --quiet --no-save -e 'devtools::check(quiet = TRUE)' 2>&1 | tail -5

# 4. Test coverage
echo -e "\n4. Running tests..."
R --quiet --no-save -e 'devtools::test()' 2>&1 | tail -3

echo -e "\n=== Check complete ==="
