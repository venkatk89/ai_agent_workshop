#!/usr/bin/env bash
# Unit tests: testthat, no bedtools required.
# Usage: ./tests/run_unit.sh
# Needs testthat: Rscript -e 'install.packages("testthat")'  (see tests/README.md)
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
Rscript -e '
  if (!requireNamespace("testthat", quietly = TRUE))
    stop("testthat is not installed: Rscript -e \"install.packages(\\\"testthat\\\")\"")
  res <- testthat::test_dir("'"$here"'/testthat", reporter = "summary", stop_on_failure = TRUE)
'
