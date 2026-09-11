#!/usr/bin/env bash
# Unit tests (testthat). Need R and the testthat package, not bedtools.
#   Rscript -e 'install.packages("testthat")'
# Usage: ./tests/run_unit.sh
set -euo pipefail
cd "$(dirname "$0")/.."
Rscript -e 'library(testthat); res <- test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)'
