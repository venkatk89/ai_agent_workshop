#!/usr/bin/env bash
# Lint every R file (the mytools entry point, R/, tests/). Exits non-zero on any lint.
# Config is .lintr at the repo root.
# Usage: ./tests/run_lint.sh
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
cd "$here/.."
Rscript - <<'RS'
files <- c("mytools", list.files(c("R", "tests"), pattern = "\\.R$", full.names = TRUE, recursive = TRUE))
lints <- do.call(c, lapply(files, function(f) lintr::lint(f, parse_settings = TRUE)))
if (length(lints) > 0) {
  print(lints)
  cat(sprintf("run_lint.sh: %d lint(s) found\n", length(lints)), file = stderr())
  quit(status = 1)
}
cat(sprintf("run_lint.sh: %d file(s) clean\n", length(files)))
RS
