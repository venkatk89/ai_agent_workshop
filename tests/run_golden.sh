#!/usr/bin/env bash
# Golden tests: diff mytools against real bedtools on the fixtures in data/.
# bedtools is the oracle -- if we differ, we are wrong.
#
# Usage: ./tests/run_golden.sh          (MYTOOLS=/path/to/build to test another binary)
set -uo pipefail

here=$(cd "$(dirname "$0")" && pwd)
MYTOOLS=${MYTOOLS:-$here/../mytools}
DATA=$here/../data
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
pass=0; fail=0

# compare <name> <got_rc> <want_rc>  -- shared reporting for both helpers
compare() {
  local name=$1 got_rc=$2 want_rc=$3
  if [[ $got_rc -ne $want_rc ]]; then
    echo "FAIL $name (exit $got_rc, bedtools gave $want_rc)"
    sed 's/^/      /' "$tmp/got.err" | head -3
    (( fail++ )); return
  fi
  if diff -q "$tmp/want" "$tmp/got" >/dev/null; then
    echo "ok   $name"; (( pass++ ))
  else
    echo "FAIL $name"
    diff -u "$tmp/want" "$tmp/got" | sed 's/^/      /' | head -20
    (( fail++ ))
  fi
}

# check <name> -- <args...>
#   runs "$MYTOOLS <args>" and "bedtools <args>", diffs stdout and exit codes
check() {
  local name=$1; shift; shift        # drop the literal --
  "$MYTOOLS" "$@" > "$tmp/got"  2>"$tmp/got.err"; local got_rc=$?
  bedtools   "$@" > "$tmp/want" 2>/dev/null;      local want_rc=$?
  compare "$name" "$got_rc" "$want_rc"
}

# check_stdin <name> <file> -- <args...>
#   same as check, but <file> is piped to both tools on stdin
check_stdin() {
  local name=$1 input=$2; shift; shift; shift
  "$MYTOOLS" "$@" < "$input" > "$tmp/got"  2>"$tmp/got.err"; local got_rc=$?
  bedtools   "$@" < "$input" > "$tmp/want" 2>/dev/null;      local want_rc=$?
  compare "$name" "$got_rc" "$want_rc"
}

# --- sort (#5) ---------------------------------------------------------------
check "sort a.bed"     -- sort -i "$DATA/a.bed"
check "sort b.bed"     -- sort -i "$DATA/b.bed"
check "sort genes.bed" -- sort -i "$DATA/genes.bed"
check "sort hg002"     -- sort -i "$DATA/hg002.highconf.bed"
check_stdin "sort stdin" "$DATA/a.bed" -- sort -i -

# merge and closest need sorted input -- sort into $tmp with bedtools first.

echo "---"
echo "$pass passed, $fail failed"
[[ $fail -eq 0 ]]
