#!/usr/bin/env bash
# Golden tests: diff mytools against real bedtools on the fixtures in data/.
# Cases are the v1 list from SPEC.md section 8. stdout and exit status must match;
# stderr is never compared; nothing is sorted before diffing (row order matters).
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

# merge and closest require sorted input; a.bed and b.bed are deliberately not.
# Pre-sort with bedtools so both tools see identical, valid input.
for f in a b genes hg002.highconf; do
  bedtools sort -i "$DATA/$f.bed" > "$tmp/$f.sorted.bed"
done

# --- sort (#5) ---------------------------------------------------------------
check "sort a.bed"     -- sort -i "$DATA/a.bed"
check "sort b.bed"     -- sort -i "$DATA/b.bed"
check "sort genes.bed" -- sort -i "$DATA/genes.bed"
check "sort hg002"     -- sort -i "$DATA/hg002.highconf.bed"
check_stdin "sort stdin" "$DATA/a.bed" -- sort -i -

# --- merge (#6) --------------------------------------------------------------
for f in a b hg002.highconf; do
  check "merge $f.bed (sorted)" -- merge -i "$tmp/$f.sorted.bed"
done

# --- intersect (#7), subtract (#8), closest (#9) -----------------------------
# NB: with a.bed as -b, bedtools intersect/subtract exit 1 with no output
# ("illegal bin number -1"): a12 `chr2 0 0` is zero-length at position 0, which
# bedtools widens to [-1, 1) and cannot index. Matching the oracle there means
# exiting 1 too; compare() checks exit status, so it will hold us to it.
for pair in "a b" "b a" "genes hg002.highconf"; do
  set -- $pair
  check "intersect $1 $2" -- intersect -a "$DATA/$1.bed" -b "$DATA/$2.bed"
  check "subtract $1 $2"  -- subtract  -a "$DATA/$1.bed" -b "$DATA/$2.bed"
  check "closest $1 $2 (sorted)" -- closest -a "$tmp/$1.sorted.bed" -b "$tmp/$2.sorted.bed"
done
# intersect also runs with the real-data pair reversed (#7): hg002's -b records
# span several bin-tree levels, which exercises bedtools' within-record hit order.
check "intersect hg002.highconf genes" -- intersect -a "$DATA/hg002.highconf.bed" -b "$DATA/genes.bed"

echo "---"
echo "$pass passed, $fail failed"
[[ $fail -eq 0 ]]
