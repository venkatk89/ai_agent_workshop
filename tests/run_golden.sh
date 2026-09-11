#!/usr/bin/env bash
# Golden tests: diff mytools against real bedtools (the oracle) on data/.
# Cases are the v1 list from SPEC.md section 8. stdout and exit status must match;
# stderr is never compared. Row order is part of the answer -- nothing is sorted
# before diffing.
#
# Usage: ./tests/run_golden.sh          (MYTOOLS=path/to/mytools to test another build)
set -uo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
MYTOOLS=${MYTOOLS:-$HERE/../mytools}
DATA=$HERE/../data
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
pass=0; fail=0

command -v bedtools >/dev/null || { echo "bedtools not found; it is the oracle" >&2; exit 1; }

# check <name> -- <args...>
#   runs "$MYTOOLS <args>" and "bedtools <args>", diffs stdout and exit status
check() {
  local name=$1; shift; shift        # drop the literal --
  "$MYTOOLS" "$@" > "$tmp/got"  2>"$tmp/got.err"
  local got_rc=$?
  bedtools   "$@" > "$tmp/want" 2>/dev/null
  local want_rc=$?

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

# merge and closest require sorted input; a.bed and b.bed are deliberately not.
# Pre-sort with bedtools so both tools see identical, valid input.
for f in a b genes hg002.highconf; do
  bedtools sort -i "$DATA/$f.bed" > "$tmp/$f.sorted.bed"
done

for f in a b genes hg002.highconf; do
  check "sort $f.bed" -- sort -i "$DATA/$f.bed"
done

for f in a b hg002.highconf; do
  check "merge $f.bed (sorted)" -- merge -i "$tmp/$f.sorted.bed"
done

# NB: with a.bed as -b, bedtools intersect/subtract exit 1 with no output
# ("illegal bin number -1"): a12 `chr2 0 0` is zero-length at position 0, which
# bedtools widens to [-1, 1) and cannot index. Matching the oracle there means
# exiting 1 too; the check() below compares exit status, so it will hold us to it.
for pair in "a b" "b a" "genes hg002.highconf"; do
  set -- $pair
  check "intersect $1 $2" -- intersect -a "$DATA/$1.bed" -b "$DATA/$2.bed"
  check "subtract $1 $2"  -- subtract  -a "$DATA/$1.bed" -b "$DATA/$2.bed"
  check "closest $1 $2 (sorted)" -- closest -a "$tmp/$1.sorted.bed" -b "$tmp/$2.sorted.bed"
done

echo "---"
echo "$pass passed, $fail failed"
[[ $fail -eq 0 ]]
