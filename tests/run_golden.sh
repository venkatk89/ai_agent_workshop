#!/usr/bin/env bash
# Golden tests: diff mytools against real bedtools on the fixtures in data/.
#
# Usage: ./tests/run_golden.sh
#   MYTOOLS=/path/to/mytools ./tests/run_golden.sh   # test a different build
#
# bedtools is the oracle (SPEC.md §8). Each case runs the same argv through both
# tools and passes only if exit codes match AND stdout is byte-identical. stderr is
# never compared -- our wording is deliberately our own (SPEC.md §7).
set -uo pipefail

here=$(cd "$(dirname "$0")" && pwd)
MYTOOLS=${MYTOOLS:-$here/../mytools}
DATA=$here/../data
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
pass=0; fail=0

if ! command -v bedtools >/dev/null; then
  echo "run_golden.sh: bedtools not found on PATH; it is the oracle, nothing to diff against" >&2
  exit 1
fi

# run <stdin-file> <tool> <args...>
#   Runs "<tool> <args>" with stdin from <stdin-file>, capturing stdout and stderr in
#   $tmp. stdin is always redirected: bedtools blocks reading a terminal when a flag
#   is missing, and a blocked oracle looks like a hung test.
run() {
  local stdin=$1 tool=$2 tag=$3; shift 3
  "$tool" "$@" < "$stdin" > "$tmp/$tag" 2> "$tmp/$tag.err"
}

# check <name> -- <args...>
#   Diffs "$MYTOOLS <args>" against "bedtools <args>" with stdin closed.
check() { check_stdin /dev/null "$@"; }

# check_stdin <file> <name> -- <args...>
#   Same, but both tools read <file> on stdin (for "-i -" / "-a stdin" cases).
check_stdin() {
  local stdin=$1 name=$2; shift 3            # drop the literal --
  run "$stdin" "$MYTOOLS" got  "$@"; local got_rc=$?
  run "$stdin" bedtools   want "$@"; local want_rc=$?

  if [[ $got_rc -ne $want_rc ]]; then
    echo "FAIL $name (exit $got_rc, bedtools gave $want_rc)"
    sed 's/^/      /' "$tmp/got.err" | head -3
    (( fail++ )); return
  fi
  if cmp -s "$tmp/want" "$tmp/got"; then
    echo "ok   $name"; (( pass++ ))
  else
    echo "FAIL $name"
    diff -u "$tmp/want" "$tmp/got" | sed 's/^/      /' | head -20
    (( fail++ ))
  fi
}

# merge and closest require sorted input (SPEC.md §4); a.bed and b.bed are
# deliberately unsorted. Sort with bedtools so both tools read the identical file.
for f in a b genes hg002.highconf; do
  bedtools sort -i "$DATA/$f.bed" > "$tmp/$f.sorted.bed"
done
gzip -c "$DATA/a.bed" > "$tmp/a.bed.gz"
# Fixtures have no header lines; make one so -header has something to echo.
{ printf '# a comment header\ntrack name=hdr\nbrowser position chr1:1-100\n'; cat "$tmp/a.sorted.bed"; } > "$tmp/a.header.bed"

# ---- sort ----------------------------------------------------------------------
check "sort a.bed"              -- sort -i "$DATA/a.bed"
check "sort b.bed"              -- sort -i "$DATA/b.bed"
check "sort genes.bed"          -- sort -i "$DATA/genes.bed"
check "sort hg002"              -- sort -i "$DATA/hg002.highconf.bed"
check "sort a.bed.gz"           -- sort -i "$tmp/a.bed.gz"
check "sort -header"            -- sort -header -i "$tmp/a.header.bed"
check_stdin "$DATA/a.bed" "sort stdin (-)"     -- sort -i -
check_stdin "$DATA/a.bed" "sort stdin (stdin)" -- sort -i stdin

# ---- merge ---------------------------------------------------------------------
check "merge a.bed"             -- merge -i "$tmp/a.sorted.bed"
check "merge b.bed"             -- merge -i "$tmp/b.sorted.bed"
check "merge hg002"             -- merge -i "$tmp/hg002.highconf.sorted.bed"
check "merge -header"           -- merge -header -i "$tmp/a.header.bed"
check_stdin "$tmp/a.sorted.bed" "merge stdin" -- merge -i -
# Unsorted input: both exit 1 with empty stdout. The fixture itself is the unsorted case.
check "merge unsorted"          -- merge -i "$DATA/a.bed"

# ---- intersect -----------------------------------------------------------------
check "intersect a b"           -- intersect -a "$DATA/a.bed" -b "$DATA/b.bed"
# Oracle behaviour, not a typo: with a.bed as -b, bedtools intersect/subtract exit 1
# with empty stdout ("illegal bin number -1") because a12 is `chr2 0 0` and bedtools
# widens zero-length -b records to (start-1, end+1). closest does not use the bin
# tree and is unaffected. We match the exit code; see tests/README.md, "Gotchas".
check "intersect b a"           -- intersect -a "$DATA/b.bed" -b "$DATA/a.bed"
check "intersect genes hg002"   -- intersect -a "$DATA/genes.bed" -b "$DATA/hg002.highconf.bed"
check "intersect -header"       -- intersect -header -a "$tmp/a.header.bed" -b "$DATA/b.bed"
check_stdin "$DATA/a.bed" "intersect stdin -a" -- intersect -a - -b "$DATA/b.bed"
check_stdin "$DATA/b.bed" "intersect stdin -b" -- intersect -a "$DATA/a.bed" -b stdin

# ---- subtract ------------------------------------------------------------------
check "subtract a b"            -- subtract -a "$DATA/a.bed" -b "$DATA/b.bed"
check "subtract b a"            -- subtract -a "$DATA/b.bed" -b "$DATA/a.bed"   # exit 1, see intersect b a
check "subtract genes hg002"    -- subtract -a "$DATA/genes.bed" -b "$DATA/hg002.highconf.bed"
check "subtract -header"        -- subtract -header -a "$tmp/a.header.bed" -b "$DATA/b.bed"
check_stdin "$DATA/a.bed" "subtract stdin -a" -- subtract -a - -b "$DATA/b.bed"

# ---- closest -------------------------------------------------------------------
check "closest a b"             -- closest -a "$tmp/a.sorted.bed" -b "$tmp/b.sorted.bed"
check "closest b a"             -- closest -a "$tmp/b.sorted.bed" -b "$tmp/a.sorted.bed"
check "closest genes hg002"     -- closest -a "$tmp/genes.sorted.bed" -b "$tmp/hg002.highconf.sorted.bed"
check "closest -header"         -- closest -header -a "$tmp/a.header.bed" -b "$tmp/b.sorted.bed"
check_stdin "$tmp/a.sorted.bed" "closest stdin -a" -- closest -a - -b "$tmp/b.sorted.bed"
check "closest unsorted"        -- closest -a "$DATA/a.bed" -b "$DATA/b.bed"

# ---- error paths: exit code must match (stdout is empty for both) ----------------
check "sort missing file"       -- sort -i "$tmp/does-not-exist.bed"
check "sort unknown flag"       -- sort -i "$DATA/a.bed" -bogus
check "intersect missing -b"    -- intersect -a "$DATA/a.bed"
printf 'chr1\t10\t5\n'                       > "$tmp/start-gt-end.bed"
printf 'chr1 10 50\n'                        > "$tmp/space-separated.bed"
printf 'chr1\t10\t50\nchr1\t60\t70\tx\n'     > "$tmp/differing-columns.bed"
check "sort start>end"          -- sort -i "$tmp/start-gt-end.bed"
check "sort space-separated"    -- sort -i "$tmp/space-separated.bed"
check "intersect differing cols" -- intersect -a "$tmp/differing-columns.bed" -b "$DATA/b.bed"

echo "---"
echo "$pass passed, $fail failed"
[[ $fail -eq 0 ]]
