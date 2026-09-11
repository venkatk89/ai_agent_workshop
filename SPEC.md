# SPEC.md — mytools

A small reimplementation of a subset of bedtools, written in R. Real `bedtools`
(v2.31.1, installed) is the oracle: where this document and bedtools disagree on the
same input, bedtools is right and this document is wrong.

Every decision below came out of the brainstorm interview. The italic line after each
is the rationale, so nobody re-litigates it by accident. The guiding answer, given
repeatedly, was **"match bedtools"** — when a question comes up that this spec does
not cover, run bedtools and encode what it does.

---

## 1. Scope

Subcommands in v1: `sort`, `merge`, `intersect`, `subtract`, `closest`.
*The full set from CLAUDE.md; all five are small enough to finish.*

Explicitly NOT in v1:

- Any flag not listed in §4. In particular no `-s`/`-S` (strand), `-f`/`-F`/`-r`
  (overlap fraction), `-v`, `-u`, `-wa`/`-wb`/`-wo`, `-c`/`-o` (merge aggregation),
  `-d` (merge distance / closest distance), `-t`/`-k` (closest ties), `-sorted`, `-g`.
- Formats other than BED: no GFF/GTF/VCF/BAM input.
- Strand awareness of any kind.
- Any subcommand beyond the five above.

*"Only what is strictly required" — scope creep is the failure mode here.*

## 2. Input formats

- **Format accepted:** BED3. Only columns 1–3 (`chrom`, `start`, `end`) are
  interpreted. Columns 4+ are never parsed; whether they appear in output is
  per-command and follows bedtools (see §5).
  *Keeps parsing trivial while still accepting real BED files.*
- **Column-count consistency:** every data line in a file must have the same number of
  columns. A change mid-file is a data error (exit 1). This is bedtools' behaviour
  (`Differing number of BED fields encountered at line: N`), not a choice.
- **Delimiter:** tab only. Space-separated lines are a data error, as in bedtools.
- **Input flags are required.** `sort` and `merge` take `-i <file>`; `intersect`,
  `subtract`, `closest` take `-a <file>` and `-b <file>`. Omitting one is a usage
  error. *bedtools-style; no silent defaulting to stdin.*
- **`-` and `stdin`** as a file value both mean standard input. *Match bedtools.*
- **Compressed input:** `.gz` files are accepted transparently for every input flag.
  *Real BED files ship gzipped; R's `gzfile()` makes this free.*
- **Non-data lines** — lines starting with `#`, `track`, or `browser`, and blank
  lines — are skipped for processing. By default they are **not** echoed to output.
  With `-header` they are printed before the results, exactly as bedtools does.
  *Match bedtools so golden tests compare byte-for-byte.*

## 3. Interval semantics

- Coordinates are **0-based, half-open**. `chr1 100 200` covers bases 100..199.
- Overlap predicate: `a.start < b.end && b.start < a.end`. Strict `<` on both sides.
- Bookended intervals (`a.end == b.start`) do **not** overlap. They **do** merge
  under `merge` (bedtools' default `-d 0`).
- Zero-length intervals (`start == end`) are legal input. Their behaviour under every
  subcommand is whatever bedtools does with them — encode the oracle, do not reason
  about it. `data/a.bed` and `data/b.bed` contain them deliberately.
- Minimum overlap to count: **1 bp, fixed**. No `-f`.
  *Keeps intersect/subtract to the single predicate above.*
- `start > end` is a **data error**, exit 1, matching bedtools.
- Negative or non-integer `start`/`end`: data error, exit 1, matching bedtools.

## 4. Flags per subcommand

Every flag name and meaning matches bedtools exactly.

| Subcommand  | Flags in v1        | Behaviour |
|-------------|--------------------|-----------|
| `sort`      | `-i`, `-header`    | Sort by chrom (**lexicographic**, bedtools' default: `chr1, chr10, chr2, …`), then start, then end. Stable. No `-sizeA`, `-g`, `-faidx`, etc. |
| `merge`     | `-i`, `-header`    | Merge overlapping and bookended intervals (`-d 0` semantics). Output is BED3 — extra columns are dropped, as bedtools does. |
| `intersect` | `-a`, `-b`, `-header` | Default bedtools mode: for each `-a` record, report the overlapping portion with each `-b` record it overlaps. Extra `-a` columns pass through. |
| `subtract`  | `-a`, `-b`, `-header` | Remove from each `-a` record the portions covered by any `-b` record; a record can split into several. Extra `-a` columns pass through. |
| `closest`   | `-a`, `-b`, `-header` | For each `-a` record, report the nearest `-b` record (overlap counts as distance 0). Output is the `-a` line followed by the full `-b` line. Ties: all tied `-b` records are reported (bedtools' default `-t all`). No `-b` record on that chrom: `-a` line followed by `.	-1	-1`. |

- Strand-aware flags: **none**.
- `-header` is the one non-input flag, on every subcommand, because header handling
  was decided as "match bedtools" and bedtools needs the flag to echo headers.
- **Sorted input:** `merge` and `closest` **require** input sorted by chrom then start
  and **error** (exit 1) on the first out-of-order record, matching bedtools.
  `intersect` and `subtract` accept input in any order, also matching bedtools.
  *Originally "sort internally"; changed after checking what bedtools actually does.*

## 5. Output

- **Byte-identical to bedtools** on the same input and flags. That is the whole
  output specification; the notes below only restate what that implies.
- Tab-separated fields, `\n` line endings, trailing newline on the last line.
- Empty result: zero bytes on stdout, exit 0.
- Extra columns: `sort`, `intersect`, `subtract` pass through `-i`/`-a` columns
  untouched; `merge` emits BED3 only; `closest` emits all `-a` columns then all `-b`
  columns.
- `-header`: header lines from `-i`/`-a` are printed first, in original order. Header
  lines from `-b` are never printed.
- stdout carries data only. Every message goes to stderr.

## 6. Memory model

Match bedtools' model, subcommand by subcommand:

| Subcommand  | Model |
|-------------|-------|
| `sort`      | Whole input in memory (the only command allowed to). |
| `merge`     | Streaming: one pass, holds only the current open interval. |
| `intersect` | `-b` fully in memory (indexed per chrom); `-a` streamed line by line. |
| `subtract`  | Same as `intersect`. |
| `closest`   | `-b` fully in memory (sorted per chrom); `-a` streamed line by line. |

- **Largest input promised:** `data/hg002.highconf.bed` (387 lines) as any input,
  on any subcommand. *The biggest real file in the workshop.*
- **Performance target:** finishes. No timing target in v1.
- Write output as it is produced; do not buffer the whole result.

## 7. Errors and exit codes

Exit codes **match bedtools**: `0` success, `1` for everything else — bad input data
*and* usage errors (bedtools has no separate usage code).

Message format on stderr: `mytools: <subcommand>: <reason>`, with the file name and
line number whenever a specific line is at fault. **Fail fast**: stop at the first bad
line; do not try to report all of them.
*Our own short format rather than bedtools' wording — bedtools' messages are
inconsistent across commands and stderr is not part of the golden comparison.*

| Situation | stderr message | exit |
|-----------|----------------|------|
| Success | — | `0` |
| Fewer than 3 columns | `mytools: <cmd>: <file>:<line>: expected at least 3 tab-separated columns` | `1` |
| Column count differs from earlier lines | `mytools: <cmd>: <file>:<line>: differing number of BED fields` | `1` |
| Non-integer or negative `start`/`end` | `mytools: <cmd>: <file>:<line>: start and end must be non-negative integers` | `1` |
| `start > end` | `mytools: <cmd>: <file>:<line>: start is greater than end` | `1` |
| Unsorted input to `merge`/`closest` | `mytools: <cmd>: <file>:<line>: input is not sorted by chrom then start` | `1` |
| Missing / unreadable input file | `mytools: <cmd>: cannot open <file>` | `1` |
| Unknown flag | `mytools: <cmd>: unrecognized option <flag>` + one-line usage | `1` |
| Required input flag missing | `mytools: <cmd>: -i is required` (or `-a`/`-b`) + one-line usage | `1` |
| Unknown subcommand / no subcommand | `mytools: unknown subcommand <x>` / usage | `1` |
| `--version` | prints `mytools <version>` on **stdout** | `0` |

## 8. Correctness

- **Oracle:** real `bedtools` on the files in `data/`. Non-negotiable.
- **Golden tests in v1** (`tests/run_golden.sh`, diffs our stdout against bedtools'):
  - `sort -i` on `a.bed`, `b.bed`, `genes.bed`, `hg002.highconf.bed`
  - `merge -i` on sorted `a.bed`, sorted `b.bed`, sorted `hg002.highconf.bed`
  - `intersect -a -b`, `subtract -a -b`, `closest -a -b` on (`a.bed`, `b.bed`),
    (`b.bed`, `a.bed`), and (`genes.bed`, `hg002.highconf.bed`) — sorted copies
    where the command requires it.
  - Golden inputs to `merge`/`closest` are pre-sorted with `bedtools sort` so both
    tools see identical, valid input.
  - Nothing beyond the above for now (no `.gz`, `-header`, or error-path goldens).
- **Unit tests** (`tests/`, `testthat`, run without bedtools), one per edge case: the
  overlap predicate itself, bookended, zero-length, nested, position 0, and the
  unsorted-input check.
- **Known deviations from bedtools**, all deliberate:
  1. stderr wording is our own — stderr is never golden-compared.
  2. Nothing else. If a golden test needs a deviation to pass, the deviation is a bug.

## 9. Language and layout

- **Language:** R (≥ 4.5). **Runtime:** standard library only, no CRAN dependencies.
  **Test time:** the ban does not apply; unit tests use `testthat`, installed from
  CRAN (`Rscript -e 'install.packages("testthat")'`, documented in `tests/README.md`).
- **Entry point:** `./mytools` at the repo root — an executable `Rscript` file that
  parses the subcommand and `source()`s the matching implementation.
- **One file per subcommand:** `R/sort.R`, `R/merge.R`, `R/intersect.R`,
  `R/subtract.R`, `R/closest.R`. Shared BED reading/validation and the overlap
  predicate live in `R/bed.R`. *Agents work on subcommands in parallel; one file each
  keeps their diffs from colliding.*
- **Tests:** `tests/run_golden.sh` (oracle diff) and `testthat` unit tests under
  `tests/`.
- **Version:** a single `VERSION` constant in `mytools`, printed by `--version`.
  Stays `0.1.0` while subcommands are being built; bump to `1.0.0` when all five pass
  their golden tests.

## 10. Follow-ups this spec creates

Decisions above that conflict with things already in the repo. Fix them when the first
subcommand lands, not before.

1. **CLAUDE.md** says usage errors exit `2`. §7 says `1` (match bedtools). Update the
   `Code` section of CLAUDE.md.
2. **`mytools` with no arguments** currently exits `2` (issue #1's acceptance
   criterion). Under §7 it exits `1`. Change it in the same commit as (1).
3. **`testthat`** is not installed on this machine yet.
