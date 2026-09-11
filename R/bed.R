# R/bed.R -- shared BED foundation: reading, validation, the overlap predicate,
# sorted-order checking and error reporting. Every subcommand builds on this.
# Base R only. See SPEC.md sections 2, 3 and 7; behaviour is bedtools' (v2.31.1).

# Process-wide settings every subcommand depends on. Set here, once, so unit tests
# that source() this file get the same behaviour as the ./mytools dispatcher.
#  - bedtools compares chromosome names bytewise (std::string), so "chrM" < "chra".
#    In a UTF-8 locale R's `<`, sort() and order() would collate differently.
#  - Coordinates are held as doubles (BED allows > 2^31; bedtools accepts
#    3000000000). Without scipen, paste(3e9) prints "3e+09" and corrupts output.
invisible(Sys.setlocale("LC_COLLATE", "C"))
options(scipen = 999)

SUBCOMMANDS <- c("sort", "merge", "intersect", "subtract", "closest")

# Which subcommand is running; the dispatcher sets it so error messages read
# `mytools: <cmd>: ...` without every helper needing cmd threaded through.
bed_cmd <- function() getOption("mytools.cmd", "mytools")

# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------

# die(): report `mytools: <cmd>: [<file>:<line>: ]<msg>` on stderr and exit with
# `status` (1 by default -- bedtools has no separate usage code, SPEC.md section 7).
#
# Implemented as a classed condition ("mytools_error") rather than a direct quit()
# so that unit tests can catch it with expect_error(class = "mytools_error"). The
# dispatcher in ./mytools is the one place that turns the condition into the stderr
# line and the exit status. Called outside the dispatcher it is an ordinary R error.
die <- function(cmd, msg, file = NULL, line = NULL, status = 1L) {
  loc <- ""
  if (!is.null(file)) {
    loc <- if (!is.null(line)) paste0(file, ":", line, ": ") else paste0(file, ": ")
  }
  full <- paste0("mytools: ", cmd, ": ", loc, msg)
  stop(structure(
    class = c("mytools_error", "error", "condition"),
    list(message = full, call = NULL, status = as.integer(status))
  ))
}

# ---------------------------------------------------------------------------
# Interval semantics
# ---------------------------------------------------------------------------

# BED is 0-based, half-open: [start, end). Two intervals overlap iff each starts
# before the other ends -- strict `<` on both sides, so bookended intervals
# (a_end == b_start) do NOT overlap. Vectorised; every off-by-one lives here.
#
# Zero-length records (start == end) are a subcommand concern, not this
# predicate's: bedtools widens them to [start-1, end+1) before testing, so
# `chr1 500 500` hits `chr1 400 500` and `chr1 500 600` but not `chr1 400 499`
# or `chr1 501 600` (verified with `bedtools intersect`; `merge` on sorted a.bed
# reports `chr1 499 600` for the same reason). Do that adjustment in the caller
# and match the oracle; keep this function the plain predicate.
overlaps <- function(a_start, a_end, b_start, b_end) {
  a_start < b_end & b_start < a_end
}

# Format a coordinate for output. Coordinates are doubles; sprintf("%.0f") is
# exact for integers below 2^53 and never falls back to scientific notation.
format_coord <- function(x) sprintf("%.0f", x)

# Write one BED record to stdout: fields joined by tab, "\n" terminated.
# Pass character fields (as read) or already-formatted coordinates.
write_bed_line <- function(fields) {
  cat(paste(fields, collapse = "\t"), "\n", sep = "")
}

# ---------------------------------------------------------------------------
# Opening input
# ---------------------------------------------------------------------------

# open_bed(): return a text-mode connection, opened for reading.
#   "-" or "stdin"  -> the process's standard input
#   "*.gz"          -> gzfile()
#   anything else   -> file()
# Missing or unreadable path -> `cannot open <file>`, exit 1 (via die()).
# The caller owns the connection and should close() it when done.
open_bed <- function(path, cmd = bed_cmd()) {
  if (identical(path, "-") || identical(path, "stdin")) {
    return(file("stdin", open = "rt"))
  }
  if (!file.exists(path) || dir.exists(path) || file.access(path, mode = 4L) != 0L) {
    die(cmd, paste0("cannot open ", path))
  }
  con <- if (grepl("\\.gz$", path)) gzfile(path) else file(path)
  ok <- tryCatch({ open(con, open = "rt"); TRUE },
                 error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok) {
    try(close(con), silent = TRUE)
    die(cmd, paste0("cannot open ", path))
  }
  con
}

# Name used in error messages for a connection: the path for file()/gzfile(),
# "stdin" for standard input.
bed_file_name <- function(con) {
  desc <- summary(con)$description
  if (is.null(desc) || !nzchar(desc)) "stdin" else desc
}

# ---------------------------------------------------------------------------
# Streaming reader
# ---------------------------------------------------------------------------

# Is this a non-data line? bedtools skips lines starting with "#", "track" or
# "browser" (prefix match, no separator required) and empty lines. A line of only
# spaces or tabs is NOT blank -- bedtools rejects it as having < 3 columns, and
# so do we (it reaches the column check below).
is_header_line <- function(line) {
  startsWith(line, "#") || startsWith(line, "track") || startsWith(line, "browser")
}

# read_bed_lines(): stream a BED connection record by record.
#
#   con        connection from open_bed()
#   on_record  function(chrom, start, end, fields) called for each data line, in
#              file order. start/end are numeric (doubles, parsed); fields is the
#              full character vector of tab-separated columns exactly as read, so
#              extra columns pass through untouched. NB: bedtools re-emits start/
#              end from the parsed value ("007" comes out as "7"); use fields for
#              columns 4+ and format_coord(start/end) for columns 2-3 if you want
#              to match that.
#   header     FALSE (default): header lines are skipped silently.
#              TRUE: header lines are kept, in order, in the returned $header and
#              handed to on_header as they are met -- they always precede the
#              data, so a streaming caller can print them before its first record.
#   on_header  optional function(line), used only when header = TRUE.
#
# What counts as the header matches bedtools: the *leading* run of lines starting
# with "#", "track" or "browser". The first line that is not one of those -- a
# blank line included -- ends it. Header-looking lines anywhere later are skipped
# silently and never echoed, even with -header (bedtools sort -header on
# "chr2 5 9 / #mid / chr1 10 50" prints no header; on "#a / <blank> / #b / data"
# it prints only "#a").
#   sorted     TRUE makes the reader die with `input is not sorted by chrom then
#              start` at the first out-of-order record (merge and closest).
#   chrom_order  passed to check_sorted() when sorted = TRUE: "lexicographic"
#              (closest) or "grouped" (merge). See check_sorted() for why.
#   cmd, file  for error messages; default to the running subcommand and the
#              connection's path.
#
# Validation is fail-fast: the first bad line raises (via die()) and on_record is
# never called for it or anything after it. Messages match SPEC.md section 7:
#   < 3 columns                      expected at least 3 tab-separated columns
#   column count != first data line  differing number of BED fields
#   start/end not /^[0-9]+$/         start and end must be non-negative integers
#   start > end                      start is greater than end
# Line numbers in messages are 1-based physical line numbers (headers/blanks
# count), as in bedtools.
#
# Tokenising follows bedtools: a trailing "\r" (CRLF files) is dropped, then the
# line is split on tab only. One trailing tab yields no extra field ("a\tb\t" is
# two fields, as bedtools does); two yield one empty field. R's strsplit() does
# exactly this, which is why it is used unadorned.
#
# Returns (invisibly) list(records = <n data lines>, ncol = <column count, or NA
# if there were none>, header = <character vector of header lines kept>). ncol
# is what `closest` needs to size its placeholder for a -b file.
read_bed_lines <- function(con, on_record, header = FALSE, on_header = NULL,
                           sorted = FALSE, chrom_order = "lexicographic",
                           cmd = bed_cmd(), file = bed_file_name(con)) {
  lineno <- 0L
  nrec <- 0L
  ncol <- NA_integer_
  headers <- character(0)
  in_header <- TRUE
  check <- if (isTRUE(sorted)) check_sorted(cmd, file, chrom_order) else NULL

  repeat {
    # Read in modest chunks for speed; memory stays bounded and no line is
    # processed past the first error.
    chunk <- readLines(con, n = 512L, warn = FALSE)
    if (length(chunk) == 0L) break
    for (line in chunk) {
      lineno <- lineno + 1L
      if (endsWith(line, "\r")) line <- substr(line, 1L, nchar(line) - 1L)
      if (in_header) {
        if (is_header_line(line)) {
          if (isTRUE(header)) {
            headers <- c(headers, line)
            if (is.function(on_header)) on_header(line)
          }
          next
        }
        in_header <- FALSE
      }
      if (!nzchar(line) || is_header_line(line)) next

      fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      if (length(fields) < 3L) {
        die(cmd, "expected at least 3 tab-separated columns", file, lineno)
      }
      if (is.na(ncol)) {
        ncol <- length(fields)
      } else if (length(fields) != ncol) {
        die(cmd, "differing number of BED fields", file, lineno)
      }
      if (!grepl("^[0-9]+$", fields[2L]) || !grepl("^[0-9]+$", fields[3L])) {
        die(cmd, "start and end must be non-negative integers", file, lineno)
      }
      start <- as.numeric(fields[2L])
      end <- as.numeric(fields[3L])
      if (start > end) {
        die(cmd, "start is greater than end", file, lineno)
      }
      if (!is.null(check)) check(fields[1L], start, lineno)

      nrec <- nrec + 1L
      on_record(fields[1L], start, end, fields)
    }
  }
  invisible(list(records = nrec, ncol = ncol, header = headers))
}

# ---------------------------------------------------------------------------
# Sorted-input check (merge, closest)
# ---------------------------------------------------------------------------

# check_sorted(): returns a stateful checker `function(chrom, start, line)` that
# dies with `input is not sorted by chrom then start` on the first record that is
# out of order. Within a chromosome, start must not decrease; equal starts are in
# order and end is never consulted (bedtools accepts `chr1 1 9` then `chr1 1 5`).
#
# chrom_order says what counts as out of order between chromosomes -- bedtools
# v2.31.1 applies two different rules, verified by running it:
#   "lexicographic"  chrom must be non-decreasing bytewise (chr1 < chr10 < chr2 <
#                    chrX). `closest` rejects `chr2` then `chr1` outright.
#   "grouped"        any chromosome order is fine, but once a chromosome ends it
#                    must not reappear. `merge` happily processes `chr2` then
#                    `chr1`, yet rejects `chr1, chr2, chr1`.
# LC_COLLATE=C (set above) is what makes `<` on chrom names bytewise.
check_sorted <- function(cmd, file, chrom_order = c("lexicographic", "grouped")) {
  chrom_order <- match.arg(chrom_order)
  prev_chrom <- NULL
  prev_start <- NULL
  seen <- character(0)
  function(chrom, start, line = NULL) {
    if (!is.null(prev_chrom)) {
      out_of_order <- if (chrom == prev_chrom) {
        start < prev_start
      } else if (chrom_order == "lexicographic") {
        chrom < prev_chrom
      } else {
        chrom %in% seen
      }
      if (out_of_order) {
        die(cmd, "input is not sorted by chrom then start", file, line)
      }
      if (chrom != prev_chrom) seen <<- c(seen, prev_chrom)
    }
    prev_chrom <<- chrom
    prev_start <<- start
    invisible(TRUE)
  }
}
