# bed.R: shared BED reading, validation, overlap predicate and error reporting.
# Base R only. See SPEC.md §2, §3, §7 and issue #4 for the contract.
#
# Minimal version written for `sort` (#5); `check_sorted()` from #4 is not here yet.

# Print `mytools: <cmd>: [<file>:<line>: ]<msg>` to stderr and exit.
die <- function(cmd, msg, file = NULL, line = NULL, status = 1L) {
  where <- if (!is.null(file)) paste0(file, ":", line, ": ") else ""
  cat("mytools: ", cmd, ": ", where, msg, "\n", sep = "", file = stderr())
  quit(save = "no", status = status)
}

# Return an open read connection for a BED path.
# `-` / `stdin` -> standard input; `*.gz` -> gzfile(); otherwise file().
open_bed <- function(path, cmd = "mytools") {
  if (path %in% c("-", "stdin")) return(file("stdin", open = "r"))
  if (!file.exists(path) || file.access(path, 4L) != 0L) {
    die(cmd, paste("cannot open", path))
  }
  con <- if (grepl("\\.gz$", path)) gzfile(path, open = "r") else file(path, open = "r")
  con
}

is_header_line <- function(line) {
  startsWith(line, "#") || startsWith(line, "track") || startsWith(line, "browser")
}

# Stream a BED connection line by line, validating as we go.
#
# Calls on_record(chrom, start, end, fields) for every data line, where `fields`
# is the full character vector of tab-separated columns so callers can pass extra
# columns through untouched. `start`/`end` are numeric (bedtools accepts
# coordinates beyond 2^31, so we do not use integer).
#
# Header handling matches bedtools: the header is the *leading* run of lines
# starting with `#`, `track` or `browser`. The first line that is not one of those
# (a blank line included) ends it; header-looking lines anywhere later are skipped
# silently and never echoed, even with -header. Blank lines are skipped throughout.
# Only truly empty lines are blank -- a whitespace-only line is a data error in
# bedtools ("less than 3 columns"), so it is one here too.
#
# Validation is fail-fast on the first bad line, exit 1, messages per SPEC.md §7.
# Returns the header lines (character vector, possibly empty) invisibly.
read_bed_lines <- function(con, on_record, header = FALSE, cmd = "mytools",
                           file = "<stdin>", chunk = 1000L) {
  header_lines <- character(0)
  in_header <- TRUE
  ncol_expected <- NA_integer_
  lineno <- 0L
  repeat {
    lines <- readLines(con, n = chunk, warn = FALSE)
    if (length(lines) == 0L) break
    for (line in lines) {
      lineno <- lineno + 1L
      if (in_header) {
        if (is_header_line(line)) {
          if (header) header_lines <- c(header_lines, line)
          next
        }
        in_header <- FALSE
      }
      if (line == "" || is_header_line(line)) next

      fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      if (length(fields) < 3L) {
        die(cmd, "expected at least 3 tab-separated columns", file, lineno)
      }
      if (is.na(ncol_expected)) {
        ncol_expected <- length(fields)
      } else if (length(fields) != ncol_expected) {
        die(cmd, "differing number of BED fields", file, lineno)
      }
      if (!grepl("^[0-9]+$", fields[2L]) || !grepl("^[0-9]+$", fields[3L])) {
        die(cmd, "start and end must be non-negative integers", file, lineno)
      }
      start <- as.numeric(fields[2L])
      end <- as.numeric(fields[3L])
      if (start > end) die(cmd, "start is greater than end", file, lineno)

      on_record(fields[1L], start, end, fields)
    }
  }
  invisible(header_lines)
}

# The overlap predicate. BED is 0-based half-open, so two intervals overlap iff
# a.start < b.end AND b.start < a.end -- strict `<` on both sides. Bookended
# intervals (a_end == b_start) do NOT overlap. Every off-by-one in this project
# lives here; do not loosen it.
overlaps <- function(a_start, a_end, b_start, b_end) {
  a_start < b_end && b_start < a_end
}
