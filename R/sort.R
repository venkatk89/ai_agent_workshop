# sort.R: `mytools sort -i <file> [-header]` -- byte-identical to `bedtools sort`.
#
# Order is chrom lexicographic (chr1 < chr10 < chr17 < chr2 < chr7 < chrX), then
# start ascending. Ties on (chrom, start) keep input order: bedtools' sort is stable
# and does NOT fall back to end (`chr1 10 50` stays before `chr1 10 20`).
#
# This is the one subcommand allowed to hold the whole input in memory (SPEC.md §6).

usage_sort <- function() {
  cat("usage: mytools sort -i <bed|-|stdin> [-header]\n", file = stderr())
}

parse_sort_args <- function(args) {
  opts <- list(input = NULL, header = FALSE)
  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (a == "-i") {
      if (i == length(args)) {
        usage_sort()
        die("sort", "-i requires a file argument")
      }
      opts$input <- args[i + 1L]
      i <- i + 2L
    } else if (a == "-header") {
      opts$header <- TRUE
      i <- i + 1L
    } else {
      usage_sort()
      die("sort", paste("unrecognized option", a))
    }
  }
  if (is.null(opts$input)) {
    usage_sort()
    die("sort", "-i is required")
  }
  opts
}

# Pure ordering helper, unit-tested on its own: returns the permutation that
# sorts records by chrom then start, stable on ties.
#
# method = "radix" matters twice over: it is stable, and it compares strings in
# the C locale (raw bytes), which is what bedtools does. The default method would
# use the session's collation and could put `Chr1` after `chr1` or reorder `-`
# against `_`. bedtools gives `Chr1 < chr1 < chr1-x < chr10 < chr1_random`.
sort_order <- function(chrom, start) {
  order(chrom, start, method = "radix")
}

main_sort <- function(args) {
  opts <- parse_sort_args(args)
  con <- open_bed(opts$input, cmd = "sort")
  on.exit(close(con), add = TRUE)

  # Grow-by-append is fine here: the largest promised input is a few hundred
  # lines, and sort is allowed to hold everything anyway.
  chrom <- character(0)
  start <- numeric(0)
  lines <- character(0)
  n <- 0L
  header_lines <- read_bed_lines(
    con,
    function(c, s, e, fields) {
      n <<- n + 1L
      chrom[n] <<- c
      start[n] <<- s
      lines[n] <<- paste(fields, collapse = "\t")
    },
    header = opts$header, cmd = "sort", file = opts$input
  )

  out <- stdout()
  if (length(header_lines)) writeLines(header_lines, out)
  if (n > 0L) writeLines(lines[sort_order(chrom, start)], out)
  0L
}
