# R/merge.R: mytools merge -i <file> [-header]
#
# Merges overlapping and bookended intervals (bedtools' default `-d 0`) and prints
# BED3. Streaming: one pass over the input, holding only the current open block.
# Input must be sorted (each chromosome contiguous, start non-decreasing within it;
# SPEC.md section 4) -- the reader dies on the first out-of-order record.

usage_merge <- function() {
  cat("usage: mytools merge -i <bed> [-header]\n", file = stderr())
}

main_merge <- function(args) {
  cmd <- "merge"
  input <- NULL
  header <- FALSE
  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (a == "-i") {
      if (i == length(args)) {
        usage_merge()
        die(cmd, "-i is required")
      }
      # Repeated -i: the first one wins, as in bedtools.
      if (is.null(input)) input <- args[i + 1L]
      i <- i + 2L
    } else if (a == "-header") {
      header <- TRUE
      i <- i + 1L
    } else {
      usage_merge()
      die(cmd, paste0("unrecognized option ", a))
    }
  }
  if (is.null(input)) {
    usage_merge()
    die(cmd, "-i is required")
  }

  con <- open_bed(input, cmd)
  on.exit(close(con))
  file <- if (input %in% c("-", "stdin")) "stdin" else input
  merge_stream(con, header = header, cmd = cmd, file = file)
  invisible(0L)
}

# merge_stream(): read `con` and write merged BED3 blocks to stdout as they close.
#
# The open block is (chrom, start, end, n, orig_start, orig_end) where n is the
# number of records folded into it. A record joins the block when it is on the
# same chromosome and its start is <= the block's end (overlap or bookended);
# the block's end becomes the larger of the two.
#
# Known deviation (documented, not encoded): bedtools ends the open block at a
# blank line mid-file (`chr1 10 20`, blank, `chr1 15 30` prints two blocks). Our
# reader skips blank lines before merge sees them (SPEC.md section 2). No fixture
# has a mid-file blank line. Also, on an unsorted record bedtools' stdout holds
# whatever its output buffer had flushed; we stream, so partial output differs.
#
# Zero-length records -- oracle behaviour, not reasoning (bedtools v2.31.1):
# bedtools widens a record with start == end to [start-1, end+1) before it takes
# part in merging, and the widened coordinates are what get printed when the
# block contains more than one record. A block of exactly one record prints the
# record as read. Verified cases:
#   chr1 500 500 + chr1 500 600   -> chr1 499 600   (sorted a.bed; the README case)
#   chr1 10 20   + chr1 21 21     -> chr1 10 22     (widened [20,22) is bookended)
#   chr1 10 20   + chr1 22 22     -> two blocks, second printed as chr1 22 22
#   chr1 10 10   + chr1 10 10     -> chr1 9 11
#   chr1 0 0     + chr1 0 10      -> chr1 -1 10     (yes, negative)
#   chr1 500 500 alone            -> chr1 500 500
merge_stream <- function(con, header = FALSE, cmd = bed_cmd(), file = bed_file_name(con)) {
  # Open-block state, in an environment so the reader callbacks can update it.
  blk <- new.env()
  blk$chrom <- NULL
  blk$n <- 0L
  blk$seen_data <- FALSE

  flush <- function() {
    if (is.null(blk$chrom)) return(invisible())
    if (blk$n == 1L) {
      write_bed_line(c(blk$chrom, format_coord(blk$orig_start), format_coord(blk$orig_end)))
    } else {
      write_bed_line(c(blk$chrom, format_coord(blk$start), format_coord(blk$end)))
    }
  }

  on_record <- function(chrom, start, end, fields) {
    blk$seen_data <- TRUE
    ws <- start
    we <- end
    if (start == end) {
      ws <- start - 1
      we <- end + 1
    }
    if (!is.null(blk$chrom) && chrom == blk$chrom && ws <= blk$end) {
      blk$end <- max(blk$end, we)
      blk$n <- blk$n + 1L
      return(invisible())
    }
    flush()
    blk$chrom <- chrom
    blk$start <- ws
    blk$end <- we
    blk$orig_start <- start
    blk$orig_end <- end
    blk$n <- 1L
  }

  # bedtools -header prints only the header lines that precede the first data
  # line; a comment mid-file is skipped like any other.
  on_header <- function(line) {
    if (!blk$seen_data) cat(line, "\n", sep = "")
  }

  read_bed_lines(con, on_record, header = header, on_header = on_header,
                 sorted = TRUE, chrom_order = "grouped", cmd = cmd, file = file)
  flush()
  invisible()
}
