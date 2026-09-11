# intersect.R: `mytools intersect -a <file> -b <file> [-header]` -- byte-identical
# to `bedtools intersect` in its default mode.
#
# For each -a record, and for each -b record on the same chrom that it overlaps,
# print the overlapping portion (chrom, max start, min end) followed by the -a
# record's columns 4+ unchanged. -b is loaded fully and indexed per chrom; -a is
# streamed and written as we go (SPEC.md section 6).
#
# Two things here are dictated by matching the oracle rather than by the plain
# predicate in SPEC.md section 3. Both were verified against bedtools v2.31.1:
#
# 1. Zero-length records (start == end). bedtools widens them to [start-1, end+1)
#    before testing for overlap -- on BOTH sides -- so `chr1 500 500` hits
#    `chr1 400 500` and `chr1 500 600` (bookended, which would not otherwise
#    overlap) but not `chr1 400 499` or `chr1 501 600`. Two zero-length records
#    at 500 and 501 hit each other; 500 and 502 do not. When REPORTING the
#    overlap, only the -b side stays widened: the -b file is stored widened in
#    bedtools' tree, while the -a record is printed from its original coordinates.
#    Hence `-a chr1 400 600` vs `-b chr1 500 500` prints `chr1 499 501`, and
#    `-a chr1 500 500` vs anything it hits prints `chr1 500 500`.
#
#    A zero-length -b record at position 0 widens to [-1, 1), which bedtools
#    cannot place in its bin tree: it exits 1 with no output at all ("illegal bin
#    number -1"). `data/a.bed` has one (a12, `chr2 0 0`), so `intersect -a b.bed
#    -b a.bed` exits 1, and the golden test holds us to that.
#
# 2. Output order within one -a record. It is NOT -b file order. bedtools keeps
#    -b in a UCSC-style bin tree (finest bin 16 kb, x8 per level, 7 levels) and,
#    for each query, walks the levels from finest to coarsest, the bins in each
#    level in increasing order, and the records in each bin in insertion order.
#    A record lives in the finest bin that contains it whole. So a 100 bp -b
#    record is printed before a 50 kb one that appears earlier in the file. For
#    inputs where every record fits in the first 16 kb bin (a.bed, b.bed) this
#    degenerates to file order, but genes.bed / hg002.highconf.bed do not.
#    We reproduce the order by sorting each chrom's -b records once by
#    (level, bin, insertion index) and scanning them in that order.

BIN_FIRST_SHIFT <- 14L   # finest bin is 2^14 = 16 kb, as in bedtools' BinTree
BIN_NEXT_SHIFT <- 3L     # each level up is 8x wider
BIN_LEVELS <- 7L
# Offset of each level's first bin in bedtools' flat bin numbering, finest first.
BIN_OFFSETS <- c(32678 + 4096 + 512 + 64 + 8 + 1, 4096 + 512 + 64 + 8 + 1,
                 512 + 64 + 8 + 1, 64 + 8 + 1, 8 + 1, 1, 0)

usage_intersect <- function() {
  cat("usage: mytools intersect -a <bed|-|stdin> -b <bed|-|stdin> [-header]\n",
      file = stderr())
}

parse_intersect_args <- function(args) {
  opts <- list(a = NULL, b = NULL, header = FALSE)
  i <- 1L
  while (i <= length(args)) {
    arg <- args[i]
    if (arg %in% c("-a", "-b")) {
      if (i == length(args)) {
        usage_intersect()
        die("intersect", paste(arg, "requires a file argument"))
      }
      opts[[substring(arg, 2L)]] <- args[i + 1L]
      i <- i + 2L
    } else if (arg == "-header") {
      opts$header <- TRUE
      i <- i + 1L
    } else {
      usage_intersect()
      die("intersect", paste("unrecognized option", arg))
    }
  }
  for (flag in c("a", "b")) {
    if (is.null(opts[[flag]])) {
      usage_intersect()
      die("intersect", paste0("-", flag, " is required"))
    }
  }
  if (opts$a %in% c("-", "stdin") && opts$b %in% c("-", "stdin")) {
    usage_intersect()
    die("intersect", "only one of -a and -b may be stdin")
  }
  opts
}

# Widen a zero-length interval to [start-1, end+1), as bedtools does before
# testing overlap. Regular intervals are returned unchanged.
widen_zero_length <- function(start, end) {
  z <- start == end
  list(start = start - z, end = end + z)
}

# bedtools' BinTree::getBin(): the (level, bin) a [start, end) record is stored
# in -- the finest level at which start and end-1 fall in the same bin. Returns
# c(level, bin) with level 0 = finest, or NULL when bedtools would call the bin
# illegal (negative start, i.e. a widened zero-length record at position 0).
# Integer division by powers of two on doubles matches C++'s arithmetic shifts,
# including for -1 (floor(-1 / 16384) is -1, as -1 >> 14 is).
bin_of <- function(start, end) {
  s <- start %/% 2^BIN_FIRST_SHIFT
  e <- (end - 1) %/% 2^BIN_FIRST_SHIFT
  for (level in seq_len(BIN_LEVELS) - 1L) {
    if (s == e) {
      if (s < 0) return(NULL)
      return(c(level, s))
    }
    s <- s %/% 2^BIN_NEXT_SHIFT
    e <- e %/% 2^BIN_NEXT_SHIFT
  }
  NULL
}

# Load the whole -b file into a per-chrom index. Each chrom maps to a list of
# equal-length vectors (start, end -- widened where zero-length) already sorted
# into bedtools' hit order: bin level (finest first), then bin, then file order.
# A record bedtools cannot bin (see bin_of) is a hard error before any output.
load_b <- function(path) {
  con <- open_bed(path, cmd = "intersect")
  on.exit(close(con), add = TRUE)
  fname <- bed_file_name(con)

  chrom <- character(0); start <- numeric(0); end <- numeric(0)
  level <- integer(0); bin <- numeric(0)
  n <- 0L
  read_bed_lines(con, function(c, s, e, fields) {
    w <- widen_zero_length(s, e)
    lb <- bin_of(w$start, w$end)
    if (is.null(lb)) {
      die("intersect",
          paste0("cannot index -b record ", c, ":", format_coord(s), "-",
                 format_coord(e), " (bedtools: illegal bin number)"), fname)
    }
    n <<- n + 1L
    chrom[n] <<- c; start[n] <<- w$start; end[n] <<- w$end
    level[n] <<- lb[1L]; bin[n] <<- lb[2L]
  }, cmd = "intersect", file = fname)

  index <- list()
  for (ch in unique(chrom)) {
    i <- which(chrom == ch)
    i <- i[order(level[i], bin[i], i, method = "radix")]
    index[[ch]] <- list(start = start[i], end = end[i])
  }
  index
}

# The hits for one -a record against the -b index: a character vector of output
# lines (possibly empty), in bedtools' order. `fields` is the full -a record.
intersect_one <- function(chrom, start, end, fields, index) {
  b <- index[[chrom]]
  if (is.null(b)) return(character(0))
  q <- widen_zero_length(start, end)
  hit <- overlaps(q$start, q$end, b$start, b$end)
  if (!any(hit)) return(character(0))
  o_start <- pmax(start, b$start[hit])
  o_end <- pmin(end, b$end[hit])
  extra <- if (length(fields) > 3L) paste0("\t", paste(fields[-(1:3)], collapse = "\t")) else ""
  paste0(chrom, "\t", format_coord(o_start), "\t", format_coord(o_end), extra)
}

main_intersect <- function(args) {
  opts <- parse_intersect_args(args)
  index <- load_b(opts$b)

  con <- open_bed(opts$a, cmd = "intersect")
  on.exit(close(con), add = TRUE)
  out <- stdout()
  read_bed_lines(
    con,
    function(chrom, start, end, fields) {
      lines <- intersect_one(chrom, start, end, fields, index)
      if (length(lines)) writeLines(lines, out)
    },
    header = opts$header,
    on_header = function(line) writeLines(line, out),
    cmd = "intersect", file = bed_file_name(con)
  )
  0L
}
