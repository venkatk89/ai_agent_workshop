# Shared setup for the unit tests: source R/bed.R (no bedtools needed).
source(testthat::test_path("..", "..", "R", "bed.R"))

MYTOOLS <- normalizePath(testthat::test_path("..", "..", "mytools"))

# Write `lines` to a temp .bed (or .bed.gz) and return its path.
write_fixture <- function(lines, gz = FALSE) {
  path <- tempfile(fileext = if (gz) ".bed.gz" else ".bed")
  con <- if (gz) gzfile(path, "wb") else file(path, "wb")
  on.exit(close(con))
  writeLines(lines, con, sep = "\n")
  path
}

# Read `lines` through read_bed_lines(), recording every on_record call.
# Returns list(records = <list of fields vectors>, error = <condition or NULL>,
# result = <reader return value or NULL>).
read_all <- function(lines, ...) {
  path <- write_fixture(lines)
  on.exit(unlink(path))
  con <- open_bed(path)
  on.exit(close(con), add = TRUE)
  seen <- list()
  err <- NULL
  res <- tryCatch(
    read_bed_lines(con, function(chrom, start, end, fields) {
      seen[[length(seen) + 1L]] <<- list(chrom = chrom, start = start, end = end, fields = fields)
    }, ...),
    mytools_error = function(e) { err <<- e; NULL })
  list(records = seen, error = err, result = res)
}

# Run ./mytools with args (and optional stdin text); capture status/stdout/stderr.
run_mytools <- function(args, input = NULL) {
  out <- tempfile(); err <- tempfile()
  on.exit(unlink(c(out, err)))
  status <- system2(MYTOOLS, args, stdout = out, stderr = err, input = input)
  list(status = status,
       stdout = readLines(out, warn = FALSE),
       stderr = readLines(err, warn = FALSE))
}
