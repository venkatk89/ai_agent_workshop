# Source the implementation under test. Tests run from any working directory:
# testthat sets the wd to tests/testthat, so R/ is two levels up.
root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
source(file.path(root, "R", "bed.R"))
source(file.path(root, "R", "sort.R"))
source(file.path(root, "R", "intersect.R"))

# Read a BED string through read_bed_lines(), returning the records and headers.
read_bed_string <- function(text, header = FALSE) {
  con <- textConnection(text)
  on.exit(close(con))
  recs <- list()
  hdr <- read_bed_lines(con, function(chrom, start, end, fields) {
    recs[[length(recs) + 1L]] <<- list(chrom = chrom, start = start, end = end,
                                       fields = fields)
  }, header = header, cmd = "test", file = "<string>")
  list(records = recs, header = hdr$header)
}
