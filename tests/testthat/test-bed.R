# Unit tests for the shared BED core in R/bed.R. Runs without bedtools.
# Validation-error tests (which need to intercept die()) belong to #4.

# --- overlap predicate --------------------------------------------------------

test_that("overlapping intervals overlap", {
  expect_true(overlaps(100, 200, 150, 250))
  expect_true(overlaps(150, 250, 100, 200))
})

test_that("disjoint intervals do not overlap", {
  expect_false(overlaps(100, 200, 300, 400))
  expect_false(overlaps(300, 400, 100, 200))
})

test_that("bookended intervals do not overlap (half-open: a.end == b.start)", {
  expect_false(overlaps(0, 100, 100, 200))
  expect_false(overlaps(100, 200, 0, 100))
})

test_that("a nested interval overlaps its container", {
  expect_true(overlaps(300, 400, 320, 350))
  expect_true(overlaps(320, 350, 300, 400))
})

test_that("a zero-length interval strictly inside an interval overlaps it", {
  expect_true(overlaps(500, 500, 400, 600))
  expect_true(overlaps(400, 600, 500, 500))
})

test_that("a zero-length interval at another's boundary does not overlap it", {
  # chr1 500 500 vs chr1 500 600: a.start < b.end (500 < 600) but b.start < a.end
  # is 500 < 500, false. Two zero-length intervals never overlap each other.
  expect_false(overlaps(500, 500, 500, 600))
  expect_false(overlaps(500, 500, 500, 500))
})

test_that("an interval at position 0 overlaps as usual", {
  expect_true(overlaps(0, 100, 50, 150))
  expect_true(overlaps(0, 100, 0, 1))
  expect_false(overlaps(0, 100, 100, 200))
})

# --- reader: header and blank-line handling (matches bedtools) -----------------

test_that("data lines yield chrom, numeric start/end, and all fields", {
  r <- read_bed_string("chr1\t0\t100\ta01\t10\t+")
  expect_length(r$records, 1L)
  expect_equal(r$records[[1]]$chrom, "chr1")
  expect_equal(r$records[[1]]$start, 0)
  expect_equal(r$records[[1]]$end, 100)
  expect_equal(r$records[[1]]$fields, c("chr1", "0", "100", "a01", "10", "+"))
})

test_that("leading #/track/browser lines form the header and are returned only with header=TRUE", {
  text <- c("#a", "track t", "browser b", "chr1\t0\t10")
  expect_equal(read_bed_string(text, header = TRUE)$header, c("#a", "track t", "browser b"))
  expect_equal(read_bed_string(text, header = FALSE)$header, character(0))
  expect_length(read_bed_string(text)$records, 1L)
})

test_that("header-looking lines after the first data line are skipped, not echoed", {
  # bedtools sort -header on "chr2 5 9 / #mid / chr1 10 50" prints no header at all.
  r <- read_bed_string(c("chr2\t5\t9", "#mid", "chr1\t10\t50", "browser x"), header = TRUE)
  expect_equal(r$header, character(0))
  expect_length(r$records, 2L)
})

test_that("a blank line ends the header; a later # line is dropped", {
  # bedtools sort -header on "#a / <blank> / #b / chr2 5 9" prints only "#a".
  r <- read_bed_string(c("#a", "", "#b", "chr2\t5\t9"), header = TRUE)
  expect_equal(r$header, "#a")
  expect_length(r$records, 1L)
})

test_that("a leading blank line means there is no header at all", {
  # bedtools sort -header on "<blank> / #top / chr2 5 9" prints no header.
  r <- read_bed_string(c("", "#top", "chr2\t5\t9"), header = TRUE)
  expect_equal(r$header, character(0))
  expect_length(r$records, 1L)
})

test_that("blank lines between data lines are skipped", {
  r <- read_bed_string(c("chr1\t0\t10", "", "chr1\t20\t30"))
  expect_length(r$records, 2L)
})

test_that("an empty input yields no records and no header", {
  r <- read_bed_string(character(0), header = TRUE)
  expect_length(r$records, 0L)
  expect_equal(r$header, character(0))
})
