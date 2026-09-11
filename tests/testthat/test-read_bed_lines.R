# read_bed_lines(): streaming reader and its fail-fast validation.
# Every error case checks two things: the message, and that no record after the
# bad line reached on_record.

good <- c("chr1\t0\t100\ta01\t10\t+", "chr1\t100\t200\ta02\t20\t-")

test_that("well-formed BED6 is delivered record by record with numeric coordinates", {
  r <- read_all(good)
  expect_null(r$error)
  expect_length(r$records, 2)
  expect_equal(r$records[[1]]$chrom, "chr1")
  expect_equal(r$records[[1]]$start, 0)
  expect_equal(r$records[[1]]$end, 100)
  expect_type(r$records[[1]]$start, "double")
  expect_equal(r$records[[2]]$fields, c("chr1", "100", "200", "a02", "20", "-"))
  expect_equal(r$result$records, 2L)
  expect_equal(r$result$ncol, 6L)
})

test_that("fewer than 3 columns is an error and stops the stream", {
  r <- read_all(c(good[1], "chr1 300 400", good[2]))   # space-separated: 1 column
  expect_s3_class(r$error, "mytools_error")
  expect_match(conditionMessage(r$error), ":2: expected at least 3 tab-separated columns$")
  expect_equal(r$error$status, 1L)
  expect_length(r$records, 1)
})

test_that("a change in column count is an error and stops the stream", {
  r <- read_all(c(good, "chr1\t300\t400", "chr1\t500\t600\tx\t0\t+"))
  expect_match(conditionMessage(r$error), ":3: differing number of BED fields$")
  expect_length(r$records, 2)
})

test_that("non-integer or negative coordinates are an error and stop the stream", {
  for (bad in c("chr1\t-1\t10", "chr1\t1.5\t10", "chr1\t1e2\t200", "chr1\t+1\t10",
                "chr1\tabc\t10", "chr1\t\t10", "chr1\t 7\t10")) {
    r <- read_all(c(bad, "chr1\t0\t1"))
    expect_match(conditionMessage(r$error), ":1: start and end must be non-negative integers$",
                 info = bad)
    expect_length(r$records, 0)
  }
})

test_that("start greater than end is an error and stops the stream", {
  r <- read_all(c(good[1], "chr1\t200\t100\tbad\t0\t+", good[2]))
  expect_match(conditionMessage(r$error), ":2: start is greater than end$")
  expect_length(r$records, 1)
})

test_that("start equal to end (zero-length) is accepted", {
  r <- read_all(c("chr1\t500\t500", "chr2\t0\t0"))
  expect_null(r$error)
  expect_length(r$records, 2)
})

test_that("error messages carry cmd and file name", {
  r <- read_all("chr1\t5\t4", cmd = "merge", file = "in.bed")
  expect_equal(conditionMessage(r$error), "mytools: merge: in.bed:1: start is greater than end")
})

test_that("blank, #, track and browser lines are skipped and do not count as data", {
  lines <- c("#comment", "track name=x", "browser position chr1", "", good[1], "#mid", good[2])
  r <- read_all(lines)
  expect_null(r$error)
  expect_length(r$records, 2)
  expect_equal(r$result$header, character(0))
})

test_that("header = TRUE keeps the leading header run and streams it to on_header", {
  # Only the leading run is the header; "#late" after data is skipped silently,
  # as bedtools does (see test-bed.R for the bedtools sort -header evidence).
  streamed <- character(0)
  path <- write_fixture(c("#h1", "track x", good, "#late"))
  con <- open_bed(path)
  on.exit(close(con))
  res <- read_bed_lines(con, function(...) NULL, header = TRUE,
                        on_header = function(l) streamed <<- c(streamed, l))
  expect_equal(res$header, c("#h1", "track x"))
  expect_equal(streamed, res$header)
})

test_that("line numbers in errors count skipped lines, as bedtools does", {
  r <- read_all(c("#h", "", good[1], "chr1\t1\t2\t3"))
  expect_match(conditionMessage(r$error), ":4: differing number of BED fields$")
})

test_that("a trailing CR is dropped and one trailing tab adds no field (bedtools tokeniser)", {
  r <- read_all(c("chr1\t1\t2\r", "chr1\t1\t3\t"))
  expect_null(r$error)
  expect_equal(r$records[[1]]$fields, c("chr1", "1", "2"))
  expect_equal(r$records[[2]]$fields, c("chr1", "1", "3"))
  # ...but two trailing tabs do add an (empty) field, so the count differs.
  r <- read_all(c("chr1\t1\t2", "chr1\t1\t3\t\t"))
  expect_match(conditionMessage(r$error), "differing number of BED fields$")
})

test_that("whitespace-only lines are not blank lines", {
  # bedtools rejects "   " and "\t" as having < 3 columns; so do we.
  r <- read_all(c(good[1], "   "))
  expect_match(conditionMessage(r$error), ":2: expected at least 3 tab-separated columns$")
})

test_that("coordinates beyond 32-bit are read exactly", {
  r <- read_all("chr1\t3000000000\t3000000001")
  expect_null(r$error)
  expect_equal(r$records[[1]]$start, 3e9)
  expect_equal(format_coord(r$records[[1]]$start), "3000000000")
})

test_that("empty input yields no records and NA ncol", {
  r <- read_all(character(0))
  expect_null(r$error)
  expect_equal(r$result$records, 0L)
  expect_true(is.na(r$result$ncol))
})
