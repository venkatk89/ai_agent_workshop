# mytools merge: one test per behaviour, no bedtools needed. Expected values that
# look odd are bedtools v2.31.1's output on the same lines (see R/merge.R).

source(testthat::test_path("..", "..", "R", "merge.R"))

# Run merge_stream() over `lines` and return stdout as a character vector.
merge_lines <- function(lines, ...) {
  path <- write_fixture(lines)
  on.exit(unlink(path))
  con <- open_bed(path)
  on.exit(close(con), add = TRUE)
  capture.output(merge_stream(con, cmd = "merge", file = "in.bed", ...))
}

test_that("bookended intervals merge into one block", {
  expect_equal(merge_lines(c("chr1\t100\t200", "chr1\t200\t300")), "chr1\t100\t300")
})

test_that("intervals separated by one base do not merge", {
  expect_equal(merge_lines(c("chr1\t100\t200", "chr1\t201\t300")),
               c("chr1\t100\t200", "chr1\t201\t300"))
})

test_that("a nested interval collapses into the outer one, keeping the outer end", {
  expect_equal(merge_lines(c("chr1\t100\t200", "chr1\t100\t150", "chr1\t120\t130")),
               "chr1\t100\t200")
})

test_that("a chromosome change flushes the open block", {
  expect_equal(merge_lines(c("chr1\t100\t200", "chr2\t100\t200")),
               c("chr1\t100\t200", "chr2\t100\t200"))
})

test_that("extra columns are dropped: output is BED3", {
  expect_equal(merge_lines(c("chr1\t100\t200\ta\t0\t+", "chr1\t150\t250\tb\t0\t-")),
               "chr1\t100\t250")
})

test_that("a lone zero-length record is printed as read", {
  expect_equal(merge_lines("chr1\t500\t500"), "chr1\t500\t500")
})

test_that("a zero-length record that merges is widened by one base each side (oracle)", {
  # bedtools merge on sorted a.bed prints chr1 499 600 for a07 + a08.
  expect_equal(merge_lines(c("chr1\t500\t500", "chr1\t500\t600")), "chr1\t499\t600")
  # Widened [20,22) is bookended with [10,20), so they merge and the end is 22.
  expect_equal(merge_lines(c("chr1\t10\t20", "chr1\t21\t21")), "chr1\t10\t22")
  # Two zero-length records at 0 print a negative start; bedtools does too.
  expect_equal(merge_lines(c("chr1\t0\t0", "chr1\t0\t0")), "chr1\t-1\t1")
})

test_that("a zero-length record two bases past the block starts its own block", {
  expect_equal(merge_lines(c("chr1\t10\t20", "chr1\t22\t22")),
               c("chr1\t10\t20", "chr1\t22\t22"))
})

test_that("unsorted input raises the error at the offending line number", {
  err <- expect_error(merge_lines(c("#h", "chr1\t100\t200", "chr1\t50\t80")),
                      class = "mytools_error")
  expect_equal(conditionMessage(err),
               "mytools: merge: in.bed:3: input is not sorted by chrom then start")
})

test_that("chromosomes may come in any order but may not reappear (bedtools merge)", {
  expect_equal(merge_lines(c("chr2\t0\t10", "chr1\t0\t10")), c("chr2\t0\t10", "chr1\t0\t10"))
  err <- expect_error(merge_lines(c("chr1\t0\t10", "chr2\t0\t10", "chr1\t20\t30")),
                      class = "mytools_error")
  expect_match(conditionMessage(err), "in.bed:3: input is not sorted")
})

test_that("-header prints only the header lines before the first record", {
  expect_equal(merge_lines(c("#h", "track x", "chr1\t10\t20", "#mid", "chr1\t15\t30"),
                           header = TRUE),
               c("#h", "track x", "chr1\t10\t30"))
  expect_equal(merge_lines(c("#h", "chr1\t10\t20"), header = FALSE), "chr1\t10\t20")
})

test_that("empty input produces no output", {
  expect_equal(merge_lines(character(0)), character(0))
})

test_that("main_merge: -i is required, unknown flags are rejected, first -i wins", {
  r <- run_mytools("merge")
  expect_equal(r$status, 1L)
  expect_true(any(grepl("^mytools: merge: -i is required$", r$stderr)))
  expect_length(r$stdout, 0)

  r <- run_mytools(c("merge", "-i", write_fixture("chr1\t1\t2"), "-d", "5"))
  expect_equal(r$status, 1L)
  expect_true(any(grepl("^mytools: merge: unrecognized option -d$", r$stderr)))

  first <- write_fixture("chr1\t1\t2")
  second <- write_fixture("chr1\t5\t6")
  r <- run_mytools(c("merge", "-i", first, "-i", second))
  expect_equal(r$status, 0L)
  expect_equal(r$stdout, "chr1\t1\t2")
})
