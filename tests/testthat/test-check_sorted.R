# check_sorted(): the unsorted-input check merge and closest rely on. bedtools
# applies two rules (see R/bed.R): closest wants lexicographic chrom order,
# merge only wants each chromosome contiguous.

feed <- function(recs, ...) {
  chk <- check_sorted("closest", "in.bed", ...)
  for (i in seq_along(recs)) chk(recs[[i]][[1]], as.numeric(recs[[i]][[2]]), i)
  invisible(TRUE)
}
unsorted <- function(recs, line, ...) {
  err <- expect_error(feed(recs, ...), class = "mytools_error")
  want <- sprintf("mytools: closest: in.bed:%d: input is not sorted by chrom then start", line)
  expect_equal(conditionMessage(err), want)
}

test_that("sorted input passes, including equal starts and ties on end", {
  expect_true(feed(list(list("chr1", 0), list("chr1", 0), list("chr1", 10),
                        list("chr10", 5), list("chr2", 0), list("chrX", 0))))
})

test_that("a start going backwards within a chromosome is caught at that line", {
  unsorted(list(list("chr1", 300), list("chr1", 0)), 2)
})

test_that("lexicographic: chr2 before chr10, or chrX before chr10, is out of order", {
  # Bytewise: "chr10" < "chr2" < "chrX". This is what bedtools closest demands.
  unsorted(list(list("chr2", 0), list("chr10", 0)), 2)
  unsorted(list(list("chrX", 0), list("chr10", 0)), 2)
  unsorted(list(list("chr2", 0), list("chr1", 0)), 2)
})

test_that("lexicographic order is bytewise regardless of locale", {
  # "chrM" < "chra" bytewise (uppercase sorts first); most UTF-8 locales disagree.
  expect_true(feed(list(list("chrM", 0), list("chra", 0))))
  unsorted(list(list("chra", 0), list("chrM", 0)), 2)
})

test_that("grouped: any chromosome order is fine but a revisit is not (bedtools merge)", {
  expect_true(feed(list(list("chr2", 0), list("chr1", 0), list("chr1", 5)),
                   chrom_order = "grouped"))
  unsorted(list(list("chr1", 0), list("chr2", 0), list("chr1", 7)), 3, chrom_order = "grouped")
  unsorted(list(list("chr1", 5), list("chr1", 0)), 2, chrom_order = "grouped")
})

test_that("read_bed_lines(sorted = TRUE) applies the check with real line numbers", {
  r <- read_all(c("#h", "chr1\t0\t10", "chr1\t5\t20", "chr1\t2\t3", "chr1\t9\t9"),
                sorted = TRUE, file = "a.bed")
  expect_match(conditionMessage(r$error), "a.bed:4: input is not sorted by chrom then start$")
  expect_length(r$records, 2)   # the out-of-order record never reaches on_record
  r <- read_all(c("chr2\t0\t1", "chr1\t0\t1"), sorted = TRUE, chrom_order = "grouped")
  expect_null(r$error)
})
