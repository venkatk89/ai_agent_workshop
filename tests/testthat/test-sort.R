# Unit tests for `sort` ordering (#5). Runs without bedtools.

test_that("chromosomes sort lexicographically: chr17 before chr7, chr10 before chr2", {
  chrom <- c("chr7", "chr17", "chr2", "chr10", "chrX", "chr1")
  start <- rep(0, length(chrom))
  expect_equal(chrom[sort_order(chrom, start)],
               c("chr1", "chr10", "chr17", "chr2", "chr7", "chrX"))
})

test_that("start is compared numerically, not as text", {
  chrom <- rep("chr1", 3)
  start <- c(1000, 20, 3)
  expect_equal(start[sort_order(chrom, start)], c(3, 20, 1000))
})

test_that("equal (chrom, start) keeps input order regardless of end", {
  # bedtools is stable on ties and does not fall back to end:
  # `chr1 10 50` before `chr1 10 20` in the input stays that way in the output.
  chrom <- c("chr1", "chr1", "chr1")
  start <- c(10, 10, 10)
  end   <- c(50, 20, 30)
  expect_equal(sort_order(chrom, start), c(1L, 2L, 3L))
  expect_equal(end[sort_order(chrom, start)], c(50, 20, 30))
})

test_that("chrom comparison is raw byte order, as in bedtools, not locale collation", {
  # bedtools gives: Chr1 < chr1 < chr1-x < chr10 < chr1_random < chrM
  chrom <- c("chr1_random", "Chr1", "chr1", "chrM", "chr10", "chr1-x")
  start <- rep(1, length(chrom))
  expect_equal(chrom[sort_order(chrom, start)],
               c("Chr1", "chr1", "chr1-x", "chr10", "chr1_random", "chrM"))
})

test_that("coordinates beyond 32-bit range sort numerically", {
  chrom <- c("chr1", "chr1")
  start <- c(3000000000, 5)
  expect_equal(sort_order(chrom, start), c(2L, 1L))
})

test_that("parse_sort_args accepts -i and -header in either order", {
  expect_equal(parse_sort_args(c("-i", "x.bed")), list(input = "x.bed", header = FALSE))
  expect_equal(parse_sort_args(c("-header", "-i", "-")), list(input = "-", header = TRUE))
  expect_equal(parse_sort_args(c("-i", "stdin", "-header")), list(input = "stdin", header = TRUE))
})
