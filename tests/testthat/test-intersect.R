# Unit tests for `intersect` (#7). Runs without bedtools. Every expectation here
# was checked against bedtools v2.31.1 first; comments say where it surprises.

# Build a -b index from BED lines the way load_b() does, without touching files.
b_index <- function(lines) {
  path <- write_fixture(lines)
  on.exit(unlink(path))
  load_b(path)
}
hits <- function(a, b_lines) {
  f <- strsplit(a, "\t", fixed = TRUE)[[1]]
  intersect_one(f[1], as.numeric(f[2]), as.numeric(f[3]), f, b_index(b_lines))
}

# --- the five cases the issue asks for -----------------------------------------

test_that("bookended intervals produce no output", {
  expect_equal(hits("chr1\t0\t100", "chr1\t100\t200"), character(0))
  expect_equal(hits("chr1\t100\t200", "chr1\t0\t100"), character(0))
})

test_that("a nested -b interval reports the inner interval", {
  expect_equal(hits("chr1\t300\t400", "chr1\t320\t350"), "chr1\t320\t350")
  # and the other way round: -a inside -b reports -a
  expect_equal(hits("chr1\t320\t350", "chr1\t300\t400"), "chr1\t320\t350")
})

test_that("a partial overlap reports the clipped portion", {
  expect_equal(hits("chr1\t100\t200", "chr1\t150\t250"), "chr1\t150\t200")
  expect_equal(hits("chr1\t150\t250", "chr1\t100\t200"), "chr1\t150\t200")
})

test_that("a chrom present only in -b yields nothing", {
  expect_equal(hits("chr1\t0\t100", "chr3\t0\t100"), character(0))
  expect_equal(hits("chr1\t0\t100", c("chr2\t0\t100", "chr3\t50\t60")), character(0))
})

test_that("extra -a columns pass through after the clipped coordinates", {
  expect_equal(hits("chr1\t100\t200\ta02\t20\t-", "chr1\t150\t250\tb03\t12\t+"),
               "chr1\t150\t200\ta02\t20\t-")
})

# --- more than one hit ---------------------------------------------------------

test_that("one -a record overlapping three -b records prints three lines", {
  out <- hits("chr1\t0\t1000\tq", c("chr1\t500\t600", "chr1\t100\t200", "chr1\t700\t800"))
  expect_length(out, 3L)
})

test_that("hits in the same 16 kb bin come out in -b file order, not by start", {
  out <- hits("chr1\t0\t1000\tq", c("chr1\t500\t600", "chr1\t100\t200", "chr1\t700\t800"))
  expect_equal(out, c("chr1\t500\t600\tq", "chr1\t100\t200\tq", "chr1\t700\t800\tq"))
})

test_that("hits are ordered by bedtools' bin tree: finest level, then bin, then file order", {
  # bedtools prints y (level 0, bin 0), x (level 0, bin 1), p (level 0, bin 9),
  # then w and v (both 128 kb-level bin 0, file order) -- NOT file order p,x,y,w,v.
  out <- hits("chr1\t0\t200000\tq",
              c("chr1\t150000\t150100\tp", "chr1\t20000\t20100\tx", "chr1\t100\t200\ty",
                "chr1\t0\t20000\tw", "chr1\t100\t20100\tv"))
  expect_equal(out, c("chr1\t100\t200\tq", "chr1\t20000\t20100\tq", "chr1\t150000\t150100\tq",
                      "chr1\t0\t20000\tq", "chr1\t100\t20100\tq"))
})

test_that("bin_of places records in the finest bin that holds them whole", {
  expect_equal(bin_of(100, 200), c(0, 0))
  expect_equal(bin_of(20000, 20100), c(0, 1))
  expect_equal(bin_of(16384, 16400), c(0, 1))
  expect_equal(bin_of(16300, 16384), c(0, 0))        # end is exclusive
  expect_equal(bin_of(0, 20000), c(1, 0))            # straddles two 16 kb bins
  expect_equal(bin_of(500000, 600000), c(2, 0))
  expect_equal(bin_of(3000000000, 3000000100), c(0, 183105))
})

# --- zero-length intervals: encode the oracle ----------------------------------

test_that("a zero-length -b record is widened by 1 bp each side and reported widened", {
  # bedtools: -a chr1 400 600 vs -b chr1 500 500 prints chr1 499 501
  expect_equal(hits("chr1\t400\t600", "chr1\t500\t500"), "chr1\t499\t501")
  expect_equal(hits("chr1\t400\t500", "chr1\t500\t500"), "chr1\t499\t500")  # bookended, yet a hit
  expect_equal(hits("chr1\t500\t600", "chr1\t500\t500"), "chr1\t500\t501")
  expect_equal(hits("chr1\t400\t499", "chr1\t500\t500"), character(0))
  expect_equal(hits("chr1\t501\t600", "chr1\t500\t500"), character(0))
})

test_that("a zero-length -a record hits bookended -b records and is printed unwidened", {
  # bedtools: -a chr1 500 500 vs -b chr1 400 500 prints chr1 500 500
  expect_equal(hits("chr1\t500\t500", "chr1\t400\t500"), "chr1\t500\t500")
  expect_equal(hits("chr1\t500\t500", "chr1\t500\t600"), "chr1\t500\t500")
  expect_equal(hits("chr1\t500\t500", "chr1\t400\t499"), character(0))
  expect_equal(hits("chr1\t500\t500", "chr1\t501\t600"), character(0))
})

test_that("two zero-length records hit when 1 bp apart but not 2", {
  expect_equal(hits("chr1\t500\t500", "chr1\t500\t500"), "chr1\t500\t500")
  expect_equal(hits("chr1\t500\t500", "chr1\t501\t501"), "chr1\t500\t500")
  expect_equal(hits("chr1\t500\t500", "chr1\t502\t502"), character(0))
})

test_that("a zero-length -a record at position 0 works; a zero-length -b record there is fatal", {
  expect_equal(hits("chr1\t0\t0", "chr1\t0\t10"), "chr1\t0\t0")
  expect_equal(hits("chr1\t0\t0", "chr1\t1\t10"), character(0))
  # bedtools widens `chr1 0 0` in -b to [-1, 1), cannot bin it, and exits 1 with no
  # output. a.bed has such a record (a12), so `-b a.bed` always fails.
  expect_null(bin_of(-1, 1))
  expect_error(b_index(c("chr1\t5\t10", "chr1\t0\t0")), class = "mytools_error")
})

test_that("an interval at position 0 overlaps as usual", {
  expect_equal(hits("chr1\t0\t100\ta01", "chr1\t0\t50\tb01"), "chr1\t0\t50\ta01")
})

# --- argument parsing ----------------------------------------------------------

test_that("parse_intersect_args needs -a and -b, accepts -header anywhere, rejects others", {
  expect_equal(parse_intersect_args(c("-a", "x", "-b", "y")),
               list(a = "x", b = "y", header = FALSE))
  expect_equal(parse_intersect_args(c("-header", "-b", "y", "-a", "-")),
               list(a = "-", b = "y", header = TRUE))
  # the usage line goes to stderr alongside the error; keep it out of the test log
  quiet <- function(expr) capture.output(expr, type = "message")
  quiet(expect_error(parse_intersect_args(c("-a", "x")), "-b is required", class = "mytools_error"))
  quiet(expect_error(parse_intersect_args(c("-b", "y")), "-a is required", class = "mytools_error"))
  quiet(expect_error(parse_intersect_args(c("-a", "x", "-b", "y", "-v")), "unrecognized option -v",
                     class = "mytools_error"))
  quiet(expect_error(parse_intersect_args(c("-a", "-", "-b", "stdin")), class = "mytools_error"))
})
