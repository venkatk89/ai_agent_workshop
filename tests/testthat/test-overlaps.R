# overlaps(): the one predicate every off-by-one hides in. BED is 0-based,
# half-open, so `chr1 100 200` covers bases 100..199.

test_that("plainly overlapping intervals overlap", {
  expect_true(overlaps(100, 200, 150, 250))
  expect_true(overlaps(150, 250, 100, 200))   # symmetric
  expect_true(overlaps(100, 200, 199, 300))   # share exactly one base (199)
})

test_that("disjoint intervals do not overlap", {
  expect_false(overlaps(100, 200, 300, 400))
  expect_false(overlaps(300, 400, 100, 200))
})

test_that("bookended intervals do not overlap (a.end == b.start)", {
  # a01 chr1 0-100 and a02 chr1 100-200 share no base: 100 belongs to a02 only.
  expect_false(overlaps(0, 100, 100, 200))
  expect_false(overlaps(100, 200, 0, 100))
})

test_that("nested intervals overlap", {
  # a06 sits inside a05; either way round.
  expect_true(overlaps(300, 400, 320, 380))
  expect_true(overlaps(320, 380, 300, 400))
  expect_true(overlaps(300, 400, 300, 400))   # identical coordinates (a09/a10)
})

test_that("zero-length interval versus an enclosing interval", {
  # The plain predicate: a07 `chr1 500 500` is inside 400-600 (500 < 600 and
  # 400 < 500) but does not touch an interval that starts at 500 or ends at 500.
  # bedtools widens zero-length records to [499, 501) before testing, so the
  # subcommands, not this predicate, are where the oracle's answer is matched.
  expect_true(overlaps(500, 500, 400, 600))
  expect_true(overlaps(400, 600, 500, 500))
  expect_false(overlaps(500, 500, 500, 600))
  expect_false(overlaps(500, 500, 400, 500))
  expect_false(overlaps(500, 500, 500, 500))
})

test_that("interval at position 0 behaves like any other", {
  expect_true(overlaps(0, 100, 0, 1))
  expect_true(overlaps(0, 100, 50, 150))
  expect_true(overlaps(0, 1, 0, 100))
  expect_false(overlaps(0, 100, 100, 101))
  expect_false(overlaps(0, 0, 0, 100))        # zero-length at 0 (a12: chr2 0 0)
})

test_that("overlaps() is vectorised", {
  expect_equal(overlaps(100, 200, c(150, 200, 50), c(250, 300, 100)),
               c(TRUE, FALSE, FALSE))
})
