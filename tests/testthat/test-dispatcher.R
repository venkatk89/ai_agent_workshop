# ./mytools dispatcher: usage and exit codes. Runs the real script.

test_that("--version prints the version on stdout and exits 0", {
  r <- run_mytools("--version")
  expect_equal(r$status, 0L)
  expect_match(r$stdout, "^mytools [0-9]+\\.[0-9]+\\.[0-9]+$")
  expect_length(r$stderr, 0)
})

test_that("no arguments gives usage on stderr and exit 1", {
  r <- run_mytools(character(0))
  expect_equal(r$status, 1L)
  expect_length(r$stdout, 0)
  expect_match(r$stderr[1], "^usage: mytools")
})

test_that("an unknown subcommand gives usage on stderr and exit 1", {
  r <- run_mytools("nope")
  expect_equal(r$status, 1L)
  expect_length(r$stdout, 0)
  expect_equal(r$stderr[1], "mytools: unknown subcommand nope")
  expect_true(any(grepl("^usage: mytools", r$stderr)))
})

test_that("die() builds `mytools: <cmd>: [<file>:<line>: ]<msg>` and carries the status", {
  e <- expect_error(die("sort", "boom", "x.bed", 7), class = "mytools_error")
  expect_equal(conditionMessage(e), "mytools: sort: x.bed:7: boom")
  expect_equal(e$status, 1L)
  e <- expect_error(die("sort", "cannot open x.bed"), class = "mytools_error")
  expect_equal(conditionMessage(e), "mytools: sort: cannot open x.bed")
  e <- expect_error(die("sort", "boom", "x.bed"), class = "mytools_error")
  expect_equal(conditionMessage(e), "mytools: sort: x.bed: boom")
  e <- expect_error(die("sort", "boom", status = 2), class = "mytools_error")
  expect_equal(e$status, 2L)
})

test_that("a known subcommand without an implementation is reported, not a crash", {
  # Only meaningful until R/<cmd>.R exists for every subcommand; skip once it does.
  missing <- SUBCOMMANDS[!file.exists(testthat::test_path("..", "..", "R", paste0(SUBCOMMANDS, ".R")))]
  skip_if(length(missing) == 0, "all subcommands implemented")
  r <- run_mytools(c(missing[1], "-i", "x.bed"))
  expect_equal(r$status, 1L)
  expect_equal(r$stderr, sprintf("mytools: %s: not implemented yet", missing[1]))
})
