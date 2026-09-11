# open_bed(): "-" / "stdin" mean standard input, .gz is transparent, anything
# unreadable is `cannot open <file>`.

lines <- c("chr1\t0\t100", "chr1\t100\t200")

test_that("a plain file opens and reads", {
  path <- write_fixture(lines); on.exit(unlink(path))
  con <- open_bed(path); on.exit(close(con), add = TRUE)
  expect_true(isOpen(con))
  expect_equal(readLines(con), lines)
})

test_that("a .gz file is decompressed transparently", {
  path <- write_fixture(lines, gz = TRUE); on.exit(unlink(path))
  con <- open_bed(path); on.exit(close(con), add = TRUE)
  expect_s3_class(con, "gzfile")
  expect_equal(readLines(con), lines)
})

test_that("'-' and 'stdin' both name standard input", {
  for (name in c("-", "stdin")) {
    con <- open_bed(name)
    expect_equal(summary(con)$description, "stdin", info = name)
    expect_equal(bed_file_name(con), "stdin", info = name)
    close(con)
  }
})

test_that("'-' and 'stdin' actually read the process's standard input", {
  # Needs a fresh Rscript so we can pipe into it. Skipped if Rscript is absent.
  rscript <- file.path(R.home("bin"), "Rscript")
  skip_if_not(file.exists(rscript))
  bed_r <- normalizePath(testthat::test_path("..", "..", "R", "bed.R"))
  for (name in c("-", "stdin")) {
    code <- sprintf(
      'source("%s"); con <- open_bed("%s"); n <- 0; read_bed_lines(con, function(c, s, e, f) n <<- n + 1); cat(n)',
      bed_r, name)
    out <- system2(rscript, c("-e", shQuote(code)), input = lines, stdout = TRUE)
    expect_equal(out, "2", info = name)
  }
})

test_that("a missing file is `cannot open <file>`", {
  err <- expect_error(open_bed("does-not-exist.bed", cmd = "sort"), class = "mytools_error")
  expect_equal(conditionMessage(err), "mytools: sort: cannot open does-not-exist.bed")
  expect_equal(err$status, 1L)
})

test_that("a directory is `cannot open <file>`", {
  err <- expect_error(open_bed(tempdir(), cmd = "sort"), class = "mytools_error")
  expect_match(conditionMessage(err), "^mytools: sort: cannot open ")
})

test_that("cmd defaults to the running subcommand", {
  old <- options(mytools.cmd = "closest"); on.exit(options(old))
  err <- expect_error(open_bed("nope.bed"), class = "mytools_error")
  expect_equal(conditionMessage(err), "mytools: closest: cannot open nope.bed")
})
