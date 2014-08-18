# Prepare GRanges objects
bed1 <- data.frame(seqnames="chr2", start=3, end=6, strand="+")
bed2 <- data.frame(seqnames="chr2", start=5, end=8, strand="+")

gr1 <- as(bed1, "GRanges")
gr2 <- as(bed2, "GRanges")

# Prepare GRangesList objects
empty_grl <- GenomicRanges::GRangesList()
grl1 <- GenomicRanges::GRangesList("gr1" = gr1)
grl2 <- GenomicRanges::GRangesList("gr1" = gr1, "gr2" = gr2)

# Prepare chip_compare objects
empty_cc <- chip_compare$new(empty_grl)
cc1 <- chip_compare$new(grl1)
cc2 <- chip_compare$new(grl2)

test_that("Initialize chip_compare - Check param - Valid", {
  expect_that({chip_compare$new(grl1); TRUE}, is_true())
  expect_that({chip_compare$new(grl1, grl2); TRUE}, is_true())
})

test_that("Initialize chip_compare - Check param - Invalid grl1", {
  expect_that(chip_compare$new(), throws_error())
  expect_that(chip_compare$new(1), throws_error())
  expect_that(chip_compare$new("a"), throws_error())
  expect_that(chip_compare$new(list(a=1, b=2)), throws_error())
})

test_that("Initialize chip_compare - Check param - Invalid grl2", {
  expect_that(chip_compare$new(grl1, 1), throws_error())
  expect_that(chip_compare$new(grl1, "a"), throws_error())
  expect_that(chip_compare$new(grl1, list(a=1, b=2)), throws_error())
})
