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

test_that("Initialize chip_compare - Check param", {
  expect_that(chip_compare$new(), throws_error())
  expect_that(chip_compare$new(1), throws_error())
  expect_that(chip_compare$new("a"), throws_error())
  expect_that(chip_compare$new(list(a=1, b=2)), throws_error())
})

test_that("getLength", {
  expect_that(empty_cc$getLength(), equals(0))
  expect_that(cc1$getLength(), equals(1))
  expect_that(cc2$getLength(), equals(2))
})