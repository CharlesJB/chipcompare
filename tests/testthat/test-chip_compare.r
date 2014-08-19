### PREPARE TESTS

# Prepare GRanges objects
bed_list <- list()
bed_list[[1]] <- data.frame(seqnames="chr2", start=3, end=6, strand="+")
bed_list[[2]] <- data.frame(seqnames="chr2", start=5, end=8, strand="+")
bed_list[[3]] <- data.frame(seqnames="chr1", start=3, end=6, strand="+")
bed_list[[4]] <- data.frame(seqnames="chr1", start=5, end=8, strand="+")
bed_list[[5]] <- data.frame(seqnames=c("chr1", "chr2"), start=c(5,20),
                            end=c(8,24), strand=c("+", "-"))

gr_list <- lapply(bed_list, function(x) as(x, "GRanges"))

# Prepare GRangesList objects
empty_grl <- GenomicRanges::GRangesList()
grl <- GenomicRanges::GRangesList(gr_list)
for (i in 1:length(grl)) {
  names(grl)[i] <- paste("gr", i, sep="")
}

### START TESTS
## initialize
test_that("chip_compare::initialize - Check param - Valid", {
  expect_that({chip_compare$new(grl[1]); TRUE}, is_true())
  expect_that({chip_compare$new(grl[1], grl[2]); TRUE}, is_true())
})

test_that("chip_compare::initialize - Check param - Invalid grl1", {
  expect_that(chip_compare$new(), throws_error())
  expect_that(chip_compare$new(1), throws_error())
  expect_that(chip_compare$new("a"), throws_error())
  expect_that(chip_compare$new(list(a=1, b=2)), throws_error())
})

test_that("chip_compare::initialize - Check param - Invalid grl2", {
  expect_that(chip_compare$new(grl[1], 1), throws_error())
  expect_that(chip_compare$new(grl[1], "a"), throws_error())
  expect_that(chip_compare$new(grl[1], list(a=1, b=2)), throws_error())
})

test_that("chip_compare::initialize - Check param - empty grl", {
  expect_that(chip_compare$new(empty_grl), throws_error())
  expect_that(chip_compare$new(grl[1], empty_grl), throws_error())
  expect_that(chip_compare$new(empty_grl, grl[1]), throws_error())
})

## get_matrix
test_that("chip_compare::get_matrix - symetric matrix", {
  expect_that(chip_compare$new(grl[1])$get_matrix(),
              is_equivalent_to(matrix(1, ncol=1, nrow=1)))
})

test_that("chip_compare::get_matrix - asymetric matrix", {
  expect_that(chip_compare$new(grl[1], grl[1:2])$get_matrix(),
              is_equivalent_to(matrix(1, ncol=2, nrow=1)))
})

test_that("chip_compare::get_matrix - no overlap matrix", {
  expect_that(chip_compare$new(grl[1], grl[3])$get_matrix(),
              is_equivalent_to(matrix(0, ncol=1, nrow=1)))
})

test_that("chip_compare::get_matrix - partial overlap matrix", {
  expect_that(chip_compare$new(grl[5], grl[3])$get_matrix(),
              is_equivalent_to(matrix(0.5, ncol=1, nrow=1)))
})
