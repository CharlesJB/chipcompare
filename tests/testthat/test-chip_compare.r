### PREPARE TESTS

# Prepare GRanges objects
bed_list <- list()
bed_list[[1]] <- data.frame(seqnames="chr2", start=3, end=6, strand="+")
bed_list[[2]] <- data.frame(seqnames="chr2", start=5, end=8, strand="+")
bed_list[[3]] <- data.frame(seqnames="chr1", start=3, end=6, strand="+")
bed_list[[4]] <- data.frame(seqnames="chr1", start=5, end=8, strand="+")
bed_list[[5]] <- data.frame(seqnames=c("chr1", "chr2"), start=c(5,20),
                            end=c(8,24), strand=c("+", "-"))
bed_list[[6]] <- data.frame(seqnames="chr3", start=1, end=5, strand="+")
bed_list[[7]] <- data.frame(seqnames="chr3", start=8, end=14, strand="-")

gr_list <- lapply(bed_list, function(x) as(x, "GRanges"))

# Prepare GRangesList objects
empty_grl <- GenomicRanges::GRangesList()
grl <- GenomicRanges::GRangesList(gr_list)
for (i in 1:length(grl)) {
  names(grl)[i] <- paste("gr", i, sep="")
}

### START TESTS
## initialize
test_that("chipcompare::initialize - Check param - Valid", {
  expect_that({chipcompare$new(grl, heatmap=FALSE); TRUE}, is_true())
  expect_that({chipcompare$new(grl, grl[1:4], heatmap=FALSE); TRUE}, is_true())
})

test_that("chipcompare::initialize - Check param - Invalid grl1", {
  expect_that(chipcompare$new(), throws_error())
  expect_that(chipcompare$new(1), throws_error())
  expect_that(chipcompare$new("a"), throws_error())
  expect_that(chipcompare$new(list(a=1, b=2)), throws_error())
})

test_that("chipcompare::initialize - Check param - Invalid grl2", {
  expect_that(chipcompare$new(grl, 1), throws_error())
  expect_that(chipcompare$new(grl, "a"), throws_error())
  expect_that(chipcompare$new(grl, list(a=1, b=2)), throws_error())
})

test_that("chipcompare::initialize - Check param - empty grl", {
  expect_that(chipcompare$new(empty_grl), throws_error())
  expect_that(chipcompare$new(grl, empty_grl), throws_error())
  expect_that(chipcompare$new(empty_grl, grl[1]), throws_error())
})

test_that("chipcompare::initialize - Check param - Not enough elements", {
  expect_that(chipcompare$new(grl[1]), throws_error())
})

## get_matrix
test_that("chipcompare::get_matrix - symetric matrix", {
  expect_that(chipcompare$new(grl[1:3], heatmap=FALSE)$get_matrix(),
              is_equivalent_to(matrix(c(1,1,0,1,1,0,0,0,1), nrow=3)))
})

test_that("chipcompare::get_matrix - asymetric matrix", {
  expect_that(chipcompare$new(grl[1:3], grl[4:5], heatmap=FALSE)$get_matrix(),
              is_equivalent_to(matrix(c(0,0,1,0,0,1), nrow=3)))
})

test_that("chipcompare::get_matrix - no overlap matrix", {
  expect_that(chipcompare$new(grl[1:3], grl[6:7], heatmap=FALSE)$get_matrix(),
              is_equivalent_to(matrix(0, ncol=2, nrow=3)))
})

test_that("chipcompare::get_matrix - partial overlap matrix", {
  expect_that(chipcompare$new(grl[5:7], grl[1:3], heatmap=FALSE)$get_matrix(),
              is_equivalent_to(matrix(c(0,0,1,0,0,0,0.5,0.0,0.0), nrow=3)))
})
