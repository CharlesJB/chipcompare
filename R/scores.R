#' Calculate the percentage between 2 GRanges objects.
#' 
#' @param qry Query. The percentage will be calculated relatively to this
#'            object. Must of class \code{GRanges}.
#' @param sbj Subject. Must be of class \code{GRanges} or \code{GIntervalTree}.
#'
#' @return The percentage of the query that overlaps the subject.
#'
#' @examples
#' gr1 <- GRanges("chr1", IRanges(1,100))
#' gr2 <- GRanges("chr1", IRanges(40,140))
#' percent_overlap(gr1, gr2)
#' gri2 <- GIntervalTree(gr2)
#' percent_overlap(gr1, gr2)
#'
#' @export
percent_overlap <- function(qry, sbj) {
  overlaps <- findOverlaps(qry, sbj)
  count_qry <- length(unique(S4Vectors::queryHits(overlaps)))
  count_qry / length(qry)
}

#' Calculate the p-value associated with a given overlap percentage.
#' 
#' @param obs The observed percentage of overlap.
#' @param overlaps The pre-calculated overlaps between baseline and current
#'                 subject. Must be a logical vector of the same length as
#'                 baseline.
#' @param sample_size The number of value to pick for each draw.
#' @param sample_count The number of draw to calculate p-value.
#'
#' @return The p-value for current subject.
#'
#' @export
# TODO: add exemples
calculate_pvalue <- function(obs, overlaps, sample_size, sample_count) {
  i <- sample(1:length(overlaps), sample_size*sample_count, replace = TRUE)
  bootstrap <- colSums(matrix(overlaps[i], nrow = sample_size)) / sample_size
  sum(bootstrap > obs)/length(bootstrap)
}
