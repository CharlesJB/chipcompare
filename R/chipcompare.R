#' A R6 class to store GRangesList and comparison matrix
#'
#' This class allows to produce an overlap matrix of multiple genomic regions
#' and it's visual representation in the heatmap format. For each pairs of
#' regions, this object will calculate the percentage of overlapping regions.
#'
#' @section Constructor:
#' \describe{
#'   \item{}{\code{cc <- chipcompare$new(grl1, grl2 = NULL, heatmap=TRUE,
#'                                       ...)}}
#'   \item{grl1:}{A \code{GRangesList} object.}
#'   \item{grl2:}{A \code{GRangesList} object. If value is \code{NULL}, all the
#'               regions in \code{grl1} will be compared against themselves.
#'               Default: \code{NULL}.}
#'   \item{heatmap:}{Print heatmap. Default: \code{TRUE}}
#'   \item{...:}{Extra options to pass to \code{heatmap.2} function.}
#' }
#'
#' @return
#'  \code{chipcompare$new} returns a \code{chipcompare} object that contains the
#'  \code{GRangesList} object(s) used to produce the overlap matrix and the
#'  overlap matrix.
#'
#' @section Methods:
#' \describe{
#'   \item{}{\code{m <- cc$get_matrix()}}
#' }
#' \describe{
#'   \item{}{\code{cc$print(...)}}
#'   \item{...:}{Extra options to pass to \code{heatmap.2} function.}
#' }
#' \describe{
#'   \item{}{\code{pval <- cc$pvalue(baseline, sample_count = 1000)}}
#'   \item{baseline:}{A \code{GRanges} object to use as backgroud.}
#'   \item{sample_count:}{The number of draw to use for bootstrap.}
#' }
#'
#' @examples
#'  # Prepare the GRangesList object
#'  gr1 <- GRanges("chr1", IRanges(c(1,100),c(10,110)))
#'  gr2 <- GRanges("chr1", IRanges(c(2,200),c(8,210)))
#'  gr3 <- GRanges(c("chr1", "chr2"), IRanges(c(1,100),c(10,110)))
#'  gr4 <- GRanges("chr3", IRanges(1,1000))
#'  grl <- GRangesList(gr1, gr2, gr3, gr4)
#'
#'  # Create a chipcompare object
#'  cc <- chipcompare$new(grl, heatmap = FALSE)
#'
#'  # Extract the overlap matrix
#'  m <- cc$get_matrix()
#'  mg <- metagene$new(regions, bam_files)
#'
#'  # Draw a heatmap
#'  \dontrun{cc$print()}
#'  \dontrun{cc}
#'
#'  # Calculate p-values. In this case we would need a GRanges object that
#'  # represent baseline values. For example, it could be a list of every
#'  # regions that were detected in the current cell line with other ChIP-Seq
#'  # experiments.
#'  \dontrun{pval <- cc$pvalue(baseline)}
#'
#' @importFrom R6 R6Class
#' @importFrom GenomicRanges findOverlaps
#' @export
#' @format An overlap calculator
chipcompare <- R6::R6Class("chipcompare",
  public = list(
    ## initialize
    initialize = function(grl1, grl2=NULL, heatmap=TRUE, ...) {
      # Check parameters
      stopifnot(class(grl1) == "GRangesList")
      stopifnot(length(grl1) > 1)
      if (!is.null(grl2)) {
        stopifnot(class(grl2) == "GRangesList")
        stopifnot(length(grl2) > 1)
      }
      # Initialize
      private$grl[[1]] <- grl1
      private$grl[[2]] <- grl2
      private$score_matrix <- private$produce_matrix()
      # Show heatmap
      if (heatmap == TRUE) {
        self$print(...)
      }
    },
    ## get_matrix
    get_matrix = function() {
      private$score_matrix
    },
    ## print
    print = function(...) {
      heatmap.2(private$score_matrix, col = redgreen(75), trace = "none", ...)
    }
  ),
  private = list(
    #### private members
    grl = list(),
    score_matrix = matrix(),

    #### private methods
    ## compute_score
    compute_score = function(qry, sbj) {
      overlaps <- findOverlaps(qry, sbj)
      count_qry <- length(unique(queryHits(overlaps)))
      count_sbj <- length(unique(subjectHits(overlaps)))
      percent_qry <- count_qry / length(qry)
      percent_sbj <- count_sbj / length(sbj)
      c(percent_qry, percent_sbj)
    },
    ## produce_matrix
    produce_matrix = function() {
      query <- private$grl[[1]]
      if (length(private$grl) == 2) {
        subject <- private$grl[[2]]
      } else {
        subject <- private$grl[[1]]
      }
      # If we convert subject to list GIntervalTree, the code is ~10X faster
      subject_names <- names(subject)
      subject <- lapply(1:length(subject),
                        function(x) GenomicRanges:::GIntervalTree(subject[[x]]))
      names(subject) <- subject_names
      # We start the calculations
      q_length <- length(query)
      s_length <- length(subject)
      min_length <- min(q_length, s_length)
      result <- matrix(nrow = q_length, ncol = s_length)
      for (i in 1:q_length) {
        for (j in 1:s_length) {
          if (is.na(result[i,j])) {
            current_result <- private$compute_score(query[[i]], subject[[j]])
            result[i,j] <- current_result[1]
            if (i != j & j <= q_length & i <= s_length) {
              result[j,i] <- current_result[2]
            }
          }
        }
      }
      rownames(result) <- names(query)
      colnames(result) <- names(subject)
      result
    }
  )
)

# This shows how it is possible to use inheritance to use and alternative
# compute_score function
chipcompare_dummy <- R6::R6Class("chipcompare_dummy",
  inherit = chipcompare,
  private = list(
    compute_score = function(qry, sbj) {
      c(0.5, 0.5)
    }
  )
)
