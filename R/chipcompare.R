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
#'   \item{cores:}{Number of cores for parallel processing. Default: 1}
#'   \item{FUN:}{The function to compute the scores. Must take 2 \code{GRanges}
#'               as input and return a numeric as output.}
#'   \item{...:}{Extra options to pass to \code{FUN}.}
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
#' @importFrom parallel mclapply
#' @export
#' @format An overlap calculator
chipcompare <- R6::R6Class("chipcompare",
  public = list(
    ## initialize
    initialize = function(grl1, grl2=NULL, heatmap=TRUE, cores = 1,
			  FUN = percent_overlap, ...) {
      # Check parameters
      stopifnot(class(grl1) == "GRangesList")
      stopifnot(length(grl1) > 1)
      if (!is.null(grl2)) {
        stopifnot(class(grl2) == "GRangesList")
        stopifnot(length(grl2) > 1)
      }
      stopifnot(is.numeric(cores))
      stopifnot(cores > 0)
      stopifnot(as.integer(cores) == cores)
      stopifnot(is.function(FUN))
      # Initialize
      private$cores <- cores
      private$grl[[1]] <- grl1
      private$grl[[2]] <- grl2
      query <- private$grl[[1]]
      if (length(private$grl) == 2) {
        subject <- private$grl[[2]]
      } else {
        subject <- private$grl[[1]]
      }
      # If we convert subject to list GIntervalTree, the code is ~10X faster
      subject_names <- names(subject)
      subject <- lapply(1:length(subject),
                        function(x) GenomicRanges::GIntervalTree(subject[[x]]))
      names(subject) <- subject_names
      # Calculate scores
      private$score_matrix <- private$produce_matrix(query, subject,
                                FUN, ...)
      # Show heatmap
      if (heatmap == TRUE) {
        self$print()
      }
    },
    ## get_matrix
    get_matrix = function() {
      private$score_matrix
    },
    ## print
    print = function(...) {
      gplots::heatmap.2(private$score_matrix, col = redgreen(75), trace = "none", ...)
    },
    pvalue = function(baseline, sample_count = 1000) {
      # 0. Check params
      stopifnot(class(baseline) == "GRanges")
      stopifnot(length(baseline) > 1)
      stopifnot(is.numeric(sample_count))
      stopifnot(sample_count > 0)

      query <- private$grl[[1]]
      if (length(private$grl) == 2) {
        subject <- private$grl[[2]]
      } else {
        subject <- private$grl[[1]]
      }

      # 1. Pre-calculate the overlaps
      add_overlaps <- function(baseline, subject_name) {
        # Make sure the seqlevels fit
        sbj <- subject[[subject_name]]
        levels <- as.character(GenomeInfoDb::seqlevels(sbj))
        levels <- c(levels, as.character(GenomeInfoDb::seqlevels(baseline)))
        levels <- unique(levels)
        GenomeInfoDb::seqlevels(baseline) <- levels

        # Calculate the overlaps
        overlaps <- GenomicRanges::overlapsAny(baseline, sbj)
        GenomicRanges::elementMetadata(baseline)[subject_name] <- overlaps
        baseline
      }
      for (name in names(subject)) {
        baseline <- add_overlaps(baseline, name)
      }

      # 2. Produce the matrix
      launch_pval_calculations <- function(i, sbj_name) {
        j <- which(names(subject) == sbj_name)
        obs <- self$get_matrix()[i,j]
        overlaps <- GenomicRanges::elementMetadata(baseline)[[sbj_name]]
        sample_size <- length(query[[i]])
        chipcompare:::calculate_pvalue(obs = obs, overlaps = overlaps,
                                       sample_size = sample_size,
                                       sample_count = sample_count)
      }
      do.call("rbind", lapply(seq_along(query), function(i) {
        unlist(mclapply(names(subject), function(sbj_name) {
          launch_pval_calculations(i, sbj_name)
        }, mc.cores = private$cores)
      )}))
    }
  ),
  private = list(
    #### private members
    grl = list(),
    score_matrix = matrix(),
    cores = NULL,

    #### private methods
    ## produce_matrix
    produce_matrix = function(query, subject, FUN, ...) {
      do.call("rbind", lapply(query, function(qry) {
        unlist(mclapply(subject, function(sbj) {
          FUN(qry, sbj, ...)
        }, mc.cores = private$cores)
      )}))
    }
  )
)
