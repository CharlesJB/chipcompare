# A R6 class to store GRangesList and comparison matrix
chip_compare <- R6::R6Class("chip_compare",
  public = list(
    ## initialize
    initialize = function(grl1, grl2=NULL, heatmap=TRUE, ...) {
      # Check parameters
      stopifnot(class(grl1) == "GRangesList")
      stopifnot(length(grl1) > 0)
      if (!is.null(grl2)) {
        stopifnot(class(grl2) == "GRangesList")
        stopifnot(length(grl2) > 0)
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
chip_compare_dummy <- R6::R6Class("chip_compare_dummy",
  inherit = chip_compare,
  private = list(
    compute_score = function(qry, sbj) {
      c(0.5, 0.5)
    }
  )
)
