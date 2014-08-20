#' A R6 class to store GRangesList and comparison matrix
#' 
chip_compare <- R6::R6Class("chip_compare",
  public = list(
    ## initialize
    initialize = function(grl1, grl2=NULL) {
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
    },
    ## get_matrix
    get_matrix = function() {
      private$score_matrix
    },
    ## print
    print = function(...) {
      message <- paste0("chip_compare Object of length: ",
                        length(private$grl), "\n")
      message <- paste0(message, "\n")
      cat(message)
      get_print_info <- function(gr_name, i) {
        cat(paste0("[", which(names(private$grl[[i]]) == gr_name), "] ",
                  gr_name, " (length: ", length(private$grl[[i]][[gr_name]]),
                  ")\n"))
      }
      for (i in 1:length(private$grl)) {
        cat(paste0("[[", i, "]] GRangesList of length ",
                   length(private$grl[[i]]),"\n"))
        invisible(lapply(names(private$grl[[i]]), get_print_info, i))
        cat("\n")
      }
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
