#' A R6 class to store GRangesList and comparison matrix
#' 
chip_compare <- R6::R6Class("chip_compare",
  public = list(
    ## initialize
    initialize = function(grl1, grl2=NULL) {
      stopifnot(class(grl1) == "GRangesList")
      stopifnot(is.null(grl2) | class(grl2) == "GRangesList")
      private$grl[[1]] <- grl1
      private$grl[[2]] <- grl2
      private$score_matrix <- private$produce_matrix()
    },
    ## set_compute_score_FUN
    set_compute_score_FUN = function(FUN) {
      stopifnot(class(FUN(private$grl[[1]], private$grl[[1]])) == "numeric")
      private$compute_score <- FUN
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
      count <- length(unique(queryHits(findOverlaps(qry, sbj))))
      count / length(qry)
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
      result <- matrix(nrow = q_length, ncol = s_length)
      for (i in 1:q_length) {
        for (j in 1:s_length) {
          result[i,j] <- private$compute_score(query[i], subject[j])
        }
      }
      rownames(result) <- names(query)
      colnames(result) <- names(subject)
      result
    }
  )
)
