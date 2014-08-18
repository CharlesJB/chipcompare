#' A R6 class to store GRangesList and comparison matrix
#' 
chip_compare <- R6::R6Class("chip_compare",
  public = list(
    initialize = function(grlist1, grlist2=NULL) {
      stopifnot(class(grlist1) == "GRangesList")
      stopifnot(is.null(grlist2) | class(grlist2) == "GRangesList")
      private$grl[[1]] <- grlist1
      private$grl[[2]] <- grlist2
    },
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
    grl = list()
  )
)