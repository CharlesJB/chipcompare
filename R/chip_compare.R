#' A R6 class to store GRangesList and comparison matrix
#' 
chip_compare <- R6::R6Class("chip_compare",
  public = list(
    initialize = function(grlist) {
      stopifnot(class(grlist) == "GRangesList")
      private$grl <- grlist
    },
    print = function(...) {
      message <- paste0("ChipCompare Object of length: ", length(private$grl), "\n")
      message <- paste0(message, "\n")
      cat(message)
      get_print_info <- function(gr_name) {
        cat(paste0("[", which(names(private$grl) == gr_name), "] ",
                  gr_name, " (length: ", length(private$grl[[gr_name]]), ")\n"))
      }
      invisible(lapply(names(private$grl), get_print_info))
    },
    getLength = function() {
      length(private$grl)
    }
  ),
  private = list(
    grl = GenomicRanges::GRangesList()
  )
)