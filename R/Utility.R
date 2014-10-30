get_collapsed_annotations <- function(grl) {
  # 1. Add GRanges name to id
  grl <- lapply(names(grl), function(x) { a <- grl[[x]]; mcols(a)$id <- paste(x, mcols(a)$id, sep = "_"); a})
  grl <- GRangesList(grl)
  
  # 2. Unlist GRangesList object
  BiocGenerics:::unlist(grl)
}

get_scores <- function(bam1, bam2, metric) {
	stopifnot(length(bam1) == length(bam2))
	scores <- numeric(length(bam1))
	factory <- MetricFactory$new(diffPosMaxTolerance=0.05)

	for (i in 1:length(bam1)) {
		vbam1 <- as.numeric(bam1[[i]])
		vbam2 <- as.numeric(bam2[[i]])
		vbam1[vbam1<0] <- 0
		vbam2[vbam2<0] <- 0
		scores[i] <- factory$createMetric(metricType=metric, profile1 = vbam1,
									profile2 = vbam2)
	}
	scores
}

calculate_similarity_matrix <- function(coverages, bam_files, metric, BPPARAM=NULL) {
  # Initialize result matrix
  ncol <- length(coverages[[1]][[1]])
  bf_count <- length(bam_files)
  nrow <- (bf_count*(bf_count-1))/2
  result <- matrix(ncol = ncol, nrow = nrow)
	rownames(result) <- rep("NA", nrow)

	# Calculate the scores
	current_row <- 1
  for (i in 1:(length(bam_files)-1)) {
    for (j in (i+1):length(bam_files)) {
      bam1 <- bam_files[i]
      bam2 <- bam_files[j]
      current_name <- file_path_sans_ext(basename(bam1))
      current_name <- paste(current_name, file_path_sans_ext(basename(bam2)),
                        sep = "_")
      print(current_name)

      bam1 <- coverages[[bam1]][[1]]
      bam2 <- coverages[[bam2]][[1]]
      if (!is.null(BPPARAM)) {
        max <- ceiling(length(bam1)/bpworkers(BPPARAM))
        bam1 <- split(bam1, ceiling(seq_along(bam1)/max))
        bam2 <- split(bam2, ceiling(seq_along(bam2)/max))
				metrics <- rep(metric, length(bam1))
				metrics <- split(metrics, seq_along(metrics))
        scores <- BiocParallel:::bpmapply(get_scores, bam1, bam2, metrics,
                    BPPARAM = BPPARAM)
				scores <- do.call("c", lapply(scores, unlist))
      } else {
        scores <- get_scores(bam1, bam2, metric)
      }
			result[current_row,] <- scores
			rownames(result)[current_row] <- current_name
      current_row <- current_row + 1
		}
	}
	result
}
