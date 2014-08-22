#' chipcompare.
#'
#' @name chipcompare
#' @docType package
NULL

#' Annotations of basic genetic elements.
#'
#' This dataset is a GRangesList of multiple genomic features.
#'   * genes
#'   * exons
#'   * TSS
#'   * TTR
#'   * enhancers
#' Built with the build_basic_annotations function.
#'
#' \itemize{
#'   \item genes. TxDb.Hsapiens.UCSC.hg19.knownGene
#'   \item exons. TxDb.Hsapiens.UCSC.hg19.knownGene
#'   \item TSS. TxDb.Hsapiens.UCSC.hg19.knownGene
#'   \item TTR. TxDb.Hsapiens.UCSC.hg19.knownGene
#'   \item enhancers. Fantom5
#'   ...
#' }
#'
#' @format A GRangesList with 5 GRanges elements of 23056, 289969, 82960, 82960
#'   and 43011 elements.
#' @source Bioconductor TxDb.Hsapiens.UCSC.hg19.knownGene and
#'   \url{http://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/hg19_enhancers.bed.gz}
#' @name basic_annotations
NULL
