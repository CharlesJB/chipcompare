#' Fetch annotations
#'
#' \code{get_basic_annotations} fetch annotations from UCSC hg19 and from Fantom
#'
#' This function is meant to be run only to generate the GRangesList
#' object that will contain the general annotation. It should only be
#' used before building the package. The goal of this function is to
#' fetch the basic annotations from UCSC (genes, exons, TSS, TTR) and
#' enhancers from Fantom.
#'
#' @return A GRangesList object with a GRanges for every basic annotation.
#'
#' @examples
#' \dontrun{
#' build_basic_annotations()
#' }
build_basic_annotations <- function() {
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

  # 1. Fetch the GRanges objects
  genes <- GenomicFeatures::genes(txdb)
  exons <- GenomicFeatures::exons(txdb)
  transcripts <- GenomicFeatures::transcripts(txdb)
  TSS <- GenomicRanges::promoters(txdb, downstream = 1, upstream = 0)
  TTR <- transcripts
  end(TTR[strand(TTR) == "+"]) <- start(TTR[strand(TTR) == "+"])
  start(TTR[strand(TTR) == "-"]) <- end(TTR[strand(TTR) == "-"])

  enhancers_url <-
    "http://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/hg19_enhancers.bed.gz"
  download.file(enhancers_url, destfile="hg19_enhancer.bed.gz")
  enhancers <- read.table(gzfile("hg19_enhancer.bed.gz"), header=FALSE,
    stringsAsFactors=FALSE)
  colnames(enhancers) <- c("chrom", "chromStart", "chromEnd", "name", "score",
    "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes",
    "blockStarts")
  enhancers$strand <- "*"
  enhancers <- as(enhancers[,c(1:4,6)], "GRanges")
  file.remove("hg19_enhancer.bed.gz")

  # 2. Make sure all metadata fits (names and count must be the same
  #    if we want to group them as a single GRangesList).
  names(mcols(genes)) <- "id"
  mcols(genes)$id <- as.character(mcols(genes)$id)
  names(mcols(exons)) <- "id"
  mcols(TSS)$tx_id <- NULL
  names(mcols(TSS)) <- "id"
  mcols(TTR)$tx_id <- NULL
  names(mcols(TTR)) <- "id"
  names(mcols(enhancers)) <- "id"

  # 3. Build the GRangesList object
  basic_annotations <- GenomicRanges::GRangesList()
  basic_annotations$genes <- genes
  basic_annotations$exons <- exons
  basic_annotations$TSS <- TSS
  basic_annotations$TTR <- TTR
  basic_annotations$enhancers <- enhancers

  # 4. Save dataset
  save(basic_annotations, file = "data/basic_annotations.rda")
  invisible(basic_annotations)
}
