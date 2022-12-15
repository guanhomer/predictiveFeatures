#' @title Sample nucleotides or motifs on a subset of genome defined by \code{GRanges}.
#'
#' @description A function used to extract the regions of a query sequence mapped to user defined sub-regions of the genome.
#'
#' @param motif a character string indicating the query sequence or motifs; Vague mapping rules in \code{\link{IUPAC_CODE_MAP}} is supported when \code{Fixed} = FALSE.
#'
#' @param region the \code{\link{GRanges}} object used to define the subset region of the genome.
#'
#' @param sequence the \code{\link{BSgenome}} object containing the sequence of the genome.
#'
#' @param fixed FALSE to support the vague mapping rule of the character string, default is FALSE.
#'
#' @param N number of ranges sampled, by default, it returns all the matched ranges without sub-sampling.
#'
#' @param replace whether sample with replacement or not, default is FALSE.
#'
#' @return a \code{GRanges} object contains the (sampled) mapped regions of the query sequence on the given subset of the genome.
#'
#' @examples
#' ## Retrieve the dimer sequence "CG" on 3 exons of hg19
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' exons_gr <- GRanges(c("chr1","chr2","chr2"),
#'                     IRanges(start=c(11874,38814,45440),
#'                     end=c(12227,41627,46588)),
#'                     strand=c("+","-","-"))
#' CpG_regions <- sampleSequence("GG", exons_gr, Hsapiens)
#'
#' \donttest{
#' ## Sample 100000 ranges of DRACH motif (the motif of m6A) on all exon regions of hg19:
#' exons_hg19 <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' Motif_exons_hg19 <- sampleSequence("DRACH", exons_hg19, Hsapiens, N = 100000)
#' }
#'
#' @importFrom GenomicRanges reduce GRanges
#' @importFrom GenomicFeatures mapFromTranscripts
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom Biostrings vmatchPattern
#' @export
sampleSequence <- function(motif, region, sequence, N = NULL, fixed = FALSE, replace = FALSE){
  stopifnot(is(region, "GRangesList")|is(region, "GRanges"))
  if(is(region, "GRangesList")) region <- unlist(region)
  region <- reduce(region)

  region_dnass <- getSeq(x=sequence,
                         names=seqnames(region),
                         start=start(region),
                         end=end(region),
                         strand=strand(region),
                         as.character=FALSE)
  indx <- paste0("reg_", seq_along(region))
  regions_GRL <- split(region, indx)
  regions_GRL <- regions_GRL[indx]
  rm(indx)
  vmp <- vmatchPattern(motif, region_dnass, fixed = fixed)
  rm(region_dnass)
  vmp_gr <- GRanges(seqnames = rep(names(regions_GRL), elementNROWS(vmp)), ranges = unlist(vmp))
  rm(vmp)
  motif_on_regions <- mapFromTranscripts(vmp_gr,regions_GRL)
  rm(vmp_gr, regions_GRL)
  mcols(motif_on_regions) = NULL
  if(is.null(N)) N <- length(motif_on_regions)
  if(replace == FALSE) N <- min(N, length(motif_on_regions))
  indx2 <- sample.int(length(motif_on_regions), N, replace = replace)
  motif_on_regions <- motif_on_regions[indx2]
  rm(indx2)
  seqlengths(motif_on_regions) <- seqlengths(region)
  return(motif_on_regions)
}
