#' @title Extract sequence-derived predictive features from interval-based data
#' @description A function to extract sequence features from the input \code{\link{GRanges}} object and the \code{\link{BSgenome}} object.
#'
#' @param x A \link{GRanges} object for the genomic ranges to be annotated, the \code{width} of x must all be equal.
#' @param sequence A \link{BSgenome} or \link{XStringSet} object for the genome sequence.
#' @param mutation A \link{GRanges} object indicating the SNP ranges (must be single-based),
#' a metadata column is required to annotate the changed nucleotide sequence.
#' @param encoding Can be one of the following:
#' \describe{
#' \item{\strong{onehot}}{From the 5' end to the 3'end of x,
#' each nucleotide position is coded by 4 indicators/dummy variables,
#' where each dummy variable indicates that the position is equal to the base "A", "T", "C", and "G", respectively.}
#'
#' \item{\strong{iRNA}}{
#' Each nucleotide position is encoded by 4 variables,
#' the first variable indicates that the nucleotide is purine (A or G),
#' the second variable indicates that the nucleotide has an amino group (A or C),
#' and the third variable indicates the formation of weak hydrogen bond (A or T),
#' the fourth variable is calculated by the cumulative frequency of nucleotides from
#' the leftmost position to the current position.
#' }
#'
#' \item{\strong{PseDNC}}{
#' This encoding method calculates the occurrence frequency of all 16 possible dinucleotides 
#' (e.g., AA, AC, AG, AT, CA, etc.) in the input nucleotide sequence. 
#' For each sequence, a feature vector is created where each element represents the 
#' frequency of a specific dinucleotide. Additionally, a correlation factor is calculated 
#' across tiers (up to 6) to capture dependencies between physicochemical properties 
#' of dinucleotides separated by various distances within the sequence. 
#' These physicochemical properties include Enthalpy, Entropy, and Free energy.
#' The resulting matrix is normalized by row sums to provide a consistent representation 
#' of the sequence characteristics.
#' }
#' }
#'
#' @param  cnn Whether report in CNN format; default is \code{FALSE}.
#'
#' @details The function first extract DNA sequence within the genomic ranges defined by \code{x}.
#' Then, the DNA strings are processed according to the selected sequence encoding method.
#'
#' @return If \code{cnn == FALSE}, the returned object is a \code{data.frame} whose number of rows is the length of \code{x}, and the number of columns is 4 times the width of \code{x}.
#' The column types in the \code{data.frame} are all numeric.
#'
#' When \code{cnn == TRUE}, the returned object is an \code{array} of n matrices corresponding to the feature vectors of each position in the sequence.
#'
#' If \code{mutation} is provided, a list of 2 matrices/arrays of equal sizes will be returned, while the first is for the normal sequences, and
#' the 2nd is for the mutated sequences. The mapping between the matrix and the rows of mutation GRanges is labeled in the rownames.
#'
#' @examples
#' library(BSgenome.Hsapiens.UCSC.hg19)
#'
#' ## Define the Granges to be annotated:
#' set.seed(01)
#'
#' X <- GRanges(rep(c("chr1", "chr2"), c(15, 15)),
#'              IRanges(c(sample(11874:12127, 15), sample(38814:41527,15)), width=5),
#'              strand=Rle(c("+", "-"), c(15, 15)))
#'
#' ## Extract onehot encoded sequence features
#' seq_onehot <- sequenceDerivedFeatures(X, Hsapiens, encoding = "onehot")
#' str(seq_onehot)
#'
#' ## Extract iRNA encoded sequence features
#' seq_iRNA <- sequenceDerivedFeatures(X, Hsapiens, encoding = "iRNA")
#' str(seq_iRNA)
#'
#' ## Extract Pseudo dinucleotide composition (PseDNC) encoded sequence features
#' seq_PseDNC <- sequenceDerivedFeatures(X, Hsapiens, encoding = "PseDNC")
#' str(seq_PseDNC)
#'
#' ## Extract mutated features
#' SNP <- resize(X[c(10,15,3,6)], 1)
#' SNP$mutateTo <- c("A","C","G","C")
#' seq_mutation <- sequenceDerivedFeatures(X, Hsapiens, mutation = SNP)
#'
#' rownames(seq_mutation[[1]]) #index for x
#' rownames(seq_mutation[[2]]) #index for mutation
#'
#' @seealso
#'
#' \itemize{
#' \item{}{The \link{genomeDerivedFeatures} for extraction of genome-derived features.}
#' \item{}{The \link{topologyOnTranscripts} for calculation of the meta-tx topologies on transcripts of genes.}
#' }
#'
#' @importFrom matrixStats rowCumsums
#' @importFrom IRanges Views
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings complement
#' @export
sequenceDerivedFeatures <- function(x,
                                    sequence,
                                    mutation = NULL,
                                    encoding = c("onehot", "iRNA", "PseDNC"),
                                    cnn = FALSE) {
  encoding = match.arg(encoding)
  stopifnot(all(width(x) == width(x)[1]))
  
  if (!is.null(mutation)) {
    stopifnot(is(mutation, "GRanges"))
    stopifnot(all(width(mutation) == 1))
    stopifnot(all(mcols(mutation)[[1]] %in% c("A", "C", "G", "T")))
  }
  
  ##Fill unknown attributes of x with the attributes of BSgenome
  if (anyNA(isCircular(x))) {
    isCircularGenome <- isCircular(sequence)
    isCircular(x) <- isCircularGenome[names(isCircular(x))]
    rm(isCircularGenome)
  }
  if (anyNA(seqlengths(x))) {
    seqLenGenome <- seqlengths(sequence)
    seqlengths(x) <- seqLenGenome[names(seqlengths(x))]
    rm(seqLenGenome)
  }
  
  N = width(x)[1]
  
  if (is.null(mutation)) {
    sequences <- as.character(DNAStringSet(Views(sequence, x)))
    sequence_M <-
      matrix(unlist(strsplit(sequences, "")), ncol =  N, byrow = TRUE)
    rm(sequences)
    
    seqfeatures <- assign_method(sequence_M, encoding)
    
    if (cnn) {
      seqfeatures <- wrap_cnn(seqfeatures, N)
    } else{
      seqfeatures <- as.data.frame(seqfeatures)
    }
    return(seqfeatures)
    
  } else{
    overlap <-
      suppressWarnings(findOverlaps(mutation, x, ignore.strand = TRUE))
    QrHit <- queryHits(overlap)
    SbHit <- subjectHits(overlap)
    rm(overlap)
    
    if (length(QrHit) == 0) {
      stop("No mutations can be mapped to x.")
    } else{
      mutHit <- mutation[QrHit]
      xHit <- x[SbHit]
      sequences_reference <- DNAStringSet(Views(sequence, xHit))
      pos_replace_plus <- start(mutHit) - start(xHit) + 1
      pos_replace_minus <- end(xHit) - start(mutHit) + 1
      nt_replace_plus <- as.character(mcols(mutHit)[[1]])
      nt_replace_minus <- as.character(complement(DNAStringSet(mcols(mutHit)[[1]])))
      indx_minus_range <- as.logical(strand(xHit) == "-")
      rm(xHit, mutHit)
      pos_plus_indx <-
        split(pos_replace_plus, seq_along(pos_replace_plus))[!indx_minus_range]
      names(pos_plus_indx) <- NULL
      sequences_mutated <- sequences_reference
      sequences_mutated[!indx_minus_range][pos_plus_indx] <-
        DNAStringSet(nt_replace_plus[!indx_minus_range])
      pos_minus_indx <-
        split(pos_replace_minus, seq_along(pos_replace_minus))[indx_minus_range]
      names(pos_minus_indx) <- NULL
      sequences_mutated[indx_minus_range][pos_minus_indx] <-
        DNAStringSet(nt_replace_minus[indx_minus_range])
      rm(
        pos_replace_plus,
        pos_replace_minus,
        indx_minus_range,
        pos_plus_indx,
        pos_minus_indx
      )
      sequence_M_ref <-
        matrix(unlist(strsplit(
          as.character(sequences_reference), ""
        )), ncol =  N, byrow = TRUE)
      seqfeatures_ref <- assign_method(sequence_M_ref, encoding)
      rm(sequence_M_ref, sequences_reference)
      sequence_M_mut <-
        matrix(unlist(strsplit(
          as.character(sequences_mutated), ""
        )), ncol =  N, byrow = TRUE)
      seqfeatures_mut <- assign_method(sequence_M_mut, encoding)
      rm(sequence_M_mut, sequences_mutated)
      feature_lst <- list(ref = seqfeatures_ref,
                          mut = seqfeatures_mut)
      if (cnn)
        feature_lst <- lapply(feature_lst, function(x)
          wrap_cnn(x, N))
      rownames(feature_lst[[1]]) <- SbHit
      rownames(feature_lst[[2]]) <- QrHit
      return(feature_lst)
    }
  }
}

onehot_encode <- function(sequence_M) {
  A_M <- sequence_M == "A"
  colnames(A_M) <- paste0("A_", seq_len(ncol(A_M)))
  T_M <- sequence_M == "T"
  colnames(T_M) <- paste0("T_", seq_len(ncol(T_M)))
  C_M <- sequence_M == "C"
  colnames(C_M) <- paste0("C_", seq_len(ncol(C_M)))
  G_M <- sequence_M == "G"
  colnames(G_M) <- paste0("G_", seq_len(ncol(G_M)))
  N <- ncol(sequence_M)
  rm(sequence_M)
  onehotfeatures <- cbind(A_M, T_M, C_M, G_M)
  rm(A_M, T_M, C_M, G_M)
  onehotfeatures <-
    as.data.frame(onehotfeatures[, rep(seq_len(N), each = 4) + rep(c(0, N, 2 *
                                                                       N, 3 * N), N)])
  onehotfeatures <- apply(onehotfeatures, 2, as.numeric)
  return(onehotfeatures)
}

iRNA_encode <- function(sequence_M) {
  N <- ncol(sequence_M)
  purine_M <- sequence_M == "A" | sequence_M == "G"
  colnames(purine_M) <-  paste0("purine_", seq_len(N))
  amino_M <- sequence_M == "A" | sequence_M == "C"
  colnames(amino_M) <-  paste0("aminoGroup_", seq_len(N))
  weakHyb_M <- sequence_M == "A" | sequence_M == "T"
  colnames(weakHyb_M) <- paste0("weakHydrogenBonds_", seq_len(N))
  
  cumFreq_A <-
    rowCumsums(matrix(
      as.numeric(sequence_M == "A"),
      ncol = N,
      byrow = FALSE
    ))
  cumFreq_T <-
    rowCumsums(matrix(
      as.numeric(sequence_M == "T"),
      ncol = N,
      byrow = FALSE
    ))
  cumFreq_C <-
    rowCumsums(matrix(
      as.numeric(sequence_M == "C"),
      ncol = N,
      byrow = FALSE
    ))
  cumFreq_G <-
    rowCumsums(matrix(
      as.numeric(sequence_M == "G"),
      ncol = N,
      byrow = FALSE
    ))
  
  cumFreq_combined <- matrix(0, ncol = N, nrow = nrow(sequence_M))
  cumFreq_combined[sequence_M == "A"] <-
    cumFreq_A[sequence_M == "A"]
  cumFreq_combined[sequence_M == "T"] <-
    cumFreq_T[sequence_M == "T"]
  cumFreq_combined[sequence_M == "C"] <-
    cumFreq_C[sequence_M == "C"]
  cumFreq_combined[sequence_M == "G"] <-
    cumFreq_G[sequence_M == "G"]
  
  cumFreq_combined <- t(t(cumFreq_combined) / seq_len(N))
  colnames(cumFreq_combined) <-
    paste0("cumulativeFrequency_", seq_len(N))
  seqfeatures <-
    as.data.frame(cbind(
      cbind(purine_M, amino_M, weakHyb_M),
      as.data.frame(cumFreq_combined)
    ))
  seqfeatures <- apply(seqfeatures, 2, as.numeric)
  return(seqfeatures)
}

wrap_cnn <- function(seqfeatures, N) {
  num_features <- ncol(seqfeatures) / N
  cnn_array <-
    array(as.vector(seqfeatures), dim = c(nrow(seqfeatures), num_features, N))
  cnn_array <- aperm(cnn_array, c(1, 3, 2))
  return(cnn_array)
}

assign_method <- function(sequence_M, encoding) {
  if (encoding == "onehot") {
    seqfeatures <- onehot_encode(sequence_M)
  } else if (encoding == "iRNA") {
    seqfeatures <- iRNA_encode(sequence_M)
  } else if (encoding == "PseDNC") {
    seqfeatures <- PseDNC_encode(sequence_M)
  }
}

get_dinucleotide_physicochemical_properties <- function() {
  # # Define possible dinucleotides
  # dinucleotides <- c("GG", "GA", "GC", "GT", "AG", "AA", "AC", "AT", 
  #                    "CG", "CA", "CC", "CT", "TG", "TA", "TC", "TT")

  # # Define physicochemical properties matrix
  # properties <- matrix(c(
  #   -12.2, -29.7, -3.26,
  #   -13.3, -35.5, -2.35,
  #   -14.2, -34.9, -3.42,
  #   -10.2, -26.2, -2.24,
  #   -7.6, -19.2, -2.08,
  #   -6.6, -18.4, -0.93,
  #   -10.2, -26.2, -2.24,
  #   -5.7, -15.5, -1.10,
  #   -8.0, -19.4, -2.36,
  #   -10.5, -27.8, -2.11,
  #   -12.2, -29.7, -3.26,
  #   -7.6, -19.2, -2.08,
  #   -7.6, -19.2, -2.11,
  #   -8.1, -22.6, -1.33,
  #   -10.2, -26.2, -2.35,
  #   -6.6, -18.4, -0.93
  # ), ncol = 3, byrow = TRUE)
  # colnames(properties) <- c("Enthalpy", "Entropy", "FreeEnergy")
  # rownames(properties) <- dinucleotides

  # # normalise
  # properties = scale(properties)

  # dput(properties)
  properties = structure(c(-1.07588314958376, -1.50235935302237, -1.85129442856305, 
-0.300471870604474, 0.707562792068599, 1.09526843155824, -0.300471870604474, 
1.44420350709892, 0.552480536272742, -0.416783562451367, -1.07588314958376, 
0.707562792068599, 0.707562792068599, 0.513709972323777, -0.300471870604474, 
1.09526843155824, -0.882196127443816, -1.82212495209119, -1.72489093574836, 
-0.314997698777298, 0.819399158555737, 0.949044513679512, -0.314997698777298, 
1.4190089260032, 0.786987819774793, -0.57428840902485, -0.882196127443816, 
0.819399158555737, 0.819399158555737, 0.268406399279691, -0.314997698777298, 
0.949044513679512, -1.45493929364262, -0.278708526544534, -1.66174909884668, 
-0.136526785466743, 0.0702830197373169, 1.55672849464149, -0.136526785466743, 
1.33699307661218, -0.291634139369787, 0.031506181261556, -1.45493929364262, 
0.0702830197373169, 0.031506181261556, 1.03970398163135, -0.278708526544534, 
1.55672849464149), dim = c(16L, 3L), dimnames = list(c("GG", 
"GA", "GC", "GT", "AG", "AA", "AC", "AT", "CG", "CA", "CC", "CT", 
"TG", "TA", "TC", "TT"), c("Enthalpy", "Entropy", "FreeEnergy"
)))

  return(properties)
}

# Pseudo dinucleotide composition (PseDNC)
# https://doi.org/10.1016/j.ab.2015.08.021
PseDNC_encode <- function(sequence_M) {
  # Define possible dinucleotides
  properties <- get_dinucleotide_physicochemical_properties()
  dinucleotides <- rownames(properties)

  # Convert the sequence matrix to a dinucleotide index matrix
  di_sequence_M <- paste0(sequence_M[, -ncol(sequence_M)], sequence_M[, -1])
  di_sequence_M <- matrix(match(di_sequence_M, dinucleotides), nrow = nrow(sequence_M))

  # Initialize a matrix to store dinucleotide frequencies
  result_matrix <- matrix(0, nrow = nrow(sequence_M), ncol = length(dinucleotides))

  # Calculate dinucleotide frequencies for each of the 16 dinucleotides
  for (di_index in seq_along(dinucleotides)) {
    result_matrix[,di_index] <- rowMeans(di_sequence_M == di_index, na.rm = T)
  }
  colnames(result_matrix) <- dinucleotides

  # Calculate correlation factor for tiers 1 to (maximum) 6 
  tier_N = min(6, nrow(sequence_M) - 1)
  result_correlation = matrix(NA, nrow = nrow(sequence_M), ncol = tier_N)
  colnames(result_correlation) = paste0("tier_", seq_len(tier_N))
  for (tier_index in seq_len(tier_N)) {
    tier_correlation = rep(0, nrow(sequence_M))
    for (property_index in seq_len(ncol(properties))) {
      property_di_sequence_M = matrix(properties[di_sequence_M, property_index], nrow = nrow(di_sequence_M))
      tem <- property_di_sequence_M[, 1:(ncol(property_di_sequence_M) - tier_index)] - 
             property_di_sequence_M[, (1 + tier_index):ncol(property_di_sequence_M)]
      tier_correlation <- tier_correlation + rowMeans(tem^2, na.rm = TRUE)
    }

    # Average over properties
    result_correlation[, tier_index] <- tier_correlation / ncol(properties)
  }

  # Weight and normalize
  result_matrix = cbind(result_matrix, result_correlation * 0.9)
  result_matrix = result_matrix / rowSums(result_matrix, na.rm = TRUE)

  return(result_matrix)
}

                              
## To Do: add more sequence feature encoding methods.
