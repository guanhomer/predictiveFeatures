#' predictiveFeatures: extract genomic features for predictive modeling in genomics
#'
#' The package can be used to extract properties of region lengths, sequence contents, relative positions, clustering effects, and conservation scores on top of the classically defined gene annotations including 5'UTR, CDS, 3'UTR, transcripts, genes, and promoters. The genome-derived features can be used to predict variety of gene related genomic assays, such as the RBP binding sites and the RNA modification sites. predictiveFeatures also provide the annotation of new attributes over user-defined genomic regions, enabling the augmentation of feature engineering under the interactive framework between genomic properties and genomic regions.
#' 
#' The package currently implements 2 types of feature extraction modules: the \strong{genome-derived features} and the \strong{sequence-derived features}.The former is encoded via the interaction between genomic properties and the genomic regions, and the later is defined through the nucleotide encoding methods, such as the one-hot encoding or the encoding of pseudo nucleotide compositions (PseTNC). The genome and sequence features can be extracted with the function \link{genomeDerivedFeatures} and \link{sequenceDerivedFeatures} respectively. Furthermore, more genome-derived features can be defined using functions with names "extractRegion" + the metric names. While the individual extractor functions can annotate new properties on the user defined region types.
#' 
#' The design of package enables the one-step extraction of hundreds of genomic features based only on the \code{GRanges} of the target genome intervals. Wide class of the Bioconductor annotation objects, including \code{Txdb}, \code{EnsDb}, \code{BSgenome}, and \code{GScores}, are supported in the package to enable the massive enumeration of features that based on novel sources of genome regions.
#' 
#' While the comprehensive feature input is vital to the success of a genomic data science project, the genome-derived features are often more informative and interpretive than the primary sequence features. Thus, the \code{predictiveFearures} package can provide general support for machine learning modeling in the domain of genomics.
#' 
#' The user's guide will firstly give an overview to the key functionalities of the package, then it will illustrate the utility of the predictive features using a case study of the full workflow in a genomic prediction task, from EDA to modeling.
#'
#' @docType package
#' @name predictiveFeatures
#' @section Key predictiveFeatures functions:
#' \link{genomeDerivedFeatures}
#' 
#' \link{sequenceDerivedFeatures}
#' 
#' \link{topologyOnTranscripts}
NULL
#> NULL