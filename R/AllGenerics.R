#' Method extractRegionProperty
#' 
#' @param x A \link{GRanges} object for the query.
#' @param region A \link{GRanges} or \link{GRangesList} object for the subject. 
#' If not provided (\code{region=NULL}), the returned properties will be calculated directly from \code{x}.
#' @param property A vector specifying the properties/attributes of the \code{region}. 
#' 
#' @param ambiguityMethod By default, for \code{logical} \code{property} input, if x overlaps multiple regions, 
#' as long as any of the properties is TRUE, the returned value will be TRUE. 
#' For other types of \code{property} input, the returned value will be the average of the properties over multiple overlapped regions.
#' If \code{ambiguityMethod} is \code{"mean"}, \code{"sum"}, \code{"min"}, or \code{"max"},
#' then the mean, sum, minimum, and maximum values of the >1 mapping will be returned in the output value.
#' 
#' @param maxgap,minoverlap,type See \code{?\link{findOverlaps}} in the \bold{IRanges} package for a description of these arguments. 
#' 
#' @param nomapValue When \code{nomapValue} is \code{"NA"}, \code{"0"}, or \code{"FALSE"}, 
#' the \code{x} that do not match the \code{region} will return \code{NA}, \code{0}, and \code{FALSE} respectively.
#' If \code{nomapValue} is \code{"nearest"}, the not matched \code{x} will be set to be the property value on its nearest \code{region}.
#'
#' @param sequence A \link{BSgenome} or \link{XStringSet} object. 
#' See the \link{available.genomes} function for how to install a genome.
#' 
#' @param gscores A \link{GScores} object. 
#' See the vignette for more details: \code{vignette("GenomicScores", package = "GenomicScores")}.
#' 
#' @param as.prob,letters See \code{?\link{letterFrequency}} in the \bold{Biostrings} package for a description of these arguments.
#'
#' @param ignore.strand When set to \code{TRUE}, the strand information is ignored in the overlap calculations.
#' 
#' @param output.logical If \code{TRUE} then the returned value is \code{logical}, otherwise \code{numeric}.
#' 
#' @param efficient \code{TRUE} if only internally extract the properties on \code{regions} overlap with \code{x}, 
#' the option makes the computation more efficient if the number of gnomic regions is much larger than the query.
#' 
#' @param missingScores Specifying the imputation methods on missing scores in \code{extractRegionScores}, 
#' default is returning zero for unknown gscores under the region.
#' 
#' @param y A \link{GRanges} or \link{GRangesList} object for the clustering annotation.
#' 
#' @param normalize If \code{TRUE}, then returned count in \code{extractRegionYCount} will be divided by the length of the region.
#' The normalized count can be understood as the density of the annotation on the area.
#' 
#' @param maxDist A value of the maximum nucleotide distance returned in \code{extractRegionNearestDistToY}, default: \code{3e+06}.
#' 
#' @param ... Additional arguments, it passes to \link{letterFrequency} in method \code{extractRegionLetterFrequency}, and it passes to \link{score} in method \code{extractRegionScores}.
#' 
#' @rdname extractRegionProperty-methods
#' @details Please check the documentation of \link{extractRegionProperty} for details.
#' @exportMethod extractRegionProperty
setGeneric("extractRegionProperty", function(x, 
                                             region, 
                                             property,
                                             ambiguityMethod = c("auto", "mean", "sum", "min", "max"),
                                             maxgap=0L, 
                                             minoverlap=0L,
                                             type=c("any", "start", "end", "within", "equal"),
                                             nomapValue=c("NA","0","FALSE","nearest"),
                                             ignore.strand=FALSE){standardGeneric("extractRegionProperty")})

#' Method extractRegionOverlap
#' @rdname extractRegionProperty-methods
#' @exportMethod extractRegionOverlap
setGeneric("extractRegionOverlap", function(x,
                                            region=NULL,
                                            ambiguityMethod=c("auto", "mean", "sum", "min", "max"),
                                            maxgap=-1L, 
                                            minoverlap=0L,
                                            type=c("any", "start", "end", "within", "equal"),
                                            ignore.strand=FALSE,
                                            output.logical=TRUE){standardGeneric("extractRegionOverlap")})

#' Method extractRegionLength
#' @rdname extractRegionProperty-methods
#' @exportMethod extractRegionLength
setGeneric("extractRegionLength", function(x,
                                           region=NULL,
                                           ambiguityMethod=c("mean", "sum", "min", "max"),
                                           maxgap=-1L, 
                                           minoverlap=0L,
                                           type=c("any", "start", "end", "within", "equal"),
                                           nomapValue=c("NA","0","nearest"),
                                           ignore.strand=FALSE){standardGeneric("extractRegionLength")})

#' Method extractRegionLetterFrequency
#' @rdname extractRegionProperty-methods
#' @exportMethod extractRegionLetterFrequency
setGeneric("extractRegionLetterFrequency", function(x,
                                                    sequence,
                                                    region=NULL,
                                                    letters = "GC",
                                                    as.prob = TRUE,
                                                    ambiguityMethod=c("mean", "sum", "min", "max"),
                                                    maxgap=-1L, 
                                                    minoverlap=0L,
                                                    type=c("any", "start", "end", "within", "equal"),
                                                    nomapValue=c("NA","0","nearest"),
                                                    ignore.strand=FALSE,
                                                    efficient=TRUE,
                                                    ...){standardGeneric("extractRegionLetterFrequency")})

#' Method extractRegionScores
#' @rdname extractRegionProperty-methods
#' @exportMethod extractRegionScores
setGeneric("extractRegionScores", function(x,
                                           gscores,
                                           region=NULL,
                                           ambiguityMethod=c("mean", "sum", "min", "max"),
                                           maxgap=-1L, 
                                           minoverlap=0L,
                                           type=c("any", "start", "end", "within", "equal"),
                                           nomapValue=c("NA","0","nearest"),
                                           missingScores=c("zero", "mean", "none"),
                                           ignore.strand=FALSE,
                                           efficient=TRUE,
                                           ...){standardGeneric("extractRegionScores")})


#' Method extractRegionYCount
#' @rdname extractRegionProperty-methods
#' @exportMethod extractRegionYCount
setGeneric("extractRegionYCount", function(x,
                                           y,
                                           region=NULL,
                                           normalize=FALSE,
                                           ambiguityMethod=c("mean", "sum", "min", "max"),
                                           maxgap=-1L, 
                                           minoverlap=0L,
                                           type=c("any", "start", "end", "within", "equal"),
                                           nomapValue=c("NA","0","nearest"),
                                           ignore.strand=FALSE,
                                           efficient=TRUE){standardGeneric("extractRegionYCount")})


#' Method extractRegionNearestDistToY
#' @rdname extractRegionProperty-methods
#' @exportMethod extractRegionNearestDistToY
setGeneric("extractRegionNearestDistToY", function(x,
                                                   y,
                                                   region=NULL,
                                                   ambiguityMethod=c("mean", "sum", "min", "max"),
                                                   maxgap=-1L, 
                                                   minoverlap=0L,
                                                   type=c("any", "start", "end", "within", "equal"),
                                                   nomapValue=c("NA","0","nearest"),
                                                   maxDist=3e6,
                                                   ignore.strand=FALSE){standardGeneric("extractRegionNearestDistToY")})

#' Method extractRegionRelativePosition
#' @rdname extractRegionProperty-methods
#' @exportMethod extractRegionRelativePosition
setGeneric("extractRegionRelativePosition", function(x,
                                                     region=NULL,
                                                     ambiguityMethod=c("mean", "sum", "min", "max"),
                                                     nomapValue=c("NA","0"),
                                                     ignore.strand=FALSE){standardGeneric("extractRegionRelativePosition")})


#' Method extractDistToRegion5end
#' @rdname extractRegionProperty-methods
#' @exportMethod extractDistToRegion5end
setGeneric("extractDistToRegion5end", function(x,
                                               region=NULL,
                                               ignore.strand=FALSE,
                                               ambiguityMethod=c("mean", "sum", "min", "max"),
                                               maxgap=-1L, 
                                               minoverlap=0L,
                                               type=c("any", "start", "end", "within", "equal"),
                                               nomapValue=c("NA","0","nearest")){standardGeneric("extractDistToRegion5end")})

#' Method extractDistToRegion3end
#' @rdname extractRegionProperty-methods
#' @exportMethod extractDistToRegion3end
setGeneric("extractDistToRegion3end", function(x,
                                               region=NULL,
                                               ignore.strand=FALSE,
                                               ambiguityMethod=c("mean", "sum", "min", "max"),
                                               maxgap=-1L, 
                                               minoverlap=0L,
                                               type=c("any", "start", "end", "within", "equal"),
                                               nomapValue=c("NA","0","nearest")){standardGeneric("extractDistToRegion3end")})
