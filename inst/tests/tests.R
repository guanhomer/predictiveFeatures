library(GenomicRanges)
library(testthat)
library(predictiveFeatures)

x_gr <- GRanges(rep(c("chr1", "chr2"), c(5, 15)),
                IRanges(c(sample(11874:12127, 5), sample(38814:41527,15)), width=100),
                strand=Rle(c("+", "-"), c(5, 15)))

exons_gr <- GRanges(c("chr1","chr2","chr2"),
                    IRanges(start=c(11874,38814,45440),end=c(12227,41627,46588)),
                    strand=c("+","-","-"))
genes_grl <- GRangesList(gene1=exons_gr[1],gene2=exons_gr[c(2,3)])

y_gr <- GRanges(rep(c("chr1", "chr2"), c(50, 50)),
                IRanges(c(sample(11874:12127, 50),
                          sample(38814:41527,50)), width=1),
                strand=Rle(c("+", "-"), c(50, 50)))

expect_that(extractRegionLength(x_gr), is_a("integer"))
expect_that(extractRegionLength(x_gr, exons_gr), is_a("numeric"))
expect_that(extractRegionLength(x_gr, exons_gr), is_a("numeric"))
exons_property <- c(1,6,8)
expect_that(extractRegionProperty(x_gr, exons_gr, exons_property), is_a("numeric"))
expect_that(extractRegionYCount(x_gr, x_gr, exons_gr), is_a("numeric"))
expect_that(extractRegionNearestDistToY(x_gr, y_gr), is_a("numeric"))
expect_that(extractDistToRegion5end(x_gr, exons_gr), is_a("numeric"))
expect_that(extractDistToRegion3end(x_gr, exons_gr), is_a("numeric"))
