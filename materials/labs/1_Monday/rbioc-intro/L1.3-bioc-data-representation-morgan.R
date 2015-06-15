## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()

## ----setup, echo=FALSE---------------------------------------------------
knitr::opts_chunk$set(cache=TRUE)
suppressPackageStartupMessages({
    library(GenomicRanges)
})

## ----iranges-------------------------------------------------------------
library(IRanges)
ir <- IRanges(start=c(10, 20, 30), width=5)
ir

## ----iranges-flank-------------------------------------------------------
flank(ir, 3)

## ----iranges-class-------------------------------------------------------
class(ir)
getClass(class(ir))

## ----iranges-flank-method, eval=FALSE------------------------------------
## ?"flank,Ranges-method"

## ----granges-------------------------------------------------------------
library(GenomicRanges)
gr <- GRanges(c("chr1", "chr1", "chr2"), ir, strand=c("+", "-", "+"))
gr

## ----granges-flank-------------------------------------------------------
flank(gr, 3)

## ----granges-class-------------------------------------------------------
class(gr)
getClass(class(gr))

## ----granges-flank-method, eval=FALSE------------------------------------
## ?"flank,GenomicRanges-method"

## ----granges-methods, eval=FALSE-----------------------------------------
## showMethods(class="GRanges", where=search())

## ----granges-man-and-vignettes, eval=FALSE-------------------------------
## help(package="GenomicRanges")
## vignette(package="GenomicRanges")
## vignette(package="GenomicRanges", "GenomicRangesHOWTOs")

## ----ranges, message=FALSE-----------------------------------------------
library(GenomicRanges)
gr <- GRanges("A", IRanges(c(10, 20, 22), width=5), "+")
shift(gr, 1)                            # 1-based coordinates!
range(gr)                               # intra-range
reduce(gr)                              # inter-range
coverage(gr)
setdiff(range(gr), gr)                  # 'introns'

## ----BSgenome-require, message=FALSE-------------------------------------
library(BSgenome.Hsapiens.UCSC.hg19)
chr14_range = GRanges("chr14", IRanges(1, seqlengths(Hsapiens)["chr14"]))
chr14_dna <- getSeq(Hsapiens, chr14_range)
letterFrequency(chr14_dna, "GC", as.prob=TRUE)

## ----bam-require---------------------------------------------------------
library(GenomicRanges)
library(GenomicAlignments)
library(Rsamtools)

## our 'region of interest'
roi <- GRanges("chr14", IRanges(19653773, width=1)) 
## sample data
library('RNAseqData.HNRNPC.bam.chr14')
bf <- BamFile(RNAseqData.HNRNPC.bam.chr14_BAMFILES[[1]], asMates=TRUE)
## alignments, junctions, overlapping our roi
paln <- readGAlignmentsList(bf)
j <- summarizeJunctions(paln, with.revmap=TRUE)
j_overlap <- j[j %over% roi]

## supporting reads
paln[j_overlap$revmap[[1]]]

## ----vcf, message=FALSE--------------------------------------------------
## input variants
library(VariantAnnotation)
fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
vcf <- readVcf(fl, "hg19")
seqlevels(vcf) <- "chr22"
## known gene model
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
coding <- locateVariants(rowRanges(vcf),
    TxDb.Hsapiens.UCSC.hg19.knownGene,
    CodingVariants())
head(coding)

## ----summarizeOverlaps-roi, message=FALSE--------------------------------
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
exByGn <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")
## only chromosome 14
seqlevels(exByGn, force=TRUE) = "chr14"

## ----summarizeOverlaps-bam, message=FALSE--------------------------------
library(RNAseqData.HNRNPC.bam.chr14)
length(RNAseqData.HNRNPC.bam.chr14_BAMFILES)

## ----summarizeOverlaps---------------------------------------------------
## next 2 lines optional; non-Windows
library(BiocParallel)
register(MulticoreParam(workers=detectCores()))
olaps <- summarizeOverlaps(exByGn, RNAseqData.HNRNPC.bam.chr14_BAMFILES)

## ----summarizeOverlaps-explore-------------------------------------------
olaps
head(assay(olaps))
colSums(assay(olaps))                # library sizes
plot(sum(width(olaps)), rowMeans(assay(olaps)), log="xy")

## ----summarizeOverlaps-gc------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg19)
sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg19, rowRanges(olaps))
gcPerExon <- letterFrequency(unlist(sequences), "GC")
gc <- relist(as.vector(gcPerExon), sequences)
gc_percent <- sum(gc) / sum(width(olaps))
plot(gc_percent, rowMeans(assay(olaps)), log="y")

