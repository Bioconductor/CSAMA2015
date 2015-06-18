## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()

## ----setup, echo=FALSE---------------------------------------------------
##knitr::opts_chunk$set(cache=TRUE)

## ----vectorize-----------------------------------------------------------
x <- 1:10
log(x)     ## NOT for (i in seq_along) x[i] <- log(x[i])

## ----pre-allocate--------------------------------------------------------
result <- numeric(10)
result[1] <- runif(1)
for (i in 2:length(result))
       result[i] <- runif(1) * result[i - 1]
result

## ----inefficient---------------------------------------------------------
f0 <- function(n, a=2) {
    ## stopifnot(is.integer(n) && (length(n) == 1) &&
    ##           !is.na(n) && (n > 0))
    result <- numeric()
    for (i in seq_len(n))
        result[[i]] <- a * log(i)
    result
}

## ----system-time---------------------------------------------------------
system.time(f0(10000))
n <- 1000 * seq(1, 20, 2)
t <- sapply(n, function(i) system.time(f0(i))[[3]])
plot(t ~ n, type="b")

## ----correct-init--------------------------------------------------------
n <- 10000
system.time(expected <- f0(n))
head(expected)

## ----hoist---------------------------------------------------------------
f1 <- function(n, a=2) {
    result <- numeric()
    for (i in seq_len(n))
        result[[i]] <- log(i)
    a * result
}
identical(expected, f1(n))

library(microbenchmark)
microbenchmark(f0(n), f1(n), times=5)

## ----preallocate-and-fill------------------------------------------------
f2 <- function(n, a=2) {
    result <- numeric(n)
    for (i in seq_len(n))
        result[[i]] <- log(i)
    a * result
}
identical(expected, f2(n))
microbenchmark(f0(n), f2(n), times=5)

## ----use-apply-----------------------------------------------------------
f3 <- function(n, a=2)
    a * sapply(seq_len(n), log)

identical(expected, f3(n))
microbenchmark(f0(n), f2(n), f3(n), times=10)

## ----use-vectorize-------------------------------------------------------
f4 <- function(n, a=2)
    a * log(seq_len(n))
identical(expected, f4(n))
microbenchmark(f0(n), f3(n), f4(n), times=10)

## ----use-compiler--------------------------------------------------------
library(compiler)
f2c <- cmpfun(f2)
n <- 10000
identical(f2(n), f2c(n))
microbenchmark(f2(n), f2c(n), times=10)

## ----vectorized-scale----------------------------------------------------
n <- 10 ^ (5:8)                         # 100x larger than f0
t <- sapply(n, function(i) system.time(f4(i))[[3]])
plot(t ~ n, log="xy", type="b")

## ----algo-init-----------------------------------------------------------
vec <- c(seq(-100,-1,length.out=1e6), rep(0,20), seq(1,100,length.out=1e6))

## ----algo-scan-----------------------------------------------------------
f0 <- function(v) sum(v < 0)

## ----algo-scan-time------------------------------------------------------
N <- seq(1, 11, 2) * 1e6
Time <- sapply(N, function(n) {
    v <- sort(rnorm(n))
    system.time(f0(v))[[3]]
})
plot(Time ~ N, type="b")

## ----algo-binary---------------------------------------------------------
f3 <- function(x, threshold=0) {
    imin <- 1L
    imax <- length(x)
    while (imax >= imin) {
        imid <- as.integer(imin + (imax - imin) / 2)
        if (x[imid] >= threshold)
            imax <- imid - 1L
        else
            imin <- imid + 1L
    }
    imax
}

## ----algo-binary-perform-------------------------------------------------
## identity
stopifnot(
    identical(f0((-2):2), f3((-2):2)),
    identical(f0(2:4), f3(2:4)),
    identical(f0(-(4:2)), f3(-(4:2))),
    identical(f0(vec), f3(vec)))

## scale
N <- 10^(1:7)

Time <- sapply(N, function(n) {
    v <- sort(rnorm(n))
    system.time(f3(v))[[3]]
})
plot(Time ~ N, type="b")

## ----algo-relative-------------------------------------------------------
## relative time
library(microbenchmark)
microbenchmark(f0(vec), f3(vec))

library(compiler)
f3c <- cmpfun(f3)
microbenchmark(f3(vec), f3c(vec))

## ----findInterval--------------------------------------------------------
f4 <- function(v, query=0)
    findInterval(query - .Machine$double.eps, v)

identical(f0(vec), f4(vec))
microbenchmark(f0(vec), f3(vec), f4(vec))

## ----findInterval-several------------------------------------------------
threshold <- rnorm(10000)
identical(sapply(threshold, f3, x=vec), f4(vec, threshold))
microbenchmark(sapply(x, f3), f4(vec, x))

## ----parallel-setup, echo=FALSE------------------------------------------
suppressPackageStartupMessages({
    library(RNAseqData.HNRNPC.bam.chr14)
    library(Biostrings)
    library(Rsamtools)
    library(GenomicFiles)
    library(GenomicAlignments)
    library(BiocParallel)
    library(ggplot2)
})

## ----bam-files-----------------------------------------------------------
library(RNAseqData.HNRNPC.bam.chr14)
fls <- RNAseqData.HNRNPC.bam.chr14_BAMFILES

## ----BamFileList---------------------------------------------------------
library(Rsamtools)
bfls <- BamFileList(fls, yieldSize=100000)

## ----yield---------------------------------------------------------------
library(GenomicAlignments)
yield <- function(bfl) # input a chunk of alignments
    readGAlignments(bfl, param=ScanBamParam(what="seq"))

## ----map-----------------------------------------------------------------
library(Biostrings)
map <- function(aln) { # GC content, bin & cummulate
    gc <- letterFrequency(mcols(aln)$seq, "GC")
    tabulate(1 + gc, 73)                # max. read length: 72
}

## ----reduce--------------------------------------------------------------
reduce <- `+`

## ----reduceByYield-------------------------------------------------------
library(GenomicFiles)
bf <- BamFile(fls[1], yieldSize=100000)
reduceByYield(bf, yield, map, reduce)

## ----bplapply------------------------------------------------------------
library(BiocParallel)
gc <- bplapply(bfls, reduceByYield, yield, map, reduce)

## ----lst2data.frame------------------------------------------------------
library(ggplot2)
df <- stack(as.data.frame(lapply(gc, cumsum)))
df$GC <- 0:72

## ----gc-plot-------------------------------------------------------------
library(ggplot2)
ggplot(df, aes(x=GC, y=values)) + geom_line(aes(colour=ind)) +
    xlab("Number of GC Nucleotides per Read") +
    ylab("Number of Reads")

