## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()
options(width=100, max.print=1000)

## ----setup, echo=FALSE, messages=FALSE, warnings=FALSE-------------------
knitr::opts_chunk$set(cache=TRUE)
suppressPackageStartupMessages({
    library(Biostrings)
    library(GenomicRanges)
})

## ----setup-objects-------------------------------------------------------
library(Biostrings)
library(GenomicRanges)

## ----Biostrings, message=FALSE-------------------------------------------
library(Biostrings)                     # Biological sequences
data(phiX174Phage)                      # sample data, see ?phiX174Phage
phiX174Phage
m <- consensusMatrix(phiX174Phage)[1:4,] # nucl. x position counts
polymorphic <- which(colSums(m != 0) > 1)
m[, polymorphic]

## ----methods, eval=FALSE-------------------------------------------------
## methods(class=class(phiX174Phage))      # 'DNAStringSet' methods

## ----phiX----------------------------------------------------------------
library(Biostrings)
data(phiX174Phage)

## ----consensusMatrix-----------------------------------------------------
m <- consensusMatrix(phiX174Phage)[1:4,]
polymorphic <- which(colSums(m != 0) > 1)
mapply(substr, polymorphic, polymorphic, MoreArgs=list(x=phiX174Phage))

## ----require-------------------------------------------------------------
library(GenomicRanges)

## ----help, eval=FALSE----------------------------------------------------
## help(package="GenomicRanges")
## vignette(package="GenomicRanges")
## vignette(package="GenomicRanges", "GenomicRangesHOWTOs")
## ?GRanges

