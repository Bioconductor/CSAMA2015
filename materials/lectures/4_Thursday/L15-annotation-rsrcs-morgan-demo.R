## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()
options(showHeadLines=3, showTailLines=2, max.print=1000)
suppressPackageStartupMessages({
    library(biomaRt)
    library(AnnotationHub)
})

## ----setup, echo=FALSE---------------------------------------------------
knitr::opts_chunk$set(cache=TRUE)

## ----AnnotationHub-takifugu----------------------------------------------
library(AnnotationHub)
hub = AnnotationHub()
query(hub, c("ensembl","release-80", "Takifugu"))
gtf <- hub[["AH47101"]]
dna <- hub[["AH47477"]]

gtf
dna
head(seqlevels(dna))

