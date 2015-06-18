## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()
options(showHeadLines=3, showTailLines=2, max.print=1000)
suppressPackageStartupMessages({
    library(biomaRt)
    library(AnnotationHub)
    library(org.Hs.eg.db)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(Rsamtools)
})

## ----setup, echo=FALSE---------------------------------------------------
#knitr::opts_chunk$set(cache=TRUE)

## ----OrgDb---------------------------------------------------------------
library(org.Hs.eg.db)
org <- org.Hs.eg.db       # convenient alias
keytypes(org)             # map from keys...
columns(org)              # ...to columns

## ----1-1-mapping---------------------------------------------------------
sym <- keys(org, "SYMBOL")
mapIds(org, c("BRCA1", "PTEN"), "GENENAME", "SYMBOL")
map <- mapIds(org, sym, "GENENAME", "SYMBOL")
length(map)
head(map)

## ----1-many-mapping------------------------------------------------------
head(select(org, keys(org), "ALIAS"))
head(mapIds(org, keys(org), "ALIAS", "ENTREZID"))
head(mapIds(org, keys(org), "ALIAS", "ENTREZID", multiVals="CharacterList"))
str(head(mapIds(org, keys(org), "ALIAS", "ENTREZID", multiVals="list")))

## ----org-and-txdb--------------------------------------------------------
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

sym <- "BRCA1"

eid <- mapIds(org.Hs.eg.db, sym, "ENTREZID", "SYMBOL")
txid <- mapIds(txdb, eid, "TXNAME", "GENEID", multiVals="list")[[eid]]
txid

cds <- cdsBy(txdb, by="tx", use.names=TRUE)
cds[names(cds) %in% txid]

## ----organismdb, eval=FALSE----------------------------------------------
## library(Homo.sapiens)
## 
## txid <- mapIds(Homo.sapiens, sym, "TXNAME", "SYMBOL", multiVals="list")[[sym]]
## cds <- cdsBy(Homo.sapiens, by="tx", use.names=TRUE)

## ----biomaRt-------------------------------------------------------------
library(biomaRt)
head(listMarts())                               # available marts, 52!
head(listDatasets(useMart("ensembl")))          # datasets in mart, 69!
ensembl <-                                      # fully specified mart
  useMart("ensembl", dataset = "hsapiens_gene_ensembl")

head(listFilters(ensembl), 3)                   # filters, 296!
myFilter <- "chromosome_name"
substr(filterOptions(myFilter, ensembl), 1, 50) # return values
myValues <- c("21", "22")
head(listAttributes(ensembl), 3)                # attributes
myAttributes <- c("ensembl_gene_id","chromosome_name")

## assemble and query the mart
res <- getBM(myAttributes, myFilter, myValues, ensembl)

## ----AnnotationHub-------------------------------------------------------
library(AnnotationHub)
hub = AnnotationHub()
hub
query(hub, c("Ensembl", "80", "gtf"))
## ensgtf = display(hub)                   # visual choice
hub["AH47107"]
gtf <- hub[["AH47107"]]
gtf
## org.* data bases available from AnnotationHub
query(hub, "OrgDb")
mcols(query(hub, "OrgDb"))[, "species", drop=FALSE]

## ----AnnotationHub-org---------------------------------------------------
library(AnnotationHub)
hub = AnnotationHub()
query(hub, "OrgDb")
hub[["AH12818"]]

## ----AnnotationHub-takifugu----------------------------------------------
library(AnnotationHub)
hub = AnnotationHub()
query(hub, c("ensembl","release-80", "Takifugu"))
gtf <- hub[["AH47101"]]
dna <- hub[["AH47477"]]

gtf
dna
head(seqlevels(dna))

## ----AnnotationHub-txdb--------------------------------------------------
library(GenomicFeatures)
txdb <- makeTxDbFromGRanges(gtf)
## saveDb(txdb, "TxDb.Takifugu.Ensembl.80.sqlite")
## loadDb("TxDb.Takifugu.Ensembl.80.sqlite")

## ----AnnotationHub-takifugi-exons----------------------------------------
library(Rsamtools)     # for getSeq,FaFile-method
exons <- exons(txdb)
getSeq(dna, exons)

