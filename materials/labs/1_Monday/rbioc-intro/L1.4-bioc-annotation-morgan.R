## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()

## ----setup, echo=FALSE, warning=FALSE------------------------------------
options(max.print=1000, width=100)
knitr::opts_chunk$set(cache=TRUE)
suppressPackageStartupMessages({
    library(org.Hs.eg.db)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(BSgenome.Hsapiens.UCSC.hg19)
    library(GenomicRanges)
    library(biomaRt)
    library(rtracklayer)
    library(Gviz)
    library(AnnotationHub)
})

## ----select-setup--------------------------------------------------------
ensids <- c("ENSG00000130720", "ENSG00000103257", "ENSG00000156414", 
            "ENSG00000144644", "ENSG00000159307", "ENSG00000144485")

## ----select--------------------------------------------------------------
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
columns(org.Hs.eg.db)
cols <- c("SYMBOL", "GENENAME")
select(org.Hs.eg.db, keys=ensids, columns=cols, keytype="ENSEMBL")

## ----biomaRt1, eval=FALSE, results="hide"--------------------------------
## ## NEEDS INTERNET ACCESS !!
## library(biomaRt)
## head(listMarts(), 3)                      ## list the marts
## head(listDatasets(useMart("ensembl")), 3) ## mart datasets
## ensembl <-                                ## fully specified mart
##     useMart("ensembl", dataset = "hsapiens_gene_ensembl")
## 
## head(listFilters(ensembl), 3)             ## filters
## myFilter <- "chromosome_name"
## substr(filterOptions(myFilter, ensembl), 1, 50) ## return values
## myValues <- c("21", "22")
## head(listAttributes(ensembl), 3)          ## attributes
## myAttributes <- c("ensembl_gene_id","chromosome_name")
## 
## ## assemble and query the mart
## res <- getBM(attributes =  myAttributes, filters =  myFilter,
##              values =  myValues, mart = ensembl)

## ----symbol-to-entrez----------------------------------------------------
library(org.Hs.eg.db)
eid <- select(org.Hs.eg.db, "BRCA1", "ENTREZID", "SYMBOL")[["ENTREZID"]]

## ----entrez-to-tx, message=FALSE-----------------------------------------
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txid <- select(txdb, eid, "TXNAME", "GENEID")[["TXNAME"]]

## ----tx-to-cds-coords----------------------------------------------------
cds <- cdsBy(txdb, by="tx", use.names=TRUE)
brca1cds <- cds[names(cds) %in% txid]
class(brca1cds)
length(brca1cds)
brca1cds[[1]]                           # exons in cds
cdswidth <- width(brca1cds)             # width of each exon
all((sum(cdswidth) %% 3) == 0)          # sum within cds, modulus 3

## ----Gviz, message=FALSE-------------------------------------------------
library(Gviz)
grt <- GeneRegionTrack(txdb)
plotTracks(list(GenomeAxisTrack(), grt), chromosome="chr17",
    from=min(start(unlist(brca1cds))),
    to=max(end(unlist(brca1cds))))

## ----cds-to-seq----------------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19
tx_seq <- extractTranscriptSeqs(genome, brca1cds)
tx_seq

## ----introns-------------------------------------------------------------
introns <- psetdiff(range(brca1cds), brca1cds)

## ----intron-seqs---------------------------------------------------------
seq <- getSeq(genome, introns)
names(seq)
seq[["uc010whl.2"]]                     # 21 introns

## ----rtracklayer-roi-----------------------------------------------------
library(GenomicRanges)
roi <- GRanges("chr10", IRanges(92106877, 112106876, names="ENSG00000099194"))

## ----rtracklayer-session, eval=FALSE-------------------------------------
## library(rtracklayer)
## session <- browserSession()

## ----rtracklayer-marks, eval=FALSE---------------------------------------
## trackName <- "wgEncodeRegTfbsClusteredV2"
## tableName <- "wgEncodeRegTfbsClusteredV2"
## trFactor <- "ERalpha_a"
## ucscTable <- getTable(ucscTableQuery(session, track=trackName,
##     range=roi, table=tableName, name=trFactor))

## ----rtracklayer-plot, fig.height=3, eval=FALSE--------------------------
## plot(score ~ chromStart, ucscTable, pch="+")
## abline(v=start(roi) + (end(roi) - start(roi) + 1) / 2, col="blue")

