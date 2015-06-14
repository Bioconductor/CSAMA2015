## ----env0, message=FALSE, echo=FALSE, warning=FALSE----------------------
library("knitr")
opts_knit$set(error = FALSE)
library("BiocInstaller")
library("RforProteomics")
library("nloptr")

## ---- env, message=FALSE, echo=TRUE, warning=FALSE-----------------------
library("mzR")
library("mzID")
library("MSnID")
library("MSnbase")
library("rpx")
library("MLInterfaces")
library("pRoloc")
library("pRolocdata")
library("MSGFplus")
library("rols")
library("hpar")

## ----r4pinstall, eval=FALSE----------------------------------------------
## library("BiocInstaller")
## biocLite("RforProteomics", dependencies = TRUE)

## ---- pk, echo=FALSE, warning=FALSE, cache=TRUE--------------------------
biocv <- as.character(biocVersion())
pp <- proteomicsPackages(biocv)
msp <- massSpectrometryPackages(biocv)
msdp <- massSpectrometryDataPackages(biocv)

## ---- pp, eval=FALSE-----------------------------------------------------
## library("RforProteomics")
## pp <- proteomicsPackages()
## display(pp)

## ---- datatab, results='asis', echo=FALSE--------------------------------
datatab <-
    data.frame(Type = c("raw", "identification", "quantitation",
                   "peak lists", "other"),
               Format = c("mzML, mzXML, netCDF, mzData",
                   "mzIdentML", "mzQuantML", "mgf", "mzTab"),
               Package = c(
                   "[`mzR`](http://bioconductor.org/packages/release/bioc/html/mzR.html) (read)",
                   paste("[`mzR`](http://bioconductor.org/packages/release/bioc/html/mzR.html) (read) and",
                         "[`mzID`](http://bioconductor.org/packages/release/bioc/html/mzID.html) (read)"),
                   "",
                   "[`MSnbase`](http://bioconductor.org/packages/release/bioc/html/MSnbase.html) (read/write)", 
                   "[`MSnbase`](http://bioconductor.org/packages/release/bioc/html/MSnbase.html) (read/write)"))
library("knitr")
kable(datatab)

## ---- rpx----------------------------------------------------------------
library("rpx")
pxannounced()

## ---- pxd, cache=TRUE----------------------------------------------------
px <- PXDataset("PXD000001")
px
pxfiles(px)

## ---- pxvar--------------------------------------------------------------
pxtax(px)
pxurl(px)
pxref(px)

## ---- pxget, eval=FALSE--------------------------------------------------
## mzf <- pxget(px, pxfiles(px)[6])

## ----pxd1----------------------------------------------------------------
library("PXD000001")
px1files()
mzf <- px1files()[5]
px1get(mzf)

## ---- rawms--------------------------------------------------------------
library("mzR")
ms <- openMSfile(mzf)
ms

## ---- hd-----------------------------------------------------------------
hd <- header(ms)
dim(hd)
names(hd)

## ------------------------------------------------------------------------
hd[1000, ]
head(peaks(ms, 1000))
plot(peaks(ms, 1000), type = "h")

## ---- msmap--------------------------------------------------------------
## a set of spectra of interest: MS1 spectra eluted
## between 30 and 35 minutes retention time
ms1 <- which(hd$msLevel == 1)
rtsel <- hd$retentionTime[ms1] / 60 > 30 &
    hd$retentionTime[ms1] / 60 < 35

## the map
M <- MSmap(ms, ms1[rtsel], 521, 523, .005, hd)

plot(M, aspect = 1, allTicks = FALSE)
plot3D(M)

## With some MS2 spectra
i <- ms1[which(rtsel)][1]
j <- ms1[which(rtsel)][2]
M2 <- MSmap(ms, i:j, 100, 1000, 1, hd)
plot3D(M2)

## ---- id, cache=TRUE-----------------------------------------------------
library("mzID")
f <- dir(system.file("extdata", package = "RforProteomics"),
         pattern = "mzid", full.names=TRUE)
basename(f)
id <- mzID(f)
id

## ----mzrvsid, eval = TRUE------------------------------------------------
library("mzR")
f <- dir(system.file("extdata", package = "RforProteomics"),
         pattern = "mzid", full.names=TRUE)

id1 <- openIDfile(f)
fid1 <- mzR::psms(id1)

head(fid1)

## ---- rtandem, eval=FALSE------------------------------------------------
## library("rTANDEM")
## ?rtandem
## library("shinyTANDEM")
## ?shinyTANDEM

## ----ex_getfas, eval=FALSE-----------------------------------------------
## fas <- pxget(px, pxfiles(px)[10])
## basename(fas)

## ----ex_getfaslocal------------------------------------------------------
fas <- px1files()[1]
px1get(fas)

## ----ex_msgfplus, message=FALSE, cache=TRUE------------------------------
library("MSGFplus")
msgfpar <- msgfPar(database = fas,
                   instrument = 'HighRes',
                   tda = TRUE,
                   enzyme = 'Trypsin',
                   protocol = 'iTRAQ')
idres <- runMSGF(msgfpar, mzf, memory=1000)
idres
## identification file (needed below)
basename(mzID::files(idres)$id)

## ----msgfgui, eval=FALSE-------------------------------------------------
## library("MSGFgui")
## MSGFgui()

## ---- msnid--------------------------------------------------------------
library("MSnID")
msnid <- MSnID(".")
msnid <- read_mzIDs(msnid,
                    basename(mzID::files(idres)$id))
show(msnid)

## ----msnvars-------------------------------------------------------------
msnid <- correct_peak_selection(msnid)
msnid$msmsScore <- -log10(msnid$`MS-GF:SpecEValue`)
msnid$absParentMassErrorPPM <- abs(mass_measurement_error(msnid))

## ----idplot, echo=FALSE--------------------------------------------------
x <- psms(msnid)$msmsScore
gr <- psms(msnid)$isDecoy
library("lattice")
densityplot(x,
            group = gr,
            auto.key=TRUE)

## ----msnidfilt-----------------------------------------------------------
filtObj <- MSnIDFilter(msnid)
filtObj$absParentMassErrorPPM <- list(comparison="<", threshold=5.0)
filtObj$msmsScore <- list(comparison=">", threshold=8.0)
filtObj

evaluate_filter(msnid, filtObj)

## ----filtopt-------------------------------------------------------------
filtObj.grid <- optimize_filter(filtObj, msnid, fdr.max=0.01,
                                method="Grid", level="PSM",
                                n.iter=50000)
filtObj.grid
evaluate_filter(msnid, filtObj.grid)

## ----applyfilt-----------------------------------------------------------
msnid <- apply_filter(msnid, filtObj.grid)
msnid

## ---- msnbase------------------------------------------------------------
library("MSnbase")
rawFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
               full.name = TRUE, pattern = "mzXML$")
basename(rawFile)
msexp <- readMSData(rawFile, verbose = FALSE)
msexp

## ------------------------------------------------------------------------
length(msexp)
msexp[1:2]
msexp[[2]]

## ---- addid--------------------------------------------------------------
fData(msexp)
## find path to a mzIdentML file
identFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                 full.name = TRUE, pattern = "dummyiTRAQ.mzid")
basename(identFile)
msexp <- addIdentificationData(msexp, identFile)
fData(msexp)

## ---- specplot-----------------------------------------------------------
msexp[[1]]
plot(msexp[[1]], full=TRUE)

## ---- specplot2----------------------------------------------------------
msexp[1:3]
plot(msexp[1:3], full=TRUE)

## ---- quanttab, echo=FALSE, results='asis'-------------------------------

qtb <- matrix(c("XIC", "Counting", "SILAC, 15N", "iTRAQ, TMT"),
              nrow = 2, ncol = 2)
dimnames(qtb) <- list(
    'MS level' = c("MS1", "MS2"),
    'Quantitation' = c("Label-free", "Labelled"))

kable(qtb)


## ---- itraq4plot---------------------------------------------------------
plot(msexp[[1]], full=TRUE, reporters = iTRAQ4)

## ---- quantitraq---------------------------------------------------------
msset <- quantify(msexp, method = "trap", reporters = iTRAQ4, verbose=FALSE)
exprs(msset)
processingData(msset)

## ---- lfms2--------------------------------------------------------------
exprs(si <- quantify(msexp, method = "SIn"))     
exprs(saf <- quantify(msexp, method = "NSAF"))

## ----mztab, cache=TRUE, eval = FALSE-------------------------------------
## mztf <- pxget(px, pxfiles(px)[2])
## (mzt <- readMzTabData(mztf, what = "PEP"))

## ----mztablocal----------------------------------------------------------
mztf <- px1files()[3]
px1get(mztf)
(mzt <- readMzTabData(mztf, what = "PEP"))

## ---- readmsnset2--------------------------------------------------------
csv <- dir(system.file ("extdata" , package = "pRolocdata"),
           full.names = TRUE, pattern = "pr800866n_si_004-rep1.csv")
getEcols(csv, split = ",")
ecols <- 7:10
res <- readMSnSet2(csv, ecols)
head(exprs(res))
head(fData(res))

## ---- pure---------------------------------------------------------------
data(itraqdata)
qnt <- quantify(itraqdata, method = "trap",
                reporters = iTRAQ4, verbose = FALSE)
impurities <- matrix(c(0.929,0.059,0.002,0.000,
                       0.020,0.923,0.056,0.001,
                       0.000,0.030,0.924,0.045,
                       0.000,0.001,0.040,0.923),
                     nrow=4, byrow = TRUE)
## or, using makeImpuritiesMatrix()
## impurities <- makeImpuritiesMatrix(4)
qnt.crct <- purityCorrect(qnt, impurities)
processingData(qnt.crct)

## ---- pureplot, echo=FALSE-----------------------------------------------
plot0 <- function(x, y, main = "") {
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    par(mar = c(4, 4, 1, 1))
    par(mfrow = c(2, 2))
    sx <- sampleNames(x)
    sy <- sampleNames(y)
    for (i in seq_len(ncol(x))) {
        plot(exprs(x)[, i], exprs(y)[, i], log = "xy",
             xlab = sx[i], ylab = sy[i])
        abline(0, 1)
        grid()
    }
}
## plot0(qnt, qnt.crct, main = "Before/after correction")

## ---- norm---------------------------------------------------------------
qnt.crct.nrm <- normalise(qnt.crct, "quantiles") 

## ---- plotnorm, echo=FALSE-----------------------------------------------
## plot0(qnt, qnt.crct.nrm)

## ---- comb---------------------------------------------------------------
## arbitraty grouping
g <- factor(c(rep(1, 25), rep(2, 15), rep(3, 15)))
g
prt <- combineFeatures(qnt.crct.nrm, groupBy = g, fun = "sum")
processingData(prt)

## ----impute0-------------------------------------------------------------
set.seed(1)
qnt0 <- qnt
exprs(qnt0)[sample(prod(dim(qnt0)), 10)] <- NA
table(is.na(qnt0))
image(qnt0)

## ----filterNA------------------------------------------------------------
## remove features with missing values
qnt00 <- filterNA(qnt0)
dim(qnt00)
any(is.na(qnt00))

## ----impute--------------------------------------------------------------
## impute missing values using knn imputation
qnt.imp <- impute(qnt0, method = "knn")
dim(qnt.imp)
any(is.na(qnt.imp))

## ---- ml-----------------------------------------------------------------
library("MLInterfaces")
library("pRoloc")
library("pRolocdata")
data(dunkley2006)
traininds <- which(fData(dunkley2006)$markers != "unknown")
ans <- MLearn(markers ~ ., data = t(dunkley2006), knnI(k = 5), traininds)
ans

## ----clust---------------------------------------------------------------
kcl <- MLearn( ~ ., data = dunkley2006, kmeansI, centers = 12)
kcl
plot(kcl, exprs(dunkley2006))

## ----nont, echo=FALSE, cache=TRUE----------------------------------------
library("rols")
nont <- nrow(ontologies())

## ----rols----------------------------------------------------------------
library("rols")
olsQuery("ESI", "MS")

## ---- si, echo=FALSE-----------------------------------------------------
print(sessionInfo(), local = FALSE)

