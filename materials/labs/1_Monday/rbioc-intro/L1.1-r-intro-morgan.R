## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()

## ----setup, echo=FALSE---------------------------------------------------
knitr::opts_chunk$set(cache=TRUE)

## ----lapply-args---------------------------------------------------------
args(lapply)

## ----lapply-eg-----------------------------------------------------------
lst <- list(a=1:2, b=2:4)
lapply(lst, log)      # 'base' argument default; natural log
lapply(lst, log, 10)  # '10' is second argument to 'log()', i.e., log base 10

## ------------------------------------------------------------------------
args(mapply)

## ----mapply-eg-----------------------------------------------------------
mapply(seq, 1:3, 4:6, SIMPLIFY=FALSE) # seq(1, 4); seq(2, 5); seq(3, 6)

## ----apply---------------------------------------------------------------
args(apply)

## ---- eval=FALSE---------------------------------------------------------
## if (test) {
##     ## code if TEST == TRUE
## } else {
##     ## code if TEST == FALSE
## }

## ----myfun---------------------------------------------------------------
fun <- function(x) {
    length(unique(x))
}
## list of length 5, each containsing a sample (with replacement) of letters
lets <- replicate(5, sample(letters, 50, TRUE), simplify=FALSE)
sapply(lets, fun)

## ------------------------------------------------------------------------
x <- rnorm(1000)                     # atomic vectors
y <- x + rnorm(1000, sd=.5)
df <- data.frame(x=x, y=y)           # object of class 'data.frame'
plot(y ~ x, df)                      # generic plot, method plot.formula
fit <- lm(y ~x, df)                  # object of class 'lm'
methods(class=class(fit))            # introspection
anova(fit)
plot(y ~ x, df)                      # methods(plot); ?plot.formula
abline(fit, col="red", lwd=3, lty=2) # a function, not generic.method

## ----lapply-user-setup ---------------------------------------
## download example data 'symgo.csv' from course repository
fl <- file.choose()      ## symgo.csv

## ----lapply--------------------------------------------------------------
symgo <- read.csv(fl, row.names=1, stringsAsFactors=FALSE)
head(symgo)
dim(symgo)
length(unique(symgo$SYMBOL))
## split-sapply
go2sym <- split(symgo$SYMBOL, symgo$GO)
len1 <- sapply(go2sym, length)          # compare with lapply, vapply
## built-in functions for common actions
len2 <- lengths(go2sym)
identical(len1, len2)
## smarter built-in functions, e.g., omiting NAs
len3 <- aggregate(SYMBOL ~ GO, symgo, length)
head(len3)
## more fun with aggregate()
head(aggregate(GO ~ SYMBOL, symgo, length))
head(aggregate(SYMBOL ~ GO, symgo, c))
## your own function -- unique, lower-case identifiers
uidfun  <- function(x) {
    unique(tolower(x))
}
head(aggregate(SYMBOL ~ GO , symgo, uidfun))
## as an 'anonymous' function
head(aggregate(SYMBOL ~ GO, symgo, function(x) {
    unique(tolower(x))
}))

## ----echo=TRUE, eval=FALSE-----------------------------------------------
## download 'ALLphenoData.tsv' from course repository
fname <- file.choose()   ## "ALLphenoData.tsv"
stopifnot(file.exists(fname))
pdata <- read.delim(fname)

## ----ALL-properties------------------------------------------------------
class(pdata)
colnames(pdata)
dim(pdata)
head(pdata)
summary(pdata$sex)
summary(pdata$cyto.normal)

## ----ALL-subset----------------------------------------------------------
pdata[1:5, 3:4]
pdata[1:5, ]
head(pdata[, 3:5])
tail(pdata[, 3:5], 3)
head(pdata$age)
head(pdata$sex)
head(pdata[pdata$age > 21,])

## ----ALL-subset-NA-------------------------------------------------------
idx <- pdata$sex == "F" & pdata$age > 40
table(idx)
dim(pdata[idx,])

## ----ALL-BCR/ABL-subset--------------------------------------------------
bcrabl <- pdata[pdata$mol.biol %in% c("BCR/ABL", "NEG"),]

## ----ALL-BCR/ABL-drop-unused---------------------------------------------
bcrabl$mol.biol <- factor(bcrabl$mol.biol)

## ----ALL-BT--------------------------------------------------------------
levels(bcrabl$BT)

## ----ALL-BT-recode-------------------------------------------------------
table(bcrabl$BT)
levels(bcrabl$BT) <- substring(levels(bcrabl$BT), 1, 1)
table(bcrabl$BT)

## ----ALL-BCR/ABL-BT------------------------------------------------------
xtabs(~ BT + mol.biol, bcrabl)

## ----ALL-aggregate-------------------------------------------------------
aggregate(age ~ mol.biol + sex, bcrabl, mean)

## ----ALL-age-------------------------------------------------------------
t.test(age ~ mol.biol, bcrabl)
boxplot(age ~ mol.biol, bcrabl)

## ----echo=TRUE, eval=FALSE-----------------------------------------------
## download 'BRFSS-subset.csv' from course repository
fname <- file.choose()   ## BRFSS-subset.csv
stopifnot(file.exists(fname))
brfss <- read.csv(fname)

## ----brfss-simple-plot---------------------------------------------------
plot(sqrt(Weight) ~ Height, brfss, main="All Years, Both Sexes")

## ----brfss-subset--------------------------------------------------------
brfss2010 <- brfss[brfss$Year == "2010", ]

## ----brfss-pair-plot-----------------------------------------------------
opar <- par(mfcol=c(1, 2))
plot(sqrt(Weight) ~ Height, brfss2010[brfss2010$Sex == "Female", ],
     main="2010, Female")
plot(sqrt(Weight) ~ Height, brfss2010[brfss2010$Sex == "Male", ],
     main="2010, Male")
par(opar)                           # reset 'par' to original value

## ----ggplot2-brfss-simple-plot-------------------------------------------
library(ggplot2)

## 'quick' plot
qplot(Height, sqrt(Weight), data=brfss)

## specify the data set and 'aesthetics', then how to plot
ggplot(brfss, aes(x=Height, y=sqrt(Weight))) +
    geom_point()

## ----ggplot2-na-in-dataset-----------------------------------------------
sum(is.na(brfss$Height))
sum(is.na(brfss$Weight))
drop <- is.na(brfss$Height) | is.na(brfss$Weight)
sum(drop)

## ----ggplot2-remove-na---------------------------------------------------
brfss <- brfss[!drop,]

## ----ggplot2-annotate----------------------------------------------------
qplot(Height, sqrt(Weight), data=brfss) +
    ylab("Square root of Weight") + 
        ggtitle("All Years, Both Sexes")

## ----ggplot2-color-------------------------------------------------------
ggplot(brfss, aes(x=Height, y=sqrt(Weight), color=Sex)) + 
    geom_point()

## ----ggplot2-shape-------------------------------------------------------
ggplot(brfss, aes(x=Height, y = sqrt(Weight), color=Sex, shape=Sex)) + 
    geom_point()

## ----ggplot2-shape-facet-------------------------------------------------
ggplot(brfss, aes(x=Height, y = sqrt(Weight), color=Sex)) + 
    geom_point() +
        facet_grid(Sex ~ .)

## ----ggplot2-subset-facet------------------------------------------------
brfss2010 <- brfss[brfss$Year == "2010", ]
ggplot(brfss2010, aes(x=sqrt(Weight), fill=Sex)) +
    geom_density(alpha=.25)

## ----ggplot2-transparent-------------------------------------------------
sp <- ggplot(brfss, aes(x=Height, y=sqrt(Weight)))
sp + geom_point(alpha=.4)

## ----ggplot2-regression--------------------------------------------------
sp + geom_point() + stat_smooth(method=lm)

## ----ggplot2-regression-2, eval=FALSE------------------------------------
## sp + geom_point() + stat_smooth(method=lm + level=0.95)
## sp + geom_point() + stat_smooth(method=lm, se=FALSE)

## ----ggplot2-regression-bygroup------------------------------------------
sps <- ggplot(brfss, aes(x=Height, y=sqrt(Weight), colour=Sex)) +
    geom_point() +
        scale_colour_brewer(palette="Set1")
sps + geom_smooth(method="lm")

