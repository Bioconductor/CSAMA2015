library(shiny)
library(RNAseqData.HNRNPC.bam.chr14)
library(Homo.sapiens)

## Get all SYMBOLs on chr14
symbols <- keys(Homo.sapiens, keytype="SYMBOL")
map <- select(Homo.sapiens, symbols, "TXCHROM", "SYMBOL")
symchoices <- sort(unique(map$SYMBOL[map$TXCHROM %in% "chr14"]))

## Possible BAM files
bamchoices <- basename(RNAseqData.HNRNPC.bam.chr14_BAMFILES)

## Define the user interface
shinyUI(fluidPage(

    ## Application title
    titlePanel("BAMSpector: Reads Supporting Gene Models"),

    sidebarLayout(
        sidebarPanel(
            ## input gene symbol (fancy: select from available)
            selectInput("symbol", "Gene Symbol", symchoices),

            ## input path to BAM file
            selectInput("bam", "BAM File", bamchoices, multiple=TRUE)),

        ## Show a plot of the generated distribution
        mainPanel(plotOutput("tracksPlot")))
    ))
