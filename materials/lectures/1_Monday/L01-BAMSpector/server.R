##  server setup

## load required libraries
library(shiny)
library(RNAseqData.HNRNPC.bam.chr14)
library(Homo.sapiens)
library(Gviz)

## where are the BAM files?
dirname <- unique(dirname(RNAseqData.HNRNPC.bam.chr14_BAMFILES))

## What are the ranges of each gene?
ranges <- genes(Homo.sapiens, columns="SYMBOL")
ranges$SYMBOL <- unlist(ranges$SYMBOL)

## Create a representation of each gene region
genes <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene,
                         chromosome="chr14")

shinyServer(function(input, output) {

    output$tracksPlot <- renderPlot({
        if (length(input$bam) > 0) {
            ## coverage on each BAM file
            bam <- file.path(dirname, input$bam)
            coverage <- Map(DataTrack,
                range = bam, name = bam,
                MoreArgs=list(type = 'histogram',
                    window = -1, genome = 'hg19',
                    chromosome = 'chr14'))
        } else {
            coverage <- list()
        }

        ## Select the correct range
        range <- ranges[match(input$symbol, ranges$SYMBOL)]

        ## plot the GeneRegionTrack and coverage
        plotTracks(c(list(genes), coverage),
                   from = start(range), to=end(range),
                   chr='chr14', windowSize = 30)
    })
})
