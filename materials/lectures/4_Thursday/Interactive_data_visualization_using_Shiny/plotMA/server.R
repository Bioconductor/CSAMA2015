### server.R : builds the plots and tables ###
### analysis.R : my additional file ### 

library(DESeq2)
library(pasilla)
library(Biobase)
data(pasillaGenes)

# DE analysis of the pasilla dataset from Bioconductor
dds <- DESeqDataSetFromMatrix(counts(pasillaGenes),
                              pData(pasillaGenes)[,2:3],
                              ~ condition)

# compare treated vs untreated
dds$condition <- relevel(dds$condition, "untreated")
dds <- dds[rowSums(counts(dds)) > 0, ]

# run DESeq
dds <- DESeq(dds)
res <- results(dds)

# this object will be used to locate points from click events.
# take the log of x, so that points are 'close' in the log x axis
data <- with(res, cbind(baseMean, log2FoldChange))

# we set the ylim so need to use
ymax <- 2.5
data[,2] <- pmin(ymax, pmax(-ymax, data[,2]))
scale <- c(diff(range(data[,1])), 2*ymax)
t.data.scaled <- t(data)/scale


library(shiny)
shinyServer(function(input, output) {
  idx = NULL
  xy = reactive(c(input$plotma_click$x, input$plotma_click$y))
  observe({
    if (!is.null(xy())) {
      ## find index of the closest point
      sqdists <- colMeans( (t.data.scaled - xy()/scale )^2 ) 
      idx <<- which.min(sqdists)  
    }
  })
  
  # MA-plot
  output$plotma <- renderPlot({
    # update on user click
    xy()
    par(mar=c(5,5,3,2),cex.lab=2,cex.main=2,cex.axis=1.5)
    # MA-plot of all genes
    plotMA( res, ylim=c( -ymax, ymax ) )
    # add circle for the selected point
    if (!is.null(idx)) points( data[idx,1], data[idx,2], col="dodgerblue", cex=3, lwd=3 )
  })

  # counts plot
  output$plotcounts <- renderPlot({
    # update on user click
    xy()
    par(mar=c(5,5,3,2),cex.lab=2,cex.main=2,cex.axis=1.5)
    # plot the counts for the selected gene
    if (!is.null(idx)) plotCounts( dds, idx )
  })
  
})
