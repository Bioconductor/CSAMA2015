### server.R : builds the plots and tables ###
library(DESeq2)

load("data.rda")
res <- results(dds)

ymax <- 5

# this object will be used to locate points from click events.
data <- with(res, cbind(baseMean, log2FoldChange))
data[,2] <- pmin(ymax, pmax(-ymax, data[,2]))
scale <- c(diff(range(data[,1])), 2*ymax)
t.data.scaled <- t(data)/scale

par(mar=c(5,5,3,2), cex.lab=2, cex.main=2, cex.axis=1.5)

shinyServer(function(input, output) {
  
  current = reactiveValues(idx = NULL)
  
  xy = reactive(c(input$plotma_click$x, input$plotma_click$y))
  
  observe({
    if (!is.null(xy())) {
      ## find index of the closest point
      sqdists <- colMeans( (t.data.scaled - xy()/scale )^2 ) 
      current$idx <- which.min(sqdists)  
    }
  })
  
  # MA-plot
  output$plotma <- renderPlot({
    # update on user click
    #xy()
    #par(mar=c(5,5,3,2), cex.lab=2, cex.main=2, cex.axis=1.5)
    # MA-plot of all genes
    plotMA(res, ylim=c( -ymax, ymax ), alpha=input$alpha)
    # add circle for the selected point
    idx = current$idx
    if (!is.null(idx)) points( data[idx,1], data[idx,2], col="dodgerblue", cex=3, lwd=3 )
  })

  # counts plot
  output$plotcounts <- renderPlot({
    par(mar=c(5,5,3,2), cex.lab=2, cex.main=2, cex.axis=1.5)
    # plot the counts for the selected gene
    # update only when idx changes
    idx = current$idx
    if (!is.null(idx)) plotCounts( dds, idx, intgroup=c("dex") )
  })
  
})
