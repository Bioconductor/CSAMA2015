### ui.R : lays out the page ###

library(shiny)
shinyUI(fluidPage(
  titlePanel("MA-plot + counts plot with shiny"),

  # flow layout fills out left to right then down
  splitLayout(

    # the MA-plot
    plotOutput("plotma", click="plotma_click", width=360, height=360),

    # the counts plot for a single gene
    plotOutput("plotcounts", width=360, height=360),
    
    # needed for proper page layout
    cellArgs = list(style="width: 360px;")    
  )
  
))
