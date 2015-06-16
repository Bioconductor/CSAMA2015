### ui.R : lays out the page ###
size = 450

shinyUI(fluidPage(
  titlePanel("MA plot explorer"),

  splitLayout(cellWidths=size,
    plotOutput("plotma", click="plotma_click", width=size, height=size),
    plotOutput("plotcounts", width=size, height=size)
  ),
  splitLayout(cellWidths=size,
    sliderInput("alpha", "Alpha", min=0, max=0.2, value=0.1, step=0.001, width=size)
  )
  
))
