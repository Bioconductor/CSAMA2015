# From: https://plot.ly/ggplot2/getting-started/
#
# library("devtools")
# install_github("ropensci/plotly")

library("ggplot2")
library("plotly")

ggiris <- ggplot(iris, aes(x = Petal.Width, y = Sepal.Length, colour = Species)) + geom_point() + coord_fixed()
print(ggiris)

set_credentials_file("DemoAccount", "lr1c37zw81")
py <- plotly()
r  <- py$ggplotly(ggiris)

print(r$response$url)


## See also
http://badhessian.org/2014/08/a-brief-introduction-to-plotly/
https://github.com/ropensci/plotly
