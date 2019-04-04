#main file
library(ggplot2)
library(dplyr)
library(tidyr)

#Include the data file
source(file = "data/data.R")


#Have a look at the plots
p <- ggplot(data, aes(y=y,x=x))
p+geom_point(aes(x,y), size = 2, colour="#CC0000")+stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)+labs(title="Quadratic regression model", x ="Number of telephone calls", y = "Duration of telephone calls",subtitle = "Model fitted wiht a quadratic approximation")

ggsave("fig/quadraticregression.png", width = 10, height = 8)
