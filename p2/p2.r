library(ggplot2)

data <- read.csv("crank.csv", header = TRUE, sep = " ")
#create bar plot with error bars
ggplot(data) +
    geom_point(aes(x=x, y=t), stat='identity', color='red')
ggsave(file = "crank.png")    
