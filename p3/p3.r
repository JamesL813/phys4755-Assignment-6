library(ggplot2)

data <- read.csv("p3.csv", header = TRUE, sep = " ")
#create bar plot with error bars
ggplot(data) +
    geom_point(aes(x=steps, y=S), stat='identity', color='steelblue') +
    geom_errorbar(aes(x=steps, , ymin=S-error, ymax=S+error), width=0.4)
ggsave(file = "p3.pdf")    
