library("scatterplot3d") # load

data_a <- read.csv("p1a.csv", header = TRUE, sep = " ")
png(filename = "p1a.png", width = 8, height = 6, units = "in", res = 300)
scatterplot3d(data_a$x, data_a$y, data_a$phi,
    # color="gray20",
    box = TRUE,
    main="(a) Results",
    xlab = "x",
    ylab = "y",
    zlab = "phi")

data_b <- read.csv("p1b.csv", header = TRUE, sep = " ")
png(filename = "p1b.png", width = 8, height = 6, units = "in", res = 300)
scatterplot3d(data_b$x, data_b$y, data_b$phi,
    # color="gray20",
    box = TRUE,
    main="(b) Results",
    xlab = "x",
    ylab = "y",
    zlab = "phi")
