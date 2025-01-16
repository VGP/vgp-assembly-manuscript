setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# library
library(ggplot2)
library(ggExtra)

rdeval_data <- read.csv("rdeval.tsv", sep="\t", header=T)

# classic plot :
p <- ggplot(rdeval_data, aes(x=Average.read.length, y=Average.read.quality)) +
  geom_point() +
  theme(legend.position="none")

# with marginal histogram
p1 <- ggMarginal(p, type="histogram")

# marginal density
p2 <- ggMarginal(p, type="density")

# marginal boxplot
p3 <- ggMarginal(p, type="boxplot")