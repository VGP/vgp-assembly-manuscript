requiredpackages <- c("this.path", "data.table", "taxonomizr", "tidyverse")

install_load <- function(packages){
     for (p in packages) {
          if (p %in% rownames(installed.packages())) {
               library(p, character.only=TRUE)
          } else {
               install.packages(p, repos = "http://cran.us.r-project.org")
               library(p,character.only = TRUE)
          }
     }
}

install_load(requiredpackages)

args <- commandArgs(trailingOnly = TRUE)

#fetch NCBI db
prepareDatabase('accessionTaxa.sql')

#load data
taxaIds <- read.csv(args[1], header=FALSE)
names(taxaIds) <- c("Taxon")
taxa <- getTaxonomy(taxaIds$Taxon,'accessionTaxa.sql')
write(taxa, file = args[2])