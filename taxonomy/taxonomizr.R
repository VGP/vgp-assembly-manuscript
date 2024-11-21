# tested with Rstudio/Rscript v==4.0

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

#fetch NCBI db
prepareDatabase('accessionTaxa.sql')

#load RefSeq data
table1 <- read.csv("accession_metadata.ls", header=FALSE, sep=",")
names(table1) <- c("Accession", "Tolid", "Species")

taxaId_Genbank<-accessionToTaxa(as.character(table1$Accession),"accessionTaxa.sql")
taxa_Genbank<-getTaxonomy(taxaId_Genbank,'accessionTaxa.sql')
table1 <- cbind(table1,taxaId_Genbank, taxa_Genbank)

#table1 %>% as_tibble() %>% print(n=Inf, width=Inf)

table2 <- read.csv("data_freeze_ids.ls", header=FALSE)
taxa_VGP<-getTaxonomy(table2[,1],'accessionTaxa.sql')
#print(taxa_VGP)

VGP_tree <- makeNewick(taxa_VGP)
write(VGP_tree, file = "VGP_tree.nwk")

descendants_all <- getId(getDescendants(7742,'accessionTaxa.sql') %>% str_subset(pattern = " sp. "), 'accessionTaxa.sql')
taxa_all<-getTaxonomy(descendants_all,'accessionTaxa.sql')

full_tree <- makeNewick(taxa_all)
write(full_tree, file = "full_tree.nwk")