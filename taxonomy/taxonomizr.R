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
taxons_Genbank<-getTaxonomy(taxaId_Genbank,'accessionTaxa.sql')
table1 <- cbind(table1,taxaId_Genbank, taxons_Genbank)

table1 %>% as_tibble() %>% print(n=Inf, width=Inf)

getDescendants(7742,'accessionTaxa.sql')

table2 <- read.csv("data_freeze_ids.ls", header=FALSE)
taxons_VGP<-getTaxonomy(table,'accessionTaxa.sql')
print(taxons_VGP)

raw <- getRawTaxonomy(7742,'accessionTaxa.sql')
normalized <- normalizeTaxa(raw, lineageOrder=c('infraorder','suborder','superorder','infraclass','subclass','class','subphylum'))

full_tree <- makeNewick(normalized)
write(full_tree, file = "full_tree.nwk")