#fetch NCBI db
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

if(!require(taxonomizr)){
    install.packages("taxonomizr")
    library(taxonomizr)
}

prepareDatabase('accessionTaxa.sql')

#load RefSeq data
table1 <- read.csv("accession_metadata.ls", header=FALSE, sep=",")
names(table1) <- c("Accession", "Tolid", "Species")

taxaId_Genbank<-accessionToTaxa(as.character(table1$Accession),"accessionTaxa.sql")
taxons_Genbank<-getTaxonomy(taxaId_Genbank,'accessionTaxa.sql')
table1 <- cbind(table1,taxaId_Genbank, taxons_Genbank)
print(table1)