library('ggplot2')
library(dplyr)
library(scales)
library(tidyr)
library(RColorBrewer)
library('svglite')

testudinata_compleasm<-read.table(file = 'testudinata_compleasmStats.txt',sep = '\t',header = T)
testudinata_asm<-read.table(file = 'testudinata_asmStats.txt',sep='\t',header=T)
testudinata_cnStatsFiltered<-testudinata_asm[,c(1,4,5,6)]
testudinata_merge<-merge(testudinata_compleasm,testudinata_cnStatsFiltered,by="Accession",all.x=TRUE)
# Create the dataframe with proper column assignment
testudinata_df <- data.frame(
  seq = testudinata_merge$Sequencing_type,
  fragmented = testudinata_merge$Fragmented_compleasm,
  frameshift = testudinata_merge$Frameshift_compleasm,
  missing = testudinata_merge$Missing_compleasm,
  complete = testudinata_merge$Complete_compleasm,
  date = testudinata_merge$Submission_Date,
  class = testudinata_merge$Class,
  superclass = testudinata_merge$Super_class,
  superorder = testudinata_merge$Super_order,
  infraclass = testudinata_merge$Infra_class,
  order = testudinata_merge$Order,
  dup = testudinata_merge$Duplicated_compleasm,
  single = testudinata_merge$Single_compleasm,
  proj='Other',
  cn50 =testudinata_merge$Contig_N50
)
crocodylia_compleasm<-read.table(file = 'crocodylia_compleasmStats.txt',sep = '\t',header = T)
crocodylia_asm<-read.table(file = 'crocodylia_asmStats.txt',sep='\t',header=T)
crocodylia_cnStatsFiltered<-crocodylia_asm[,c(1,4,5,6)]
crocodylia_merge<-merge(crocodylia_compleasm,crocodylia_cnStatsFiltered,by="Accession",all.x=TRUE)
# Create the dataframe with proper column assignment
crocodylia_df <- data.frame(
  seq = crocodylia_merge$Sequencing_type,
  fragmented = crocodylia_merge$Fragmented_compleasm,
  frameshift = crocodylia_merge$Frameshift_compleasm,
  missing = crocodylia_merge$Missing_compleasm,
  complete = crocodylia_merge$Complete_compleasm,
  date = crocodylia_merge$Submission_Date,
  class = crocodylia_merge$Class,
  superclass = crocodylia_merge$Super_class,
  superorder = crocodylia_merge$Super_order,
  infraclass = crocodylia_merge$Infra_class,
  order = crocodylia_merge$Order,
  dup = crocodylia_merge$Duplicated_compleasm,
  single = crocodylia_merge$Single_compleasm,
  proj='Other',
  cn50 =crocodylia_merge$Contig_N50
)
lepidosauria_compleasm<-read.table(file = 'lepidosauria_compleasmStats.txt',sep = '\t',header = T)
lepidosauria_asm<-read.table(file = 'lepidosauria_asmStats.txt',sep='\t',header=T)
lepidosauria_cnStatsFiltered<-lepidosauria_asm[,c(1,4,5,6)]
lepidosauria_merge<-merge(lepidosauria_compleasm,lepidosauria_cnStatsFiltered,by="Accession",all.x=TRUE)
# Create the dataframe with proper column assignment
lepidosauria_df <- data.frame(
  seq = lepidosauria_merge$Sequencing_type,
  fragmented = lepidosauria_merge$Fragmented_compleasm,
  frameshift = lepidosauria_merge$Frameshift_compleasm,
  missing = lepidosauria_merge$Missing_compleasm,
  complete = lepidosauria_merge$Complete_compleasm,
  date = lepidosauria_merge$Submission_Date,
  class = lepidosauria_merge$Class,
  superclass = lepidosauria_merge$Super_class,
  superorder = lepidosauria_merge$Super_order,
  infraclass = lepidosauria_merge$Infra_class,
  order = lepidosauria_merge$Order,
  dup = lepidosauria_merge$Duplicated_compleasm,
  single = lepidosauria_merge$Single_compleasm,
  proj='Other',
  cn50 =lepidosauria_merge$Contig_N50
)

reptile_df<-rbind(testudinata_df,crocodylia_df,lepidosauria_df)
reptile_submitter<-c(testudinata_compleasm$Submitter,crocodylia_compleasm$Submitter,lepidosauria_compleasm$Submitter)
for (i in c(1:length(reptile_submitter))){
  if (reptile_submitter[i]=="B10K Consortium"){
    reptile_df$proj[i]='B10K'
  }else if(reptile_submitter[i]=="G10K"){
    reptile_df$proj[i]='VGP'
  }else if(reptile_submitter[i]=="Iridian Genomes"){
    reptile_df$proj[i]='Iridian'
  }else if(reptile_submitter[i]=="IRIDIAN GENOMES"){
    reptile_df$proj[i]='Iridian'
  }else if(reptile_submitter[i]=="Vertebrate Genomes Project"){
    reptile_df$proj[i]='VGP'
  }else if(reptile_submitter[i]=="WELLCOME SANGER INSTITUTE"){
    reptile_df$proj[i]='VGP'
  }else if(reptile_submitter[i]=="Wellcome Sanger Institute"){
    reptile_df$proj[i]='VGP'
  }else if(reptile_submitter[i]=="The Max Planck Institute of Molecular Cell Biology and Genetics"){
    reptile_df$proj[i]='VGP'
  }else if(reptile_submitter[i]=="Zhejiang University"){
    reptile_df$proj[i]='VGP'
  }else if(reptile_submitter[i]=="BGI"){
    reptile_df$proj[i]='BGI'
  }else if(reptile_submitter[i]=="BGI-Shenzhen"){
    reptile_df$proj[i]='BGI'
  }
}


ggplot(data = reptile_df, aes(x = order, y = frameshift, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +ylim(0,500)+
  ggtitle('Compleasm genes containing frameshifts') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = reptile_df, aes(x = order, y = complete, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  scale_fill_manual(values = colors) +
  ggtitle('Compleasm complete genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = reptile_df, aes(x = order, y = single, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  scale_fill_manual(values = colors) +
  ggtitle('Compleasm single-copy genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = reptile_df, aes(x = order, y = single, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  scale_fill_manual(values = colors) +
  ggtitle('Compleasm duplicate genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = reptile_df, aes(x = order, y = fragmented, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  ggtitle('Compleasm fragmented genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##############################################################################
reptile_dfFiltered<-reptile_df[reptile_df$proj %in% c('VGP','B10K','BGI','Iridian','Other'),]
reptile_dfFiltered<-reptile_df[reptile_df$proj %in% c('VGP','B10K','BGI','Iridian','Other'),]
reptile_dfFiltered$proj <- factor(reptile_dfFiltered$proj, 
                               levels = c('VGP','BGI','B10K','Other', 'Iridian'))
dark_palette <- c('Other' = 'orange', 
                  'Iridian' = 'lightblue', 
                  'VGP' = 'darkgreen', 
                  'B10K' = 'darkred', 
                  'BGI' = 'black',
                  'UCLA'='darkred')

ggplot(data = reptile_df, 
       aes(x = log10(cn50), y = complete, size = fragmented)) +
  # Plot 'Other' and 'Iridian' first
  geom_point(data = subset(reptile_df, proj %in% c('Other', 'Iridian')), 
             aes( color = proj), 
             alpha = 0.8) +
  # Plot the other groups on top
  geom_point(data = subset(reptile_df, !proj %in% c('Other', 'Iridian')), 
             aes( color = proj), 
             alpha = 0.8) +
  scale_color_manual(values = dark_palette) +
  theme_bw() +
  xlab('Contig N50 (log10)') +
  ylab('Compleasm Complete') +
  ggtitle('Reptile Reference Genomes (Testudinata, N=47; Crocodylia, N=6; Lepidosauria, N=212)')
#ggsave(filename = 'testudinata_compleasm.svg',plot = p)
