library('ggplot2')
library(dplyr)
library(scales)
library(tidyr)
library(RColorBrewer)
library('svglite')

amphibia_compleasm<-read.table(file = 'amphibia_compleasmStats.txt',sep = '\t',header = T)
amphibia_asm<-read.table(file = 'amphibia_asmStats.txt',sep='\t',header=T)
amphibia_cnStatsFiltered<-amphibia_asm[,c(1,4,5,6)]
amphibia_merge<-merge(amphibia_compleasm,amphibia_cnStatsFiltered,by="Accession",all.x=TRUE)
# Create the dataframe with proper column assignment
amphibia_df <- data.frame(
  seq = amphibia_merge$Sequencing_type,
  fragmented = amphibia_merge$Fragmented_compleasm,
  frameshift = amphibia_merge$Frameshift_compleasm,
  missing = amphibia_merge$Missing_compleasm,
  complete = amphibia_merge$Complete_compleasm,
  date = amphibia_merge$Submission_Date,
  class = amphibia_merge$Class,
  superclass = amphibia_merge$Super_class,
  superorder = amphibia_merge$Super_order,
  infraclass = amphibia_merge$Infra_class,
  order = amphibia_merge$Order,
  dup = amphibia_merge$Duplicated_compleasm,
  single = amphibia_merge$Single_compleasm,
  proj='Other',
  cn50 =amphibia_merge$Contig_N50
)
for (i in c(1:length(amphibia_compleasm$Submitter))){
  if (amphibia_compleasm$Submitter[i]=="B10K Consortium"){
    amphibia_df$proj[i]='B10K'
  }else if(amphibia_compleasm$Submitter[i]=="G10K"){
    amphibia_df$proj[i]='VGP'
  }else if(amphibia_compleasm$Submitter[i]=="Iridian Genomes"){
    amphibia_df$proj[i]='Iridian'
  }else if(amphibia_compleasm$Submitter[i]=="IRIDIAN GENOMES"){
    amphibia_df$proj[i]='Iridian'
  }else if(amphibia_compleasm$Submitter[i]=="SC"){
    amphibia_df$proj[i]='VGP'
  }else if(amphibia_compleasm$Submitter[i]=="Vertebrate Genomes Project"){
    amphibia_df$proj[i]='VGP'
  }else if(amphibia_compleasm$Submitter[i]=="WELLCOME SANGER INSTITUTE"){
    amphibia_df$proj[i]='VGP'
  }else if(amphibia_compleasm$Submitter[i]=="Wellcome Sanger Institute"){
    amphibia_df$proj[i]='VGP'
  }else if(amphibia_compleasm$Submitter[i]=="The Max Planck Institute of Molecular Cell Biology and Genetics"){
    amphibia_df$proj[i]='VGP'
  }else if(amphibia_compleasm$Submitter[i]=="Zhejiang University"){
    amphibia_df$proj[i]='VGP'
  }else if(amphibia_compleasm$Submitter[i]=="BGI"){
    amphibia_df$proj[i]='BGI'
  }else if(amphibia_compleasm$Submitter[i]=="BGI-Shenzhen"){
    amphibia_df$proj[i]='BGI'
  }else if(amphibia_compleasm$Submitter[i]=="McDonnell Genome Institute - Washington University School of Medicine"){
    amphibia_df$proj[i]='McDonnell'
  }
}


ggplot(data = amphibia_df, aes(x = order, y = frameshift, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +ylim(0,500)+
  ggtitle('Compleasm genes containing frameshifts') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = amphibia_df, aes(x = order, y = complete, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  scale_fill_manual(values = colors) +
  ggtitle('Compleasm complete genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = amphibia_df, aes(x = order, y = single, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  scale_fill_manual(values = colors) +
  ggtitle('Compleasm single-copy genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = amphibia_df, aes(x = order, y = single, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  scale_fill_manual(values = colors) +
  ggtitle('Compleasm duplicate genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = amphibia_df, aes(x = order, y = fragmented, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  ggtitle('Compleasm fragmented genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##############################################################################
amphibia_dfFiltered<-amphibia_df[amphibia_df$proj %in% c('VGP','B10K','BGI','Iridian','Other'),]
amphibia_dfFiltered<-amphibia_df[amphibia_df$proj %in% c('VGP','B10K','BGI','Iridian','Other'),]
amphibia_dfFiltered$proj <- factor(amphibia_dfFiltered$proj, 
                               levels = c('VGP','BGI','B10K','Other', 'Iridian'))
dark_palette <- c('Other' = 'orange', 
                  'Iridian' = 'lightblue', 
                  'VGP' = '#16a085', 
                  'B10K' = 'darkred', 
                  'BGI' = 'black',
                  'McDonnell'='darkred')

ggplot(data = amphibia_df, 
       aes(x = log10(cn50), y = complete, size = fragmented)) +
  # Plot 'Other' and 'Iridian' first
  geom_point(data = subset(amphibia_df, proj %in% c('Other', 'Iridian')), 
             aes( color = proj), 
             alpha = 0.8) +
  # Plot the other groups on top
  geom_point(data = subset(amphibia_df, !proj %in% c('Other', 'Iridian')), 
             aes( color = proj), 
             alpha = 0.8) +
  scale_color_manual(values = dark_palette) +
  theme_bw() +
  xlab('Contig N50 (log10)') +
  ylab('Compleasm Complete') +
  ggtitle('Amphibia Reference Genomes (N=151)')
#ggsave(filename = 'amphibia_compleasm.svg',plot = p)
