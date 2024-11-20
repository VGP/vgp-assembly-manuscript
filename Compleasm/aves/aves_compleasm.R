library('ggplot2')
library(dplyr)
library(scales)
library(tidyr)
library(RColorBrewer)
library('svglite')

aves_compleasm<-read.table(file = 'aves_compleasmStats.txt',sep = '\t',header = T)
aves_asm<-read.table(file = 'aves_asmStats.txt',sep='\t',header=T)
aves_cnStatsFiltered<-aves_asm[,c(1,4,5,6)]
aves_merge<-merge(aves_compleasm,aves_cnStatsFiltered,by="Accession",all.x=TRUE)
# Create the dataframe with proper column assignment
aves_df <- data.frame(
  seq = aves_merge$Sequencing_type,
  fragmented = aves_merge$Fragmented_compleasm,
  frameshift = aves_merge$Frameshift_compleasm,
  missing = aves_merge$Missing_compleasm,
  complete = aves_merge$Complete_compleasm,
  date = aves_merge$Submission_Date,
  class = aves_merge$Class,
  superclass = aves_merge$Super_class,
  superorder = aves_merge$Super_order,
  infraclass = aves_merge$Infra_class,
  order = aves_merge$Order,
  dup = aves_merge$Duplicated_compleasm,
  single = aves_merge$Single_compleasm,
  proj='Other',
  cn50 =aves_merge$Contig_N50
)
for (i in c(1:length(aves_compleasm$Submitter))){
  if (aves_compleasm$Submitter[i]=="B10K Consortium"){
    aves_df$proj[i]='B10K'
  }else if(aves_compleasm$Submitter[i]=="G10K"){
    aves_df$proj[i]='VGP'
  }else if(aves_compleasm$Submitter[i]=="Iridian Genomes"){
    aves_df$proj[i]='Iridian'
  }else if(aves_compleasm$Submitter[i]=="IRIDIAN GENOMES"){
    aves_df$proj[i]='Iridian'
  }else if(aves_compleasm$Submitter[i]=="Vertebrate Genomes Project"){
    aves_df$proj[i]='VGP'
  }else if(aves_compleasm$Submitter[i]=="WELLCOME SANGER INSTITUTE"){
    aves_df$proj[i]='VGP'
  }else if(aves_compleasm$Submitter[i]=="Wellcome Sanger Institute"){
    aves_df$proj[i]='VGP'
  }else if(aves_compleasm$Submitter[i]=="The Max Planck Institute of Molecular Cell Biology and Genetics"){
    aves_df$proj[i]='VGP'
  }else if(aves_compleasm$Submitter[i]=="Zhejiang University"){
    aves_df$proj[i]='VGP'
  }else if(aves_compleasm$Submitter[i]=="BGI"){
    aves_df$proj[i]='BGI'
  }else if(aves_compleasm$Submitter[i]=="BGI-Shenzhen"){
    aves_df$proj[i]='BGI'
  }
}


ggplot(data = aves_df, aes(x = order, y = frameshift, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +ylim(0,500)+
  ggtitle('Compleasm genes containing frameshifts') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = aves_df, aes(x = order, y = complete, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  scale_fill_manual(values = colors) +
  ggtitle('Compleasm complete genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = aves_df, aes(x = order, y = single, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  scale_fill_manual(values = colors) +
  ggtitle('Compleasm single-copy genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = aves_df, aes(x = order, y = single, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  scale_fill_manual(values = colors) +
  ggtitle('Compleasm duplicate genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = aves_df, aes(x = order, y = fragmented, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  ggtitle('Compleasm fragmented genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##############################################################################
aves_dfFiltered<-aves_df[aves_df$proj %in% c('VGP','B10K','BGI','Iridian','Other'),]
aves_dfFiltered<-aves_df[aves_df$proj %in% c('VGP','B10K','BGI','Iridian','Other'),]
aves_dfFiltered$proj <- factor(aves_dfFiltered$proj, 
                               levels = c('VGP','BGI','B10K','Other', 'Iridian'))
dark_palette <- c('Other' = 'orange', 
                  'Iridian' = 'lightblue', 
                  'VGP' = '#16a085', 
                  'B10K' = 'darkred', 
                  'BGI' = 'black')

ggplot(data = aves_dfFiltered, 
       aes(x = log10(cn50), y = complete, size = fragmented)) +
  # Plot 'Other' and 'Iridian' first
  geom_point(data = subset(aves_dfFiltered, proj %in% c('Other', 'Iridian')), 
             aes( color = proj), 
             alpha = 0.8) +
  # Plot the other groups on top
  geom_point(data = subset(aves_dfFiltered, !proj %in% c('Other', 'Iridian')), 
             aes( color = proj), 
             alpha = 0.8) +
  scale_color_manual(values = dark_palette) +
  theme_bw() +
  xlab('Contig N50 (log10)') +
  ylab('Compleasm Complete') +
  ggtitle('Aves Reference Genomes (N=1,296)')
ggsave(filename = 'aves_compleasm.svg',plot = p)
