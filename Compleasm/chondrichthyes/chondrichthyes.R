library('ggplot2')
library(dplyr)
library(scales)
library(tidyr)
library(RColorBrewer)
library('svglite')

chondrichthyes_compleasm<-read.table(file = 'chondrichthyes_compleasmStats.txt',sep = '\t',header = T)
chondrichthyes_asm<-read.table(file = 'chondrichthyes_asmStats.txt',sep='\t',header=T)
chondrichthyes_cnStatsFiltered<-chondrichthyes_asm[,c(1,4,5,6)]
chondrichthyes_merge<-merge(chondrichthyes_compleasm,chondrichthyes_cnStatsFiltered,by="Accession",all.x=TRUE)
# Create the dataframe with proper column assignment
chondrichthyes_df <- data.frame(
  seq = chondrichthyes_merge$Sequencing_type,
  fragmented = chondrichthyes_merge$Fragmented_compleasm,
  frameshift = chondrichthyes_merge$Frameshift_compleasm,
  missing = chondrichthyes_merge$Missing_compleasm,
  complete = chondrichthyes_merge$Complete_compleasm,
  date = chondrichthyes_merge$Submission_Date,
  class = chondrichthyes_merge$Class,
  superclass = chondrichthyes_merge$Super_class,
  superorder = chondrichthyes_merge$Super_order,
  infraclass = chondrichthyes_merge$Infra_class,
  order = chondrichthyes_merge$Order,
  dup = chondrichthyes_merge$Duplicated_compleasm,
  single = chondrichthyes_merge$Single_compleasm,
  proj='Other',
  cn50 =chondrichthyes_merge$Contig_N50
)
for (i in c(1:length(chondrichthyes_compleasm$Submitter))){
  if (chondrichthyes_compleasm$Submitter[i]=="B10K Consortium"){
    chondrichthyes_df$proj[i]='B10K'
  }else if(chondrichthyes_compleasm$Submitter[i]=="G10K"){
    chondrichthyes_df$proj[i]='VGP'
  }else if(chondrichthyes_compleasm$Submitter[i]=="Iridian Genomes"){
    chondrichthyes_df$proj[i]='Iridian'
  }else if(chondrichthyes_compleasm$Submitter[i]=="IRIDIAN GENOMES"){
    chondrichthyes_df$proj[i]='Iridian'
  }else if(chondrichthyes_compleasm$Submitter[i]=="Vertebrate Genomes Project"){
    chondrichthyes_df$proj[i]='VGP'
  }else if(chondrichthyes_compleasm$Submitter[i]=="WELLCOME SANGER INSTITUTE"){
    chondrichthyes_df$proj[i]='VGP'
  }else if(chondrichthyes_compleasm$Submitter[i]=="Wellcome Sanger Institute"){
    chondrichthyes_df$proj[i]='VGP'
  }else if(chondrichthyes_compleasm$Submitter[i]=="The Max Planck Institute of Molecular Cell Biology and Genetics"){
    chondrichthyes_df$proj[i]='VGP'
  }else if(chondrichthyes_compleasm$Submitter[i]=="Zhejiang University"){
    chondrichthyes_df$proj[i]='VGP'
  }else if(chondrichthyes_compleasm$Submitter[i]=="BGI"){
    chondrichthyes_df$proj[i]='BGI'
  }else if(chondrichthyes_compleasm$Submitter[i]=="BGI-Shenzhen"){
    chondrichthyes_df$proj[i]='BGI'
  }else if(chondrichthyes_compleasm$Submitter[i]=="McDonnell Genome Institute - Washington University School of Medicine"){
    chondrichthyes_df$proj[i]='McDonnell'
  }
}


ggplot(data = chondrichthyes_df, aes(x = order, y = frameshift, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +ylim(0,500)+
  ggtitle('Compleasm genes containing frameshifts') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = chondrichthyes_df, aes(x = order, y = complete, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  scale_fill_manual(values = colors) +
  ggtitle('Compleasm complete genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = chondrichthyes_df, aes(x = order, y = single, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  scale_fill_manual(values = colors) +
  ggtitle('Compleasm single-copy genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = chondrichthyes_df, aes(x = order, y = single, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  scale_fill_manual(values = colors) +
  ggtitle('Compleasm duplicate genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = chondrichthyes_df, aes(x = order, y = fragmented, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  ggtitle('Compleasm fragmented genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##############################################################################
chondrichthyes_dfFiltered<-chondrichthyes_df[chondrichthyes_df$proj %in% c('VGP','B10K','BGI','Iridian','Other'),]
chondrichthyes_dfFiltered<-chondrichthyes_df[chondrichthyes_df$proj %in% c('VGP','B10K','BGI','Iridian','Other'),]
chondrichthyes_dfFiltered$proj <- factor(chondrichthyes_dfFiltered$proj, 
                               levels = c('VGP','BGI','B10K','Other', 'Iridian'))
dark_palette <- c('Other' = 'orange', 
                  'Iridian' = 'lightblue', 
                  'VGP' = '#16a085', 
                  'B10K' = 'darkred', 
                  'BGI' = 'black',
                  'McDonnell'='darkred')

ggplot(data = chondrichthyes_df, 
       aes(x = log10(cn50), y = complete, size = fragmented)) +
  # Plot 'Other' and 'Iridian' first
  geom_point(data = subset(chondrichthyes_df, proj %in% c('Other', 'Iridian')), 
             aes( color = proj), 
             alpha = 0.8) +
  # Plot the other groups on top
  geom_point(data = subset(chondrichthyes_df, !proj %in% c('Other', 'Iridian')), 
             aes( color = proj), 
             alpha = 0.8) +
  scale_color_manual(values = dark_palette) +
  theme_bw() +
  xlab('Contig N50 (log10)') +
  ylab('Compleasm Complete') +
  ggtitle('chondrichthyes Reference Genomes (N=47)')
#ggsave(filename = 'chondrichthyes_compleasm.svg',plot = p)
