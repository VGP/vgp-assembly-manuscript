library('ggplot2')
library(dplyr)
library(scales)
library(tidyr)
library(RColorBrewer)
library('svglite')

vgp_gcas<-read.table(file = '../all_gcasVGP.txt',header = F,quote = "")

testudinata_compleasm<-read.table(file = 'testudinata_compleasmStats.txt',sep = '\t',header = T)
testudinata_asm<-read.table(file = 'testudinata_asmStats.txt',sep='\t',header=T)
testudinata_cnStatsFiltered<-testudinata_asm[,c(1,4,5,6)]
testudinata_merge<-merge(testudinata_compleasm,testudinata_cnStatsFiltered,by="Accession",all.x=TRUE)
# Create the dataframe with proper column assignment
testudinata_df <- data.frame(
  acc = testudinata_merge$Accession,
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
  cn50 =testudinata_merge$Contig_N50,
  cn90 =testudinata_merge$Contig_N90,
  Submitter=testudinata_merge$Submitter
)
crocodylia_compleasm<-read.table(file = 'crocodylia_compleasmStats.txt',sep = '\t',header = T)
crocodylia_asm<-read.table(file = 'crocodylia_asmStats.txt',sep='\t',header=T)
crocodylia_cnStatsFiltered<-crocodylia_asm[,c(1,4,5,6)]
crocodylia_merge<-merge(crocodylia_compleasm,crocodylia_cnStatsFiltered,by="Accession",all.x=TRUE)
# Create the dataframe with proper column assignment
crocodylia_df <- data.frame(
  acc = crocodylia_merge$Accession,
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
  cn50 =crocodylia_merge$Contig_N50,
  cn90 =crocodylia_merge$Contig_N90,
  Submitter=crocodylia_merge$Submitter
)
lepidosauria_compleasm<-read.table(file = 'lepidosauria_compleasmStats.txt',sep = '\t',header = T)
lepidosauria_asm<-read.table(file = 'lepidosauria_asmStats.txt',sep='\t',header=T)
lepidosauria_cnStatsFiltered<-lepidosauria_asm[,c(1,4,5,6)]
lepidosauria_merge<-merge(lepidosauria_compleasm,lepidosauria_cnStatsFiltered,by="Accession",all.x=TRUE)
# Create the dataframe with proper column assignment
lepidosauria_df <- data.frame(
  acc = lepidosauria_merge$Accession,
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
  cn50 =lepidosauria_merge$Contig_N50,
  cn90 =lepidosauria_merge$Contig_N90,
  Submitter=lepidosauria_merge$Submitter
)

reptile_df<-rbind(testudinata_df,crocodylia_df,lepidosauria_df)
reptile_submitter<-c(testudinata_compleasm$Submitter,crocodylia_compleasm$Submitter,lepidosauria_compleasm$Submitter)
for (i in c(1:length(reptile_submitter))){
  if(sum(vgp_gcas==reptile_df$acc[i])>0){
    reptile_df$proj[i]="VGP"
  }else if (reptile_submitter[i]=="B10K Consortium"){
    reptile_df$proj[i]='B10K'
  }else if(reptile_submitter[i]=="Iridian Genomes"){
    reptile_df$proj[i]='Iridian'
  }else if(reptile_submitter[i]=="IRIDIAN GENOMES"){
    reptile_df$proj[i]='Iridian'
  }else if(reptile_submitter[i]=="BGI"){
    reptile_df$proj[i]='BGI'
  }else if(reptile_submitter[i]=="BGI-Shenzhen"){
    reptile_df$proj[i]='BGI'
  }
}


dark_palette <- c('Cichlid~X' = 'darkred',
                  'Iridian' = 'darkgreen', 
                  'Other' = 'lightblue',
                  'B10K'= 'red',
                  'DNAZoo'= 'yellow',
                  'Broad' = 'blue',
                  'VGP' = 'purple'
)

p<-ggplot(data = reptile_df, 
       aes(x = log10(cn50), y = 100*(complete/3354), size = 100*(fragmented/3354))) +
  # Plot 'Other' and 'Iridian' first
  geom_point(data = subset(reptile_df, proj %in% c('Other', 'Iridian')), 
             aes(color = proj), 
             alpha = 0.5) +
  # Add circles around 'Other' and 'Iridian'
  geom_point(data = subset(reptile_df, proj %in% c('Other', 'Iridian')), 
             aes(color = proj), # Adjust the size for the outline
             shape = 21, fill = NA, stroke = 0.5) + # Create outline
  # Plot the other groups on top
  geom_point(data = subset(reptile_df, !proj %in% c('Other', 'Iridian')), 
             aes(color = proj), 
             alpha = 0.5) +
  # Add circles around the other groups
  geom_point(data = subset(reptile_df, !proj %in% c('Other', 'Iridian')), 
             aes(color = proj), # Adjust the size for the outline
             shape = 21, fill = NA, stroke = 0.5) + # Create outline
  annotate("rect", xmin = 5, xmax = 8.5, ymin = 80, ymax = 100, 
           fill = NA, color = "black", linewidth = 0.2) +
  scale_color_manual(values = dark_palette) +
  theme_bw() +xlim(2,9)+ylim(0,100)+
  scale_size_continuous(range = c(0, 8), 
                        limits = c(0, 55),  # Map the size scale to a fixed max of 50
                        breaks = c(10, 25, 50), 
                        labels = c("10", "25", "50")) +
  labs(color = "Submitter",  # Custom legend heading for 'color'
       size = "Fragmented (%)") +  # Custom legend heading for 'size'
  xlab('Contig N50 (log10)') +
  ylab('Compleasm Complete (%)') +
  ggtitle('Reptile Reference Genomes (N=264)')
ggsave(filename = 'reptiles_compleasm20241206.svg',plot = p,device = 'svg')


ggplot(data = reptile_df, 
       aes(x = log10(cn90), y = 100*(complete/3354), size = 100*(fragmented/3354))) +
  # Plot 'Other' and 'Iridian' first
  geom_point(data = subset(reptile_df, proj %in% c('Other', 'Iridian')), 
             aes(color = proj), 
             alpha = 0.5) +
  # Add circles around 'Other' and 'Iridian'
  geom_point(data = subset(reptile_df, proj %in% c('Other', 'Iridian')), 
             aes(color = proj), # Adjust the size for the outline
             shape = 21, fill = NA, stroke = 0.5) + # Create outline
  # Plot the other groups on top
  geom_point(data = subset(reptile_df, !proj %in% c('Other', 'Iridian')), 
             aes(color = proj), 
             alpha = 0.5) +
  # Add circles around the other groups
  geom_point(data = subset(reptile_df, !proj %in% c('Other', 'Iridian')), 
             aes(color = proj), # Adjust the size for the outline
             shape = 21, fill = NA, stroke = 0.5) + # Create outline
  annotate("rect", xmin = 4, xmax = 8.5, ymin = 80, ymax = 100, 
           fill = NA, color = "black", size = 0.2) +
  scale_color_manual(values = dark_palette) +
  theme_bw() +xlim(2,9)+ylim(0,100)+
  labs(color = "Submitter",  # Custom legend heading for 'color'
       size = "Fragmented (%)") +  # Custom legend heading for 'size'
  xlab('Contig N90 (log10)') +
  ylab('Compleasm Complete') +
  ggtitle('Reptile Reference Genomes (N=264)')
#ggsave(filename = 'crocodylia_compleasm.svg',plot = p)
