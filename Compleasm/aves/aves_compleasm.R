library('ggplot2')
library(dplyr)
library(scales)
library(tidyr)
library(RColorBrewer)
library('svglite')

vgp_gcas<-read.table(file = '../all_gcasVGP.txt',header = F,quote = "")

aves_compleasm<-read.delim(file = 'aves_compleasmStats.txt',sep = '\t',header = T)
aves_asm<-read.table(file = 'aves_asmStats.txt',sep='\t',header=T)
aves_cnStatsFiltered<-aves_asm[,c(1,4,5,6)]
aves_merge<-merge(aves_compleasm,aves_cnStatsFiltered,by="Accession",all.x=TRUE)
# Create the dataframe with proper column assignment
aves_df <- data.frame(
  acc = aves_merge$Accession,
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
  cn50 =aves_merge$Contig_N50,
  cn90 =aves_merge$Contig_N90,
  Submitter=aves_merge$Submitter
)
for (i in c(1:length(aves_df$Submitter))){
  if(sum(vgp_gcas==aves_df$acc[i])>0){
    aves_df$proj[i]="VGP"
  }else if (aves_df$Submitter[i]=="CICHLID~X"){
    aves_df$proj[i]='Cichlid~X'
  }else if (aves_df$Submitter[i]=="B10K Consortium"){
    aves_df$proj[i]='B10K'
  }else if(aves_df$Submitter[i]=="Iridian Genomes"){
    aves_df$proj[i]='Iridian'
  }else if(aves_df$Submitter[i]=="IRIDIAN GENOMES"){
    aves_df$proj[i]='Iridian'
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
#aves_df$complete=100*aves_df$complete/3354
#aves_df$fragmented=100*aves_df$fragmented/3354
p<-ggplot(data = aves_df, 
       aes(x = log10(cn50), y = 100*(complete/3354), size = 100*(fragmented/3354))) +
  # Plot 'Other' and 'Iridian' first
  geom_point(data = subset(aves_df, proj %in% c('Other', 'Iridian')), 
             aes(color = proj), 
             alpha = 0.5,stroke=0) +
  # Add circles around 'Other' and 'Iridian'
  geom_point(data = subset(aves_df, proj %in% c('Other', 'Iridian')), 
             aes(color = proj), # Adjust the size for the outline
             shape = 21, fill = NA, stroke = 0.5) + # Create outline
  # Plot the other groups on top
  geom_point(data = subset(aves_df, !proj %in% c('Other', 'Iridian')), 
             aes(color = proj), 
             alpha = 0.5,stroke=0) +
  # Add circles around the other groups
  geom_point(data = subset(aves_df, !proj %in% c('Other', 'Iridian')), 
             aes(color = proj), # Adjust the size for the outline
             shape = 21, fill = NA, stroke = 0.5) + # Create outline
  annotate("rect", xmin = 5, xmax = 8.5, ymin = 80, ymax = 100, 
           fill = NA, color = "black", size = 0.2) +
  scale_color_manual(values = dark_palette) +
  scale_size_continuous(range = c(0, 8), 
                        limits = c(0, 55),  # Map the size scale to a fixed max of 50
                        breaks = c(10, 25, 50), 
                        labels = c("10", "25", "50")) +
  labs(color = "Submitter",  # Custom legend heading for 'color'
       size = "Fragmented (%)") +  # Custom legend heading for 'size'
  theme_bw() +xlim(2,9)+ylim(0,100)+
  xlab('Contig N50 (log10)') +
  ylab('Compleasm Complete (%)') +
  ggtitle('Aves Reference Genomes (N=1,296)')
ggsave(filename = 'aves_compleasm20241206.svg',plot = p,device = 'svg')


ggplot(data = aves_df, 
       aes(x = log10(cn90), y = 100*(complete/3354), size = 100*(fragmented/3354))) +
  # Plot 'Other' and 'Iridian' first
  geom_point(data = subset(aves_df, proj %in% c('Other', 'Iridian')), 
             aes(color = proj), 
             alpha = 0.5) +
  # Add circles around 'Other' and 'Iridian'
  geom_point(data = subset(aves_df, proj %in% c('Other', 'Iridian')), 
             aes(color = proj), # Adjust the size for the outline
             shape = 21, fill = NA, stroke = 0.5) + # Create outline
  # Plot the other groups on top
  geom_point(data = subset(aves_df, !proj %in% c('Other', 'Iridian')), 
             aes(color = proj), 
             alpha = 0.5) +
  # Add circles around the other groups
  geom_point(data = subset(aves_df, !proj %in% c('Other', 'Iridian')), 
             aes(color = proj), # Adjust the size for the outline
             shape = 21, fill = NA, stroke = 0.5) + # Create outline
  annotate("rect", xmin = 4, xmax = 8.5, ymin = 80, ymax = 100, 
           fill = NA, color = "black", size = 0.2) +
  scale_color_manual(values = dark_palette) +
  labs(color = "Submitter",  # Custom legend heading for 'color'
       size = "Fragmented (%)") +  # Custom legend heading for 'size'
  theme_bw() +xlim(2,9)+ylim(0,100)+
  xlab('Contig N90 (log10)') +
  ylab('Compleasm Complete') +
  ggtitle('Aves Reference Genomes (N=1,296)')
#ggsave(filename = 'aves_compleasm.svg',plot = p)
