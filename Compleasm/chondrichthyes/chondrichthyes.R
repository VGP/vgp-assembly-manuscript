library('ggplot2')
library(dplyr)
library(scales)
library(tidyr)
library(RColorBrewer)
library('svglite')

vgp_gcas<-read.table(file = '../all_gcasVGP.txt',header = F,quote = "")

chondrichthyes_compleasm<-read.delim(file = 'chondrichthyes_compleasmStats.txt',sep = '\t',header = T)
chondrichthyes_asm<-read.table(file = 'chondrichthyes_asmStats.txt',sep='\t',header=T)
chondrichthyes_cnStatsFiltered<-chondrichthyes_asm[,c(1,4,5,6)]
chondrichthyes_merge<-merge(chondrichthyes_compleasm,chondrichthyes_cnStatsFiltered,by="Accession",all.x=TRUE)
# Create the dataframe with proper column assignment
chondrichthyes_df <- data.frame(
  acc = chondrichthyes_merge$Accession,
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
  cn50 =chondrichthyes_merge$Contig_N50,
  cn90 =chondrichthyes_merge$Contig_N90,
  Submitter=chondrichthyes_merge$Submitter
)
for (i in c(1:length(chondrichthyes_df$Submitter))){
  if(sum(vgp_gcas==chondrichthyes_df$acc[i])>0){
    chondrichthyes_df$proj[i]="VGP"
  }else if (chondrichthyes_df$Submitter[i]=="CICHLID~X"){
    chondrichthyes_df$proj[i]='Cichlid~X'
  }else if(chondrichthyes_df$Submitter[i]=="Iridian Genomes"){
    chondrichthyes_df$proj[i]='Iridian'
  }else if(chondrichthyes_df$Submitter[i]=="IRIDIAN GENOMES"){
    chondrichthyes_df$proj[i]='Iridian'
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
#chondrichthyes_df$complete=100*chondrichthyes_df$complete/3354
#chondrichthyes_df$fragmented=100*chondrichthyes_df$fragmented/3354
p<-ggplot(data = chondrichthyes_df, 
       aes(x = log10(cn50), y = 100*(complete/3354), size = 100*(fragmented/3354))) +
  # Plot 'Other' and 'Iridian' first
  geom_point(data = subset(chondrichthyes_df, proj %in% c('Other', 'Iridian')), 
             aes(color = proj), 
             alpha = 0.5,stroke=0) +
  # Add circles around 'Other' and 'Iridian'
  geom_point(data = subset(chondrichthyes_df, proj %in% c('Other', 'Iridian')), 
             aes(color = proj), # Adjust the size for the outline
             shape = 21, fill = NA, stroke = 0.5) + # Create outline
  # Plot the other groups on top
  geom_point(data = subset(chondrichthyes_df, !proj %in% c('Other', 'Iridian')), 
             aes(color = proj), 
             alpha = 0.5,stroke=0) +
  # Add circles around the other groups
  geom_point(data = subset(chondrichthyes_df, !proj %in% c('Other', 'Iridian')), 
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
  ggtitle('Chondrichthyes Reference Genomes (N=41)')
ggsave(filename = 'chondrichthyes_compleasm20241206.svg',plot = p,device = 'svg')


ggplot(data = chondrichthyes_df, 
       aes(x = log10(cn90), y = 100*(complete/3354), size = 100*(fragmented/3354))) +
  # Plot 'Other' and 'Iridian' first
  geom_point(data = subset(chondrichthyes_df, proj %in% c('Other', 'Iridian')), 
             aes(color = proj), 
             alpha = 0.5) +
  # Add circles around 'Other' and 'Iridian'
  geom_point(data = subset(chondrichthyes_df, proj %in% c('Other', 'Iridian')), 
             aes(color = proj), # Adjust the size for the outline
             shape = 21, fill = NA, stroke = 0.5) + # Create outline
  # Plot the other groups on top
  geom_point(data = subset(chondrichthyes_df, !proj %in% c('Other', 'Iridian')), 
             aes(color = proj), 
             alpha = 0.5) +
  # Add circles around the other groups
  geom_point(data = subset(chondrichthyes_df, !proj %in% c('Other', 'Iridian')), 
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
  ggtitle('Chondrichthyes Reference Genomes (N=41)')
#ggsave(filename = 'chondrichthyes_compleasm.svg',plot = p)
