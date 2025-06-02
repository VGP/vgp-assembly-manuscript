library('ggplot2')
library(dplyr)
library(scales)
library(tidyr)
library(RColorBrewer)
library('svglite')

vgp_gcas<-read.table(file = '../all_gcasVGP.txt',header = F,quote = "")

amphibia_compleasm<-read.delim(file = 'amphibia_compleasmStats.txt',sep = '\t',header = T)
amphibia_asm<-read.table(file = 'amphibia_asmStats.txt',sep='\t',header=T)
amphibia_cnStatsFiltered<-amphibia_asm[,c(1,4,5,6)]
amphibia_merge<-merge(amphibia_compleasm,amphibia_cnStatsFiltered,by="Accession",all.x=TRUE)
# Create the dataframe with proper column assignment
amphibia_df <- data.frame(
  acc = amphibia_merge$Accession,
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
  cn50 =amphibia_merge$Contig_N50,
  cn90 =amphibia_merge$Contig_N90,
  Submitter=amphibia_merge$Submitter
)
for (i in c(1:length(amphibia_df$Submitter))){
  if(sum(vgp_gcas==amphibia_df$acc[i])>0){
    amphibia_df$proj[i]="VGP"
  }else if (amphibia_df$Submitter[i]=="CICHLID~X"){
    amphibia_df$proj[i]='Cichlid~X'
  }else if(amphibia_df$Submitter[i]=="Iridian Genomes"){
    amphibia_df$proj[i]='Iridian'
  }else if(amphibia_df$Submitter[i]=="IRIDIAN GENOMES"){
    amphibia_df$proj[i]='Iridian'
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
#amphibia_df$complete=100*amphibia_df$complete/3354
#amphibia_df$fragmented=100*amphibia_df$fragmented/3354
p<-ggplot(data = amphibia_df, 
       aes(x = log10(cn50), y = 100*(complete/3354), size = 100*(fragmented/3354))) +
  # Plot 'Other' and 'Iridian' first
  geom_point(data = subset(amphibia_df, proj %in% c('Other', 'Iridian')), 
             aes(color = proj), 
             alpha = 0.5,stroke=0) +
  # Add circles around 'Other' and 'Iridian'
  geom_point(data = subset(amphibia_df, proj %in% c('Other', 'Iridian')), 
             aes(color = proj), # Adjust the size for the outline
             shape = 21, fill = NA, stroke = 0.5) + # Create outline
  # Plot the other groups on top
  geom_point(data = subset(amphibia_df, !proj %in% c('Other', 'Iridian')), 
             aes(color = proj), 
             alpha = 0.5,stroke=0) +
  # Add circles around the other groups
  geom_point(data = subset(amphibia_df, !proj %in% c('Other', 'Iridian')), 
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
  ylab('Compleasm Complete') +
  ggtitle('Amphibia Reference Genomes (N=151)')
ggsave(filename = 'amphibia_compleasm20241206.svg',plot = p,device = 'svg')

ggplot(data = amphibia_df, 
       aes(x = log10(cn90), y = 100*(complete/3354), size = 100*(fragmented/3354))) +
  # Plot 'Other' and 'Iridian' first
  geom_point(data = subset(amphibia_df, proj %in% c('Other', 'Iridian')), 
             aes(color = proj), 
             alpha = 0.5) +
  # Add circles around 'Other' and 'Iridian'
  geom_point(data = subset(amphibia_df, proj %in% c('Other', 'Iridian')), 
             aes(color = proj), # Adjust the size for the outline
             shape = 21, fill = NA, stroke = 0.5) + # Create outline
  # Plot the other groups on top
  geom_point(data = subset(amphibia_df, !proj %in% c('Other', 'Iridian')), 
             aes(color = proj), 
             alpha = 0.5) +
  # Add circles around the other groups
  geom_point(data = subset(amphibia_df, !proj %in% c('Other', 'Iridian')), 
             aes(color = proj), # Adjust the size for the outline
             shape = 21, fill = NA, stroke = 0.5) + # Create outline
  annotate("rect", xmin = 4, xmax = 8.5, ymin = 80, ymax = 100, 
           fill = NA, color = "black", size = 0.2) +
  labs(color = "Submitter",  # Custom legend heading for 'color'
       size = "Fragmented (%)") +  # Custom legend heading for 'size'
  scale_color_manual(values = dark_palette) +
  theme_bw() +xlim(2,9)+ylim(0,100)+
  xlab('Contig N90 (log10)') +
  ylab('Compleasm Complete') +
  ggtitle('amphibia Reference Genomes (N=151)')
#ggsave(filename = 'amphibia_compleasm.svg',plot = p)
