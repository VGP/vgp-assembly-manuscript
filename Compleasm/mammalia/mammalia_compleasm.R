library('ggplot2')
library(dplyr)
library(scales)
library(tidyr)
library(RColorBrewer)
library('svglite')

vgp_gcas<-read.table(file = '../all_gcasVGP.txt',header = F,quote = "")

mammalia_compleasm<-read.delim(file = 'mammalia_compleasmStats.txt',sep = '\t',header = T)
mammalia_asm<-read.table(file = 'mammalia_asmStats.txt',sep='\t',header=T)
mammalia_cnStatsFiltered<-mammalia_asm[,c(1,4,5,6)]
mammalia_merge<-merge(mammalia_compleasm,mammalia_cnStatsFiltered,by="Accession",all.x=TRUE)
# Create the dataframe with proper column assignment
mammalia_df <- data.frame(
  acc = mammalia_merge$Accession,
  seq = mammalia_merge$Sequencing_type,
  fragmented = mammalia_merge$Fragmented_compleasm,
  frameshift = mammalia_merge$Frameshift_compleasm,
  missing = mammalia_merge$Missing_compleasm,
  complete = mammalia_merge$Complete_compleasm,
  date = mammalia_merge$Submission_Date,
  class = mammalia_merge$Class,
  superclass = mammalia_merge$Super_class,
  superorder = mammalia_merge$Super_order,
  infraclass = mammalia_merge$Infra_class,
  order = mammalia_merge$Order,
  dup = mammalia_merge$Duplicated_compleasm,
  single = mammalia_merge$Single_compleasm,
  proj='Other',
  cn50 =mammalia_merge$Contig_N50,
  cn90 =mammalia_merge$Contig_N90,
  Submitter=mammalia_merge$Submitter
)
for (i in c(1:length(mammalia_df$Submitter))){
  if(sum(vgp_gcas==mammalia_df$acc[i])>0){
    mammalia_df$proj[i]="VGP"
  }else if (mammalia_df$Submitter[i]=="CICHLID~X"){
    mammalia_df$proj[i]='Cichlid~X'
  }else if (mammalia_df$Submitter[i]=="B10K Consortium"){
    mammalia_df$proj[i]='B10K'
  }else if (mammalia_df$Submitter[i]=="DNA Zoo"){
    mammalia_df$proj[i]='DNAZoo'
  }else if (mammalia_df$Submitter[i]=="Broad Institute"){
    mammalia_df$proj[i]='Broad'
  }else if(mammalia_df$Submitter[i]=="Iridian Genomes"){
    mammalia_df$proj[i]='Iridian'
  }else if(mammalia_df$Submitter[i]=="IRIDIAN GENOMES"){
    mammalia_df$proj[i]='Iridian'
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
#mammalia_df$complete=100*mammalia_df$complete/3354
#mammalia_df$fragmented=100*mammalia_df$fragmented/3354
p<-ggplot(data = mammalia_df, 
       aes(x = log10(cn50), y = 100*(complete/3354), size = 100*(fragmented/3354))) +
  # Plot 'Other' and 'Iridian' first
  geom_point(data = subset(mammalia_df, proj %in% c('Other', 'Iridian')), 
             aes(color = proj), 
             alpha = 0.5,stroke=0) +
  # Add circles around 'Other' and 'Iridian'
  geom_point(data = subset(mammalia_df, proj %in% c('Other', 'Iridian')), 
             aes(color = proj), # Adjust the size for the outline
             shape = 21, fill = NA, stroke = 0.5) + # Create outline
  # Plot the other groups on top
  geom_point(data = subset(mammalia_df, !proj %in% c('Other', 'Iridian')), 
             aes(color = proj), 
             alpha = 0.5,stroke=0) +
  # Add circles around the other groups
  geom_point(data = subset(mammalia_df, !proj %in% c('Other', 'Iridian')), 
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
  ggtitle('Mammalia Reference Genomes (N=1,296)')
ggsave(filename = 'mammalia_compleasm20241206.svg',plot = p,device = 'svg')


ggplot(data = mammalia_df, 
       aes(x = log10(cn90), y = 100*(complete/3354), size = 100*(fragmented/3354))) +
  # Plot 'Other' and 'Iridian' first
  geom_point(data = subset(mammalia_df, proj %in% c('Other', 'Iridian')), 
             aes(color = proj), 
             alpha = 0.5) +
  # Add circles around 'Other' and 'Iridian'
  geom_point(data = subset(mammalia_df, proj %in% c('Other', 'Iridian')), 
             aes(color = proj), # Adjust the size for the outline
             shape = 21, fill = NA, stroke = 0.5) + # Create outline
  # Plot the other groups on top
  geom_point(data = subset(mammalia_df, !proj %in% c('Other', 'Iridian')), 
             aes(color = proj), 
             alpha = 0.5) +
  # Add circles around the other groups
  geom_point(data = subset(mammalia_df, !proj %in% c('Other', 'Iridian')), 
             aes(color = proj), # Adjust the size for the outline
             shape = 21, fill = NA, stroke = 0.5) + # Create outline
  annotate("rect", xmin = 4, xmax = 8.5, ymin = 80, ymax = 100, 
           fill = NA, color = "black", size = 0.2) +
  scale_color_manual(values = dark_palette) +
  labs(color = "Submitter",  # Custom legend heading for 'color'
       size = "Fragmented (%)") +  # Custom legend heading for 'size'
  theme_bw() +xlim(2,9)+ylim(0,100)+
  xlab('Contig N90 (log10)') +
  ylab('Compleasm Complete (%)') +
  ggtitle('Mammalia Reference Genomes (N=1,296)')
#ggsave(filename = 'mammalia_compleasm.svg',plot = p)
