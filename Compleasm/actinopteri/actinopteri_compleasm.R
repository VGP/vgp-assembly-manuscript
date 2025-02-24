library('ggplot2')
library(dplyr)
library(scales)
library(tidyr)
library(RColorBrewer)
library('svglite')

vgp_gcas<-read.table(file = '../all_gcasVGP.txt',header = F,quote = "")

actinopteri_compleasm<-read.delim(file = 'actinopteri_compleasmStats.txt',sep = '\t',header = T)
actinopteri_asm<-read.table(file = 'actinopteri_asmStats.txt',sep='\t',header=T)
actinopteri_cnStatsFiltered<-actinopteri_asm[,c(1,4,5,6)]
actinopteri_merge<-merge(actinopteri_compleasm,actinopteri_cnStatsFiltered,by="Accession",all.x=TRUE)
# Create the dataframe with proper column assignment
actinopteri_df <- data.frame(
  acc = actinopteri_merge$Accession,
  seq = actinopteri_merge$Sequencing_type,
  fragmented = actinopteri_merge$Fragmented_compleasm,
  frameshift = actinopteri_merge$Frameshift_compleasm,
  missing = actinopteri_merge$Missing_compleasm,
  complete = actinopteri_merge$Complete_compleasm,
  date = actinopteri_merge$Submission_Date,
  class = actinopteri_merge$Class,
  superclass = actinopteri_merge$Super_class,
  superorder = actinopteri_merge$Super_order,
  infraclass = actinopteri_merge$Infra_class,
  order = actinopteri_merge$Order,
  dup = actinopteri_merge$Duplicated_compleasm,
  single = actinopteri_merge$Single_compleasm,
  proj='Other',
  cn50 =actinopteri_merge$Contig_N50,
  cn90 =actinopteri_merge$Contig_N90,
  Submitter=actinopteri_merge$Submitter
)
for (i in c(1:length(actinopteri_df$Submitter))){
  if(sum(vgp_gcas==actinopteri_df$acc[i])>0){
    actinopteri_df$proj[i]="VGP"
  }else if (actinopteri_df$Submitter[i]=="CICHLID~X"){
    actinopteri_df$proj[i]='Cichlid~X'
  }else if(actinopteri_df$Submitter[i]=="Iridian Genomes"){
    actinopteri_df$proj[i]='Iridian'
  }else if(actinopteri_df$Submitter[i]=="IRIDIAN GENOMES"){
    actinopteri_df$proj[i]='Iridian'
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
#actinopteri_df$complete=100*actinopteri_df$complete/3354
#actinopteri_df$fragmented=100*actinopteri_df$fragmented/3354
p<-ggplot(data = actinopteri_df, 
          aes(x = log10(cn50), y = 100*(complete/3354), size = 100*(fragmented/3354))) +
  # Plot 'Other' and 'Iridian' first
  geom_point(data = subset(actinopteri_df, proj %in% c('Other', 'Iridian')), 
             aes(color = proj), 
             alpha = 0.5,stroke=0) +
  # Add circles around 'Other' and 'Iridian'
  geom_point(data = subset(actinopteri_df, proj %in% c('Other', 'Iridian')), 
             aes(color = proj), # Adjust the size for the outline
             shape = 21, fill = NA, stroke = 0.5) + # Create outline
  # Plot the other groups on top
  geom_point(data = subset(actinopteri_df, !proj %in% c('Other', 'Iridian')), 
             aes(color = proj), 
             alpha = 0.5,stroke=0) +
  # Add circles around the other groups
  geom_point(data = subset(actinopteri_df, !proj %in% c('Other', 'Iridian')), 
             aes(color = proj), # Adjust the size for the outline
             shape = 21, fill = NA, stroke = 0.5) + # Create outline
  annotate("rect", xmin = 5, xmax = 8.5, ymin = 80, ymax = 100, 
           fill = NA, color = "black", size = 0.2) +
  labs(color = "Submitter",  # Custom legend heading for 'color'
       size = "Fragmented (%)") +  # Custom legend heading for 'size'
  scale_color_manual(values = dark_palette) +
  theme_bw() +xlim(2,9)+ylim(0,100)+
  scale_size_continuous(range = c(0, 8), 
                        limits = c(0, 55),  # Map the size scale to a fixed max of 50
                        breaks = c(10, 25, 50), 
                        labels = c("10", "25", "50")) +
  xlab('Contig N50 (log10)') +
  ylab('Compleasm Complete (%)') +
  ggtitle('Actinopteri Reference Genomes (N=1,634)')
ggsave(filename = 'actinopteri_compleasm20241206.svg',plot = p,device = 'svg')


ggplot(data = actinopteri_df, 
       aes(x = log10(cn90), y = 100*(complete/3354), size = 100*(fragmented/3354))) +
  # Plot 'Other' and 'Iridian' first
  geom_point(data = subset(actinopteri_df, proj %in% c('Other', 'Iridian')), 
             aes(color = proj), 
             alpha = 0.5) +
  # Add circles around 'Other' and 'Iridian'
  geom_point(data = subset(actinopteri_df, proj %in% c('Other', 'Iridian')), 
             aes(color = proj), # Adjust the size for the outline
             shape = 21, fill = NA, stroke = 0.5) + # Create outline
  # Plot the other groups on top
  geom_point(data = subset(actinopteri_df, !proj %in% c('Other', 'Iridian')), 
             aes(color = proj), 
             alpha = 0.5) +
  # Add circles around the other groups
  geom_point(data = subset(actinopteri_df, !proj %in% c('Other', 'Iridian')), 
             aes(color = proj), # Adjust the size for the outline
             shape = 21, fill = NA, stroke = 0.5) + # Create outline
  annotate("rect", xmin = 4, xmax = 8.5, ymin = 80, ymax = 100, 
           fill = NA, color = "black", size = 0.2) +
  scale_color_manual(values = dark_palette) +
  theme_bw() +xlim(2,9)+ylim(0,100)+
  labs(color = "Submitter",  # Custom legend heading for 'color'
       size = "Fragmented (%)") +  # Custom legend heading for 'size'
  xlab('Contig N90 (log10)') +
  ylab('Compleasm Complete (%)') +
  ggtitle('Actinopteri Reference Genomes (N=1,634)')
#ggsave(filename = 'actinopteri_compleasm.svg',plot = p)
