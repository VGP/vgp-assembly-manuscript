
library('ggplot2')
library(dplyr)
library(scales)
library(tidyr)
library(stringr)

genbank_dataDF<-read.table(file = 'genbank_genomesData.tsv',header = T,sep = '\t',colClasses = c(accession='character'))
my_colors <- c(
  "#B0B0B0","#0072B2","" # highlight
)

#Show genome completeness vs contig n50 for vertebrate genomes on genbank with those from the VGP datafreeze coloured in blue
p<-ggplot(data=genbank_dataDF,aes(x=log10(contig_n50),y=100*compleasm_complete/3390,color=type))+
  geom_point(shape=1, stroke=0.7, alpha=0.1) +  # Hollow colored rings
  geom_point(alpha=0.1) +  # Filled bubbles
  theme_bw()+
  annotate('rect',xmin=6,ymin=90,xmax=8.5,ymax=100,fill=NA,size=0.1)+
  labs(color="")+scale_color_manual(values=my_colors)+
  ylab('Compleasm Complete (%)')+xlab('log10 Contig N50')+
  ggtitle(paste0('Vertebrate Genomes N=',dim(genbank_dataDF)[1],', Vertebrata odb12 n=3390'))+
  ylim(c(00,100))+
  theme(
    legend.title = element_text(size = 14),  # Legend title font size
    legend.text = element_text(size = 12)    # Legend items font size
  )
ggsave('Fig2e_bubbleCompleteness.svg',plot=p,dpi=1200)