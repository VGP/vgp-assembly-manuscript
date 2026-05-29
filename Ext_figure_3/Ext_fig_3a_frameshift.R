
library('ggplot2')
library(dplyr)
library(scales)
library(tidyr)
library(stringr)

genbank_dataDF<-read.table(file = 'genbank_genomesData.tsv',header = T,sep = '\t',colClasses = c(accession='character'))
my_colors <- c(
  "#B0B0B0","#0072B2","" # highlight
)
#Show percentage of BUSCO genes detected in each genome containing frameshift mutations for VGP genomes and other vertebrate genomes on GenBank
p<-ggplot(data=genbank_dataDF,aes(x=type,y=100*compleasm_frameshift/3390,color=type,fill=type))+
  geom_jitter(height=0,alpha=0.1)+ylab("Genes with frameshift (%)")+xlab('')+scale_color_manual(values=my_colors)+scale_fill_manual(values=my_colors)+
  theme_bw()+geom_violin(alpha=0.4,scale = 'width')+ylim(c(0,20))
ggsave('Ext3a_frameshiftGenbank',plot=p,dpi=1200)

#Show genome completeness vs contig N50, coloured by "large" sequencing project
p<-ggplot(data=genbank_dataDF,aes(x=log10(contig_n50),y=100*compleasm_complete/3390,color=project))+
  geom_point(shape=1, stroke=0.7, alpha=0.4) +  # Hollow colored rings
  geom_point(alpha=0.3) +  # Filled bubbles
  theme_bw()+#annotate('rect',xmin=6,ymin=90,xmax=8.5,ymax=100,fill=NA,color='black',size=0.1)+
  ylab('Compleasm Complete (%)')+xlab('log10 Contig N50')+
  ylim(c(00,100))+scale_color_manual(values=c('yellow','red','purple','grey','#0072B2','darkgreen'))+
  theme(
    legend.title = element_text(size = 14),  # Legend title font size
    legend.text = element_text(size = 12)    # Legend items font size
  )
ggsave('Ext3g_bubblePerProject',plot=p,dpi=1200)

#Show the same plot, but this time only colouring the VGP set and those from the same species used for comparisons in other analysis

ext_comp<-data.frame(read.delim(file = 'ext_compar.txt',header = T))
ext_acc<-c()
for (i in c(1:dim(ext_comp)[1])){
  ext_acc<-c(ext_acc,unlist(strsplit((ext_comp[i,]),"_"))[2])
}
for (i in c(1:dim(ext_comp)[1])){
  ext_acc[i]<-unlist(strsplit((ext_acc[i]),"\\."))[1]
}
my_colors <- c(
  "#E69F00","#B0B0B0","#0072B2","" # highlight
)
genbank_dataDF$type[genbank_dataDF$accession %in% ext_acc]="Comparison"
p<-ggplot(genbank_dataDF, aes(
  x = log10(contig_n50),
  y = 100 * compleasm_complete / 3390,
  color = type,
  alpha = type  # map alpha to Type
)) +
  geom_point(shape = 16, stroke = 0.7) +
  theme_bw() +
  labs(color = "", alpha = "") +
  ylab("Compleasm Complete (%)") +
  xlab("log10 Contig N50") +
  ylim(c(0, 100)) +
  scale_color_manual(values = my_colors) +
  scale_alpha_manual(values = c(
    "Other" = 0.1,       # background grey → very faint
    "Comparison" = 0.5,  # orange → stronger
    "VGP" = 0.5          # blue → strongest
  )) +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
ggsave('Ext3f_bubbleWithComparators',plot=p,dpi=1200)