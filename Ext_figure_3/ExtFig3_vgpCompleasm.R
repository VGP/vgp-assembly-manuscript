
library('ggplot2')
library(dplyr)
library(scales)
library(tidyr)
library(lubridate)

vgp_compleasmStats<-read.table(file = 'vgp_compleasmStats.tsv',header = T,sep = '\t')
head(vgp_compleasmStats)

vgp_compleasmDF<-data.frame(vgp_compleasmStats)

#Plot number of genes containing frameshift mutations separated by long-read sequencing technology
p<-ggplot(data = vgp_compleasmDF,
       aes(x = technology, y = 100*compleasm_frameshift/3390, fill = technology)) +
  geom_violin(aes(color = technology), alpha = 0.2, scale = 'width') +scale_color_manual(
    values = c(
      "CLR"  = "#E69F00",  # orange
      "HiFi" = "#56B4E9",  # sky blue
      "ONT"  = "#009E73"   # green
    )
  ) +scale_fill_manual(
    values = c(
      "CLR"  = "#E69F00",  # orange
      "HiFi" = "#56B4E9",  # sky blue
      "ONT"  = "#009E73"   # green
    )
  ) +
  #scale_fill_manual(values = c('darkgreen','lightblue','pink')) +
  #scale_color_manual(values = c('darkgreen','lightblue','pink')) +
  geom_jitter(aes(color = technology), height = 0, alpha = 0.3) +
  theme_bw() + ylab('Genes with Frameshift (%)')+xlab('Sequencing Technology')+
  scale_x_discrete(guide = guide_axis(angle = 45))
ggsave('Ext3b_frameshiftsPerTech.svg',plot=p,dpi=1200)

my_colors <- c(
  "Mammalia"       = "#E69F00",  # blue
  "Aves"           = "#00796B",  # orange
  "Actinopteri" = "#56B4E9",  # green
  "Cladistia" = "#56B4E9",  # green
  "Lepidosauria"       = "#80CBC4",  # red
  "Amphibia"       = "#984EA3",  # purple
  "Turtles"      = "#4DB6AC",  # brown
  "Crocodiles"      = "#009688",  # brown
  "Chondrichthyes"  = "#0072B2",  # cyan
  "Coelacanth"     = "#A6761D",
  "Hyperoartia"     = "#CC79A7",
  "Lungfish"     = "#A6761D",
  "Myxini"     = "#CC79A7"
)

custom_order <- c(
  "Myxini",
  "Hyperoartia",
  "Chondrichthyes",
  "Actinopteri",
  "Cladistia",
  "Coelacanth",
  "Lungfish",
  "Amphibia",
  "Mammalia",
  "Lepidosauria",
  "Turtles",
  "Crocodiles",
  "Aves"
)
vgp_compleasmDF$plot_class<-factor(vgp_compleasmDF$plot_class,levels=c("Mammals","Birds","Crocodiles","Turtles","Lepidosauria",
                                                           "Amphibians","Lobe-finned Fishes","Ray-finned Fishes",
                                                           "Cartilaginous Fishes","Cyclostomes"))

#Plot Compleasm completness of each genome separated by class - with and without the crocodiles (only 2 genomes), cyclostomes and lobe-finned fishes (both classes exhibit low completeness likely due to large genome sizes or divergence from vertebrate busco lineage)
p<-ggplot(data = vgp_compleasmDF[vgp_compleasmDF$plot_class %in% c('Cartilaginous Fishes','Ray-finned Fishes','Amphibians','Lepidosauria','Turtles','Birds','Mammals'),], aes(x = plot_class, fill=plot_class,y = 100*(compleasm_complete)/3390, color = plot_class)) +
  geom_violin(aes(color = plot_class), alpha = 0.2, scale = 'width') +
  scale_color_manual(values = plot_colors) +
  scale_fill_manual(values = plot_colors) +
  geom_jitter(aes(color = plot_class), height = 0, alpha = 0.2) +
  theme_bw() +ylab('Compleasm Complete (%)')+xlab('')+
  scale_x_discrete(guide = guide_axis(angle = 45))
ggsave('Ext3c_completenessClassReduced',plot = p,dpi=1200)

p<-ggplot(data = vgp_compleasmDF, aes(x = plot_class, fill=plot_class,y = 100*(compleasm_complete)/3390, color = plot_class)) +
  geom_violin(aes(color = plot_class), alpha = 0.2, scale = 'width') +
  scale_color_manual(values = plot_colors) +
  scale_fill_manual(values = plot_colors) +
  geom_jitter(aes(color = plot_class), height = 0, alpha = 0.2) +
  theme_bw() +ylab('Compleasm Complete (%)')+xlab('')+
  scale_x_discrete(guide = guide_axis(angle = 45))
ggsave('Ext3d_completenessClass.svg',plot = p,dpi=1200)

#Plot Contig N50 vs Compleasm completeness for all VGP vertebrate genomes, coloured by extended lineage class and the size of each bubble defined by the percentage of BUSCO genes found fragmented in the genome
p<-ggplot(data=vgp_compleasmDF,aes(x=log10(contig_n50),y=100*compleasm_complete/3390,color=class,size=100*(compleasm_fragmented/3390)))+
  geom_point(alpha=0.5)+
  theme_bw()+
  annotate('rect',xmin=6,ymin=90,xmax=8.5,ymax=100,fill=NA,color='black',size=0.1)+
  ggtitle(paste0('VGP Genomes N=',dim(vgp_compleasmDF)[1],', Vertebrata odb12 n=3,390'))+
  labs(size="Compleasm Fragmented (%)",color="Category")+
  ylab('Compleasm Complete (%)')+xlab('log10 Contig N50')+
  scale_color_manual(values = my_colors) +#ylim(c(65,100))+
  theme(
    legend.title = element_text(size = 14),  # Legend title font size
    legend.text = element_text(size = 12)    # Legend items font size
  )
ggsave('Ext3e_vgpBubblePerClass.svg',plot=p,dpi=1200)