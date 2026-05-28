
library('ggplot2')
library(dplyr)
library(scales)
library(tidyr)
library(lubridate)
library('svglite')
asmStats_df<-read.table(file = 'vgp_asmMetaData.tsv',header = T,sep = '\t')

#Show length of each assembly split by class - define the plot groups for e.g. turtles and lunfish
asmStats_df$class[asmStats_df$species=="Protopterus annectens"]<-'Lungfish'
asmStats_df$class[asmStats_df$species=="Latimeria chalumnae"]<-'Coelacanth'
asmStats_df$class[asmStats_df$group=='turtles']<-'Turtles'
asmStats_df$class[asmStats_df$group=='vertebrates']<-'Crocodiles'
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
asmStats_df$class <- factor(asmStats_df$class, levels = custom_order)
asmStats_df$plot_class<-as.character(asmStats_df$class)
asmStats_df$plot_class[asmStats_df$class=='Lepidosauria']='Lepidosauria'
asmStats_df$plot_class[asmStats_df$class=='Turtles']='Turtles'
asmStats_df$plot_class[asmStats_df$class=='Crocodiles']='Crocodiles'
asmStats_df$plot_class[asmStats_df$class=='Actinopteri']='Ray-finned Fishes'
asmStats_df$plot_class[asmStats_df$class=='Cladistia']='Ray-finned Fishes'
asmStats_df$plot_class[asmStats_df$class=='Coelacanth']='Lobe-finned Fishes'
asmStats_df$plot_class[asmStats_df$class=='Hyperoartia']='Cyclostomes'
asmStats_df$plot_class[asmStats_df$class=='Myxini']='Cyclostomes'
asmStats_df$plot_class[asmStats_df$class=='Lungfish']='Lobe-finned Fishes'
asmStats_df$plot_class[asmStats_df$class=='Amphibia']='Amphibians'
asmStats_df$plot_class[asmStats_df$class=='Aves']='Birds'
asmStats_df$plot_class[asmStats_df$class=='Chondrichthyes']='Cartilaginous Fishes'
asmStats_df$plot_class[asmStats_df$class=='Mammalia']='Mammals'
asmStats_df$plot_class<-factor(asmStats_df$plot_class,levels=c("Mammals","Birds","Crocodiles","Turtles","Lepidosauria",
                                                               "Amphibians","Lobe-finned Fishes","Ray-finned Fishes",
                                                               "Cartilaginous Fishes","Cyclostomes"))
plot_colors <- c(
  "Mammals"       = "#E69F00", 
  "Birds"           = "#00796B",
  "Crocodiles"      = "#009688",
  "Turtles"      = "#4DB6AC", 
  "Lepidosauria"       = "#80CBC4",
  "Amphibians"       = "#984EA3",
  "Lobe-finned Fishes"     = "#A6761D",
  "Ray-finned Fishes" = "#56B4E9",
  "Cartilaginous Fishes"  = "#0072B2",
  "Cyclostomes"     = "#CC79A7"
)
#Add assembly lengths from non-vertebrates

genome_sizeDF<-data.frame(asmStats_df$accession,asmStats_df$plot_class,as.numeric(asmStats_df$assembly_size))
colnames(genome_sizeDF)<-c('Accession','plot_class','size')
extra_genomes<-read.delim('OtherDeuterostomes_genomeSizes.txt',header = F)
extra_genomesdf<-data.frame(Accession=extra_genomes[,1],plot_class=rep('Other Deuterostomes',13),size=extra_genomes[,2])
genome_sizesDFNew<-rbind(genome_sizeDF,extra_genomesdf)

p<-ggplot(
  data = genome_sizesDFNew,
  aes(y = size, x = plot_class, fill = plot_class, color = plot_class)
) +
  geom_violin(position = 'identity', alpha = 0.5, scale = 'width') +
  geom_jitter(height = 0) +
  scale_fill_manual(values = plot_colors) +   # <- controls violin fill
  scale_color_manual(values = plot_colors) +  # <- controls outline + points
  theme_bw() +
  theme(axis.text.y = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
  ylab('Genome Size (Log10)') + xlab('') +
  labs(fill = NULL) +
  scale_y_log10(breaks = c(1e9, 1e10), labels = c("1 Gb", "10 Gb")) +
  theme(legend.position = "none") +
  coord_flip() 
ggsave('Fig2a_assemblySizes.svg',plot=p,dpi=1200)