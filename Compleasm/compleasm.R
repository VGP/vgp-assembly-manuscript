library('ggplot2')
library(dplyr)
library(scales)
library(tidyr)
library(RColorBrewer)

vgp_stats<-read.table(file = 'vgp_compleasmStats.txt',sep = '\t',header = T)
vgp_cnStats<-read.table(file = 'vgp_asmStats.txt',sep='\t',header=T)
vgp_cnStatsFiltered<-vgp_cnStats[,c(1,4,5,6,7,8,9,10)]
vgp_merge<-merge(vgp_stats,vgp_cnStatsFiltered,by="Accession",all.x=TRUE)
# Create the dataframe with proper column assignment
vgp_df <- data.frame(
  seq = vgp_merge$Sequencing_type,
  fragmented = vgp_merge$Fragmented_compleasm,
  frameshift = vgp_merge$Frameshift_compleasm,
  missing = vgp_merge$Missing_compleasm,
  complete = vgp_merge$Complete_compleasm,
  date = vgp_merge$Submission_Date,
  class = vgp_merge$Class,
  superclass = vgp_merge$Super_class,
  superorder = vgp_merge$Super_order,
  infraclass = vgp_merge$Infra_class,
  order = vgp_merge$Order,
  dup = vgp_merge$Duplicated_compleasm,
  single = vgp_merge$Single_compleasm,
  proj = 'vgp',
  cn50 =vgp_merge$Contig_N50,
  asm_size=vgp_merge$Assembly_Size
)

dnazoo_stats<-read.table(file = 'dnazoo_compleasmStats.txt',sep = '\t',header = T)
dnazoo_cnStats<-read.table(file = 'dnazoo_asmStats.txt',sep='\t',header=T)
dnazoo_cnStatsFiltered<-dnazoo_cnStats[,c(1,4,5,6,7,8,9,10)]
dnazoo_merge<-merge(dnazoo_stats,dnazoo_cnStatsFiltered,by="Accession",all.x=TRUE)
# Create the dataframe with proper column assignment
dnazoo_df <- data.frame(
  seq = dnazoo_merge$Sequencing_type,
  fragmented = dnazoo_merge$Fragmented_compleasm,
  frameshift = dnazoo_merge$Frameshift_compleasm,
  missing = dnazoo_merge$Missing_compleasm,
  complete = dnazoo_merge$Complete_compleasm,
  date = dnazoo_merge$Submission_Date,
  class = dnazoo_merge$Class,
  superclass = dnazoo_merge$Super_class,
  superorder = dnazoo_merge$Super_order,
  infraclass = dnazoo_merge$Infra_class,
  order = dnazoo_merge$Order,
  dup = dnazoo_merge$Duplicated_compleasm,
  single = dnazoo_merge$Single_compleasm,
  proj = 'dnazoo',
  cn50 =dnazoo_merge$Contig_N50,
  asm_size=dnazoo_merge$Assembly_Size
)

zoonomia_stats<-read.table(file = 'zoonomia_compleasmStats.txt',sep = '\t',header = T)
zoonomia_cnStats<-read.table(file = 'zoonomia_asmStats.txt',sep='\t',header=T)
zoonomia_cnStatsFiltered<-zoonomia_cnStats[,c(1,4,5,6,7,8,9,10)]
zoonomia_merge<-merge(zoonomia_stats,zoonomia_cnStatsFiltered,by="Accession",all.x=TRUE)
# Create the dataframe with proper column assignment
zoonomia_df <- data.frame(
  seq = zoonomia_merge$Sequencing_type,
  fragmented = zoonomia_merge$Fragmented_compleasm,
  frameshift = zoonomia_merge$Frameshift_compleasm,
  missing = zoonomia_merge$Missing_compleasm,
  complete = zoonomia_merge$Complete_compleasm,
  date = zoonomia_merge$Submission_Date,
  class = zoonomia_merge$Class,
  superclass = zoonomia_merge$Super_class,
  superorder = zoonomia_merge$Super_order,
  infraclass = zoonomia_merge$Infra_class,
  order = zoonomia_merge$Order,
  dup = zoonomia_merge$Duplicated_compleasm,
  single = zoonomia_merge$Single_compleasm,
  proj = 'zoonomia',
  cn50 =zoonomia_merge$Contig_N50,
  asm_size=zoonomia_merge$Assembly_Size
)
# Define the target classes and filter the dataframe
target_classes <- c("Actinopteri", "Amphibia", "Aves", "Chondrichthyes",
                    "Lepidosauria", "Mammalia", "Testudinata","Crocodylia",
                    "Coelacanth","Lungfish","Cladistiaa")
filtered_df <- vgp_df[vgp_df$class %in% target_classes, ]

mammalia_df <- vgp_df[vgp_df$class=="Mammalia",]
aves_df <- vgp_df[vgp_df$class=="Aves",]

fish_df<-vgp_df[vgp_df$class %in% c("Actinopteri","Chondrichthyes","Coelacanth","Lungfish","Cladistia"),]
colors <- brewer.pal(8, "Set2")
# Plot the vgp dataframe
ggplot(data = filtered_df, aes(x = class, y = frameshift, color = class)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  ggtitle('Compleasm genes containing frameshifts') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(data = filtered_df, aes(x = class, y = frameshift, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  scale_fill_manual(values = colors) +
  ggtitle('Compleasm genes containing frameshifts') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(data = filtered_df, aes(x = seq, y = frameshift, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  ggtitle('Compleasm genes containing frameshifts') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(data = filtered_df, aes(x = seq, y = frameshift, color = class)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  ggtitle('Compleasm genes containing frameshifts') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot the complete BUSCOs
ggplot(data = filtered_df, aes(x = class, y = complete, color = class)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  ggtitle('Compleasm complete genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Plot the complete BUSCOs
ggplot(data = filtered_df, aes(x = class, y = complete, color = class)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +ylim(3100,3350)+
  ggtitle('Compleasm complete genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot the complete BUSCOs for mammals
ggplot(data = mammalia_df, aes(x = order, y = complete, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  ggtitle('Compleasm complete genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(data = mammalia_df, aes(x = order, y = frameshift, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  ggtitle('Compleasm genes with frameshifts') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot the complete BUSCOs for birds
ggplot(data = aves_df, aes(x = order, y = complete, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  ggtitle('Compleasm complete genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(data = aves_df, aes(x = order, y = frameshift, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  ggtitle('Compleasm genes with frameshifts') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot the complete BUSCOs for fish
ggplot(data = fish_df, aes(x = order, y = complete, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  ggtitle('Compleasm complete genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(data = fish_df, aes(x = order, y = frameshift, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  ggtitle('Compleasm genes with frameshifts') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##############################################################################

ggplot(data = dnazoo_df, aes(x = order, y = frameshift, color = order)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  ggtitle('Compleasm genes containing frameshifts') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = zoonomia_df, aes(x = order, y = frameshift, color = order)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  ggtitle('Compleasm genes containing frameshifts') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

comb_df<-rbind(mammalia_df,dnazoo_df,zoonomia_df)
ggplot(data = comb_df, aes(x = order, y = frameshift, color = proj)) +
  geom_jitter(width = 0.3, height = 0) +
  theme_bw() +ylim(50,300)+
  ggtitle('Compleasm genes containing frameshifts') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(data = comb_df, aes(x = order, y = complete, color = proj)) +
  geom_jitter(width = 0.3, height = 0) +
  theme_bw() +ylim(2500,3350)+
  ggtitle('Compleasm complete genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(data = comb_df, aes(x = order, y = missing, color = proj)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +ylim(0,300)+
  ggtitle('Compleasm missing genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(data = comb_df, aes(x = order, y = fragmented, color = proj)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +ylim(0,1000)+
  ggtitle('Compleasm fragmented genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data=comb_df,aes(x=log10(cn50),y=complete,color=proj,size=fragmented))+
  geom_point()+theme_bw()+xlab('Contig N50 (log10)')+ylab('Compleasm Complete')
