vgp_stats<-read.table(file = 'vgp_assemblyStats.tsv',sep = '\t',header = T)
library('ggplot2')
library(dplyr)
library(scales)
library(tidyr)
# Create the dataframe with proper column assignment
df <- data.frame(
  seq = vgp_stats$Sequencing_type,
  cn50 = vgp_stats$Contig_N50,
  cn90 = vgp_stats$Contig_N90,
  fragmented = vgp_stats$Fragmented_BUSCOs,
  nc = vgp_stats$Number_Contigs,
  gs = vgp_stats$Assembly_Size,
  stop = vgp_stats$Internal_stop_codons,
  missing = vgp_stats$Missing_BUSCOs,
  busco = vgp_stats$Complete_BUSCOs,
  date = vgp_stats$Submission_Date,
  class = vgp_stats$Class,
  superclass = vgp_stats$Super_class,
  superorder = vgp_stats$Super_order,
  infraclass = vgp_stats$Infra_class,
  order = vgp_stats$Order,
  dup = vgp_stats$Duplicate_BUSCOs,
  proj = 'vgp'
)

dnazoo_stats<-read.table(file = 'dnazoo_assemblyStats.txt',sep = '\t',header = T)
# Create the dataframe with proper column assignment
dnazoo_df <- data.frame(
  seq = dnazoo_stats$Sequencing_type,
  cn50 = dnazoo_stats$Contig_N50,
  cn90 = dnazoo_stats$Contig_N90,
  fragmented = dnazoo_stats$Fragmented_BUSCOs,
  nc = dnazoo_stats$Number_Contigs,
  gs = dnazoo_stats$Assembly_Size,
  stop = dnazoo_stats$Internal_stop_codons,
  missing = dnazoo_stats$Missing_BUSCOs,
  busco = dnazoo_stats$Complete_BUSCOs,
  date = dnazoo_stats$Submission_Date,
  class = dnazoo_stats$Class,
  superclass = dnazoo_stats$Super_class,
  superorder = dnazoo_stats$Super_order,
  infraclass = dnazoo_stats$Infra_class,
  order = dnazoo_stats$Order,
  dup = dnazoo_stats$Duplicate_BUSCOs,
  proj='dnazoo'
)

zoonomia_stats<-read.table(file = 'zoonomia_assemblyStats.txt',sep = '\t',header = T)
# Create the dataframe with proper column assignment
zoonomia_df <- data.frame(
  seq = zoonomia_stats$Sequencing_type,
  cn50 = zoonomia_stats$Contig_N50,
  cn90 = zoonomia_stats$Contig_N90,
  fragmented = zoonomia_stats$Fragmented_BUSCOs,
  nc = zoonomia_stats$Number_Contigs,
  gs = zoonomia_stats$Assembly_Size,
  stop = zoonomia_stats$Internal_stop_codons,
  missing = zoonomia_stats$Missing_BUSCOs,
  busco = zoonomia_stats$Complete_BUSCOs,
  date = zoonomia_stats$Submission_Date,
  class = zoonomia_stats$Class,
  superclass = zoonomia_stats$Super_class,
  superorder = zoonomia_stats$Super_order,
  infraclass = zoonomia_stats$Infra_class,
  order = zoonomia_stats$Order,
  dup = zoonomia_stats$Duplicate_BUSCOs,
  proj='zoonomia'
)
# Define the target classes and filter the dataframe
target_classes <- c("Actinopteri", "Amphibia", "Aves", "Chondrichthyes",
                    "Lepidosauria", "Mammalia", "Testudinata","Crocodylia",
                    "Coelacanth","Lungfish","Cladistia")
filtered_df <- df[df$class %in% target_classes, ]

mammalia_df <- df[df$class=="Mammalia",]

# Plot the filtered dataframe
ggplot(data = mammalia_df, aes(x = order, y = stop, color = order)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +ylim(c(0,360))+
  ggtitle('BUSCO genes containing STOP codons - Mammalia') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot the filtered dataframe
ggplot(data = mammalia_df, aes(x = order, y = nc, color = superorder)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  ggtitle('Number of contigs') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

aves_df <- df[df$class=="Aves",]

# Plot the filtered dataframe
ggplot(data = aves_df, aes(x = date, y = nc, color = seq)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  ggtitle('BUSCO genes containing STOP codons - Aves') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fish_df <- df[df$class=="Actinopteri",]
ggplot(data = fish_df, aes(x = order, y = stop, color = infraclass)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  ggtitle('BUSCO genes containing STOP codons - Actinopteri') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = filtered_df, aes(x = class, y = (busco), color = class)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() + ylim(c(3100,3350))+
  ggtitle('BUSCO completeness') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = filtered_df, aes(x = class, y = (missing), color = class)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +ylim(c(0,200))+
  ggtitle('Missing BUSCO genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = filtered_df, aes(x = class, y = dup, color = class)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() + 
  ggtitle('Duplicated BUSCO genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = filtered_df, aes(x = class, y = fragmented, color = class)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  ggtitle('Fragmented BUSCO genes') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = filtered_df, aes(x = class, y = stop, color = class)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  ggtitle('BUSCO genes containing STOP codons') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = filtered_df,aes(x=seq,y=log10(nc),color=seq))+
  geom_jitter(width=0.2,height=0)+theme_bw()+
  ggtitle('No. contigs per assembly')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = filtered_df,aes(x=class,y=log10(nc/gs),color=class))+
  geom_jitter(width=0.2,height=0)+theme_bw()+
  ggtitle('No. contigs per assembly per Gb')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = filtered_df,aes(x=seq,y=log10(cn50),color=seq))+
  geom_jitter(width=0.2,height=0)+theme_bw()+
  ggtitle('Contig N50')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = filtered_df,aes(x=class,y=log10(cn50/gs),color=seq))+
  geom_jitter(width=0.2,height=0)+theme_bw()+
  ggtitle('Contig N50 per Gb')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = filtered_df,aes(x=class,y=busco,color=class))+
  geom_jitter(width=0.2,height=0)+theme_bw()+
  ylim(c(2750,3350))+
  ggtitle('Complete BUSCO genes (vertebrata_odb10)')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = filtered_df,aes(x=seq,y=busco,color=class))+
  geom_jitter(width=0.2,height=0)+theme_bw()+
  ylim(c(2750,3350))+
  ggtitle('Complete BUSCO genes (vertebrata_odb10)')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = filtered_df,aes(x=class,y=dup,color=class))+
  geom_jitter(width=0.2,height=0)+theme_bw()+
  ylim(c(0,350))+
  ggtitle('Duplicated BUSCO genes (vertebrata_odb10)')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = filtered_df,aes(x=seq,y=dup,color=class))+
  geom_jitter(width=0.2,height=0)+theme_bw()+
  ylim(c(0,350))+
  ggtitle('Duplicated BUSCO genes (vertebrata_odb10)')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = filtered_df,aes(x=class,y=fragmented,color=class))+
  geom_jitter(width=0.2,height=0)+theme_bw()+
  ylim(c(0,250))+
  ggtitle('Fragmented BUSCO genes (vertebrata_odb10)')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = filtered_df,aes(x=seq,y=fragmented,color=class))+
  geom_jitter(width=0.2,height=0)+theme_bw()+
  ylim(c(0,250))+
  ggtitle('Fragmented BUSCO genes (vertebrata_odb10)')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = filtered_df,aes(x=class,y=log10(cn50),color=class))+
  geom_jitter(width = 0.2, height = 0)+theme_bw()+
  ggtitle('Log10-Contig N50')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = filtered_df,aes(x=seq,y=log10(cn50),color=class))+
  geom_jitter(width = 0.2, height = 0)+theme_bw()+
  ggtitle('Log10-Contig N50 by Sequencing platform and class')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = filtered_df,aes(x=class,y=log10(cn50/gs),color=class))+
  geom_jitter(width = 0.2, height = 0)+theme_bw()+
  ggtitle('Log-contigN50 divided by asm size')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = filtered_df,aes(x=seq,y=log10(cn50/gs),color=class))+
  geom_jitter(width = 0.2, height = 0)+theme_bw()+
  ggtitle('Log-contigN50 divided by asm size')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = filtered_df,aes(x=seq,y=log10(cn90),color=class))+
  geom_jitter(width = 0.2, height = 0)+theme_bw()+
  ggtitle('Log10-Contig N90 by Sequencing platform and class')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = filtered_df,aes(x=class,y=log10(cn90),color=class))+
  geom_jitter(width = 0.2, height = 0)+theme_bw()+
  ggtitle('Log10-Contig N90')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = filtered_df,aes(x=seq,y=log10(cn90/gs),color=class))+
  geom_jitter(width = 0.2, height = 0)+theme_bw()+
  ggtitle('Log10-Contig N90 divided by asm size')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = filtered_df,aes(x=class,y=log10(cn90/gs),color=class))+
  geom_jitter(width = 0.2, height = 0)+theme_bw()+
  ggtitle('Log10-Contig N90 divided by asm size')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

########################

df_comb<-rbind(filtered_df,dnazoo_df,zoonomia_df)

ggplot(data = df_comb[df_comb$class=="Mammalia",], aes(x = order, y = busco, color = proj)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw()+
  ggtitle('BUSCO Completeness') +
  ylab('Completeness')+
  xlab('Class')+scale_color_manual(values = c("#5DA5DA", "#FAA43A", "#60BD68")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = df_comb[df_comb$class=="Mammalia",], aes(x = order, y = stop, color = proj)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw() +
  ggtitle('BUSCO genes containing STOP codons') +scale_color_manual(values = c("#5DA5DA", "#FAA43A", "#60BD68")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot the filtered dataframe
ggplot(data = df_comb[df_comb$class=="Mammalia",], aes(x = order, y = log10(cn50), color = proj)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw()+
  ggtitle('Contig N50') +
  ylab('Log10(Contig N50)')+scale_color_manual(values = c("#5DA5DA", "#FAA43A", "#60BD68")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = df_comb[df_comb$class=="Mammalia",], aes(x = order, y = log10(cn50/gs), color = proj)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw()+
  ggtitle('Contig N50 / Genome Size') +
  ylab('Log10(Contig N50/Assembly Size)')+scale_color_manual(values = c("#5DA5DA", "#FAA43A", "#60BD68")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = df_comb[df_comb$class=="Mammalia",], aes(x = order, y = log10(nc), color = proj)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw()+
  ggtitle('Number of Contigs') +
  ylab('Log10(No. Contigs)')+scale_color_manual(values = c("#5DA5DA", "#FAA43A", "#60BD68")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(data = df_comb[df_comb$class=="Mammalia",], aes(x = order, y = log10(nc/gs), color = proj)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw()+
  ggtitle('Number of Contigs per bp') +
  ylab('Log10(No. Contigs/Assembly Size)')+scale_color_manual(values = c("#5DA5DA", "#FAA43A", "#60BD68")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


########################

df_cum<-data.frame(df$class,df$seq,as.Date(df$date),df$gs)
colnames(df_cum)<-c('class','seq','date','gs')
df_cum<-df_cum[order(df$date),]
df_cum$cum_size<-cumsum(df_cum$gs)
ggplot(data=df_cum,aes(x=date,y=cum_size,color=seq))+geom_point()+
  scale_y_continuous(labels = label_number(scale = 1e-9, suffix = " GB")) +
  scale_x_date(breaks = seq(from = as.Date("2019-01-01"), to = as.Date("2024-01-01"), by = "1 year"),
               date_labels = "%Y")+
  ylab('Cumulative Genome Span')+theme_bw()
  xlab('Date')
  
df_cum$g_num<-c(1:472)
ggplot(data=df_cum,aes(x=date,y=g_num,color=seq))+geom_point()+
  scale_x_date(breaks = seq(from = as.Date("2018-11-01"), 
                            to = as.Date("2024-11-01"), by = "1 year"),date_labels = "%Y")+
  ylab('Number of genomes')+theme_bw()+xlab('Date')
ggplot(data=df_cum,aes(x=date,y=g_num,color=class))+geom_point()+
  scale_x_date(breaks = seq(from = as.Date("2019-01-01"), to = as.Date("2024-01-01"), by = "1 year"),
    date_labels = "%Y")+
  ylab('Number of genomes')+theme_bw()+xlab('Date')

# Ensure date is in Date format, then extract the year
df <- df %>%
  mutate(date = as.Date(date),
         year = format(date, "%Y"))  # Extracting the year as a new column
filtered_df <- filtered_df %>%
  mutate(date = as.Date(date),
         year = format(date, "%Y"))  # Extracting the year as a new column

# Count entries per year and sequencing type
df_counts <- df %>%
  group_by(year, seq) %>%
  summarise(count = n()) %>%
  ungroup()

df_counts <- df_counts %>%
  complete(year, seq, fill = list(count = 0))

ggplot(df_counts, aes(x = year, y = count, fill = seq)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Number of Genomes per Year by Sequencing Type",
       x = "Year",
       y = "Number of Genomes",
       fill = "Sequencing Type") +
  scale_y_continuous(labels = comma) +  # Use comma format for large numbers
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Count entries per year and sequencing type
df_counts <- filtered_df %>%
  group_by(year, class) %>%
  summarise(count = n()) %>%
  ungroup()
df_counts<-df_counts[df_counts$class %in% c('Actinopteri','Amphibia','Aves','Lepidosauria','Mammalia'),]
df_counts <- df_counts %>%
  complete(year, class, fill = list(count = 0))

ggplot(df_counts, aes(x = year, y = count, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Number of Genomes per Year by Class",
       x = "Year",
       y = "Number of Genomes",
       fill = "Class") +
  scale_y_continuous(labels = comma) +  # Use comma format for large numbers
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

df_class_counts <- df %>%
  group_by(class) %>%
  summarise(count = n()) %>%
  ungroup()

# Plot the bar chart
ggplot(df_class_counts, aes(x = reorder(class, -count), y = count, fill = class)) +
  geom_bar(stat = "identity") +
  labs(title = "Total Number of Genomes per Class",
       x = "Class",
       y = "Number of Genomes") +
  scale_y_continuous(labels = comma) +  # Use comma format for large numbers
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")  # Remove legend as it's redundant
