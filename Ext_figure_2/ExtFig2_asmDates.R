
library('ggplot2')
library(dplyr)
library(scales)
library(tidyr)
library(lubridate)
library('svglite')
asmStats_df<-read.table(file = 'vgp_asmMetaData.tsv',header = T,sep = '\t')

#Plot assemblies produced per year coloured by long-read sequencing technology
p<-asmStats_df %>%
  mutate(date = as.Date(date)) %>%
  arrange(date) %>%
  mutate(cum_n = row_number()) %>%
  ggplot(aes(x = date, y = cum_n)) +
  geom_point(aes(color = technology), alpha = 0.7, size = 2) + # points by tech
  scale_color_manual(
    values = c(
      "CLR"  = "#E69F00",  # orange
      "HiFi" = "#56B4E9",  # sky blue
      "ONT"  = "#009E73"   # green
    )
  ) +
  theme_bw() +
  labs(
    x = "Date",
    y = "Cumulative number of genomes",
    color = "Technology"
  )
ggsave(filename = 'Ext2a_cumulative_genomes.svg',plot = p,dpi = 1200)

#Plot assemblies produced per year coloured by long-read sequencing technology, with the y-axis showing cumulative assembly length
p<-asmStats_df %>%

  mutate(date = as.Date(date)) %>%
  arrange(date) %>%
  mutate(assembly_size = assembly_size/1e9, cum_size = cumsum(assembly_size)) %>%
  ggplot(aes(x = date, y = cum_size)) +
  geom_point(aes(color = technology), alpha = 0.8, size = 2) + # points by tech
  scale_color_manual(
    values = c(
      "CLR"  = "#E69F00",  # orange
      "HiFi" = "#56B4E9",  # sky blue
      "ONT"  = "#009E73"   # green
    )
  ) +
  theme_bw() +
  labs(
    x = "Date",
    y = "Cumulative genome length (Gb)",
    color = "Technology"
  )
ggsave(filename = 'Ext2b_cumulative_genomeSpan.svg',plot = p,dpi = 1200)
#Plot assemblies produced per year coloured by long-read sequencing technology - this time split by sequencing technology
p<-asmStats_df %>%
  mutate(date = as.Date(date)) %>%  # ensure it's a Date
  group_by(technology, date) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(technology, date) %>%
  group_by(technology) %>%
  mutate(cum_n = cumsum(n)) %>%
  ggplot(aes(x = date, y = cum_n, color = technology)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(
    values = c(
      "CLR"  = "#E69F00",  # orange
      "HiFi" = "#56B4E9",  # sky blue
      "ONT"  = "#009E73"   # green
    )
  ) +
  labs(
    x = "Date",
    y = "Cumulative number of genomes",
    color = "Technology"
  )
ggsave(filename = 'Ext2c_cumulative_genomesSplitTech.svg',plot = p,dpi = 1200)
#Plot assemblies produced per year coloured by long-read sequencing technology, with the y-axis showing cumulative assembly length, again split by sequencing technology
p<-asmStats_df %>%

  mutate(date = as.Date(date)) %>%  # ensure it's a Date
  group_by(technology, date) %>%
  summarise(assembly_size = sum(assembly_size) / 1e9, .groups = "drop") %>%  # total size per day, in Gb
  arrange(technology, date) %>%
  group_by(technology) %>%
  mutate(cum_size_gb = cumsum(assembly_size)) %>%
  ggplot(aes(x = date, y = cum_size_gb, color = technology)) +
  geom_point(size = 1) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "CLR"  = "#E69F00",  # orange
      "HiFi" = "#56B4E9",  # sky blue
      "ONT"  = "#009E73"   # green
    )
  ) +
  labs(
    x = "Date",
    y = "Cumulative genome size (Gb)",
    color = "Technology"
  )
ggsave(filename = 'Ext2d_cumulative_genomeSpanSplitTech.svg',plot = p,dpi = 1200)

#Show number of genomes publisehd per year split by class - define the plot groups for e.g. turtles and lunfish
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
p<-asmStats_df %>%
  mutate(
    date = as.Date(date),
    year = year(date)      # extract year
  ) %>%
  group_by(year, plot_class) %>%
  summarise(n = n(), .groups = "drop") %>%
  ggplot(aes(x = year, y = n, fill = plot_class)) +
  geom_col() +
  scale_fill_manual(values = plot_colors) +
  theme_bw() +
  labs(
    x = "Year",
    y = "Number of genomes",
    fill = "Class"
  )
ggsave(filename = 'Ext2e_asms_perYearByClass.svg',plot = p,dpi = 1200)
