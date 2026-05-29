library(ggplot2)
library(matrixStats)
library(dplyr)

my_colors <- c(
  "Mammals"              = "#E69F00",
  "Birds"                = "#00796B",
  "Ray-finned Fishes"    = "#56B4E9",
  "Lepidosauria"         = "#80CBC4",
  "Amphibians"           = "#984EA3",
  "Turtles"              = "#4DB6AC",
  "Crocodiles"           = "#009688",
  "Cartilaginous Fishes" = "#0072B2",
  "Lobe-finned Fishes"   = "#A6761D",
  "Cyclostomes"          = "#CC79A7"
)

compute_summary <- function(files, group, type){
  stats <- c()
  for (f in files){
    stats_in <- read.delim(f, sep='\t', header=TRUE)
    stats <- cbind(stats, stats_in[,2])
  }
  if (ncol(stats) == 0) return(NULL) # skip empty groups
  quantiles <- rowQuantiles(stats, probs=c(0.05, 0.95))
  data.frame(
    Nx     = (1:1000)/1000,
    mean   = rowMeans(stats),
    median = rowMedians(stats),
    q5     = quantiles[,1],
    q95    = quantiles[,2],
    Group  = group,
    Type   = type
  )
}

# ---- Compute scaffold summaries ----
vgp_scaffolds <- bind_rows(
  compute_summary(list.files('./ucsc/reptiles/',       pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Lepidosauria", "Scaffold"),
  compute_summary(list.files('./ucsc/cyclostomes/',    pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Cyclostomes", "Scaffold"),
  compute_summary(list.files('./ucsc/lobe-finned_fishes/', pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Lobe-finned Fishes", "Scaffold"),
  compute_summary(list.files('./ucsc/crocodiles/',     pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Crocodiles", "Scaffold"),
  compute_summary(list.files('./ucsc/turtles/',        pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Turtles", "Scaffold"),
  compute_summary(list.files('./ucsc/birds/',          pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Birds", "Scaffold"),
  compute_summary(list.files('./ucsc/mammals/',        pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Mammals", "Scaffold"),
  compute_summary(list.files('./ucsc/sharks/',         pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Cartilaginous Fishes", "Scaffold"),
  compute_summary(list.files('./ucsc/fish/',           pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Ray-finned Fishes", "Scaffold"),
  compute_summary(list.files('./ucsc/amphibians/',           pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Amphibians", "Scaffold")
)

# ---- Compute contig summaries ----
vgp_contigs <- bind_rows(
  compute_summary(list.files('./ucsc/reptiles/',       pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Lepidosauria", "Contig"),
  compute_summary(list.files('./ucsc/cyclostomes/',    pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Cyclostomes", "Contig"),
  compute_summary(list.files('./ucsc/lobe-finned_fishes/', pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Lobe-finned Fishes", "Contig"),
  compute_summary(list.files('./ucsc/crocodiles/',     pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Crocodiles", "Contig"),
  compute_summary(list.files('./ucsc/turtles/',        pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Turtles", "Contig"),
  compute_summary(list.files('./ucsc/birds/',          pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Birds", "Contig"),
  compute_summary(list.files('./ucsc/mammals/',        pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Mammals", "Contig"),
  compute_summary(list.files('./ucsc/sharks/',         pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Cartilaginous Fishes", "Contig"),
  compute_summary(list.files('./ucsc/fish/',           pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Ray-finned Fishes", "Contig"),
  compute_summary(list.files('./ucsc/amphibians/',           pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Amphibians", "Contig")
)

# ---- Compute scaffold summaries ----
other_scaffolds <- bind_rows(
  compute_summary(list.files('./all/reptiles/',       pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Lepidosauria", "Scaffold"),
  compute_summary(list.files('./all/cyclostomes/',    pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Cyclostomes", "Scaffold"),
  compute_summary(list.files('./all/lobe-finned_fishes/', pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Lobe-finned Fishes", "Scaffold"),
  compute_summary(list.files('./all/crocodiles/',     pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Crocodiles", "Scaffold"),
  compute_summary(list.files('./all/turtles/',        pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Turtles", "Scaffold"),
  compute_summary(list.files('./all/birds/',          pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Birds", "Scaffold"),
  compute_summary(list.files('./all/mammals/',        pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Mammals", "Scaffold"),
  compute_summary(list.files('./all/sharks/',         pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Cartilaginous Fishes", "Scaffold"),
  compute_summary(list.files('./all/fish/',           pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Ray-finned Fishes", "Scaffold"),
  compute_summary(list.files('./all/amphibians/',           pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Amphibians", "Scaffold")
)

# ---- Compute contig summaries ----
other_contigs <- bind_rows(
  compute_summary(list.files('./all/reptiles/',       pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Lepidosauria", "Contig"),
  compute_summary(list.files('./all/cyclostomes/',    pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Cyclostomes", "Contig"),
  compute_summary(list.files('./all/lobe-finned_fishes/', pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Lobe-finned Fishes", "Contig"),
  compute_summary(list.files('./all/crocodiles/',     pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Crocodiles", "Contig"),
  compute_summary(list.files('./all/turtles/',        pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Turtles", "Contig"),
  compute_summary(list.files('./all/birds/',          pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Birds", "Contig"),
  compute_summary(list.files('./all/mammals/',        pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Mammals", "Contig"),
  compute_summary(list.files('./all/sharks/',         pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Cartilaginous Fishes", "Contig"),
  compute_summary(list.files('./all/fish/',           pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Ray-finned Fishes", "Contig"),
  compute_summary(list.files('./all/amphibians/',           pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Amphibians", "Contig")
)

# ---- Combine everything ----
plot_df <- bind_rows(vgp_scaffolds, vgp_contigs)
plot_df_other <- bind_rows(other_scaffolds, other_contigs)
plot_df$Project="VGP"
plot_df_other$Project="Other"
comb_df<-rbind(plot_df,plot_df_other)
# ---- Plot ----

ggplot(comb_df[comb_df$Type == "Scaffold", ], 
       aes(x = Nx, y = log10(median), color = Group, linetype = Project, alpha = Project)) +
  geom_line() +
  scale_color_manual(values = my_colors) +
  scale_linetype_manual(values = c("VGP" = "solid", "Other" = "dotted")) +
  scale_alpha_manual(values = c("VGP" = 1, "Other" = 1)) +  # reduce transparency for "Other"
  labs(
    x = "Genome fraction",
    y = "log10(Scaffold length)",
    title = "Median Scaffold Nx curves"
  ) +
  scale_y_continuous(
    breaks = log10(c(1e2,1e3,1e4, 1e5, 1e6, 1e7, 1e8,1e9)),   # tick positions in log space
    labels = c("100 bp","1 kb","10 kb", "100 kb", "1 Mb", "10 Mb", "100 Mb","1 Gb")
  ) +
  theme_minimal()
#ggsave(filename = 'scaffold_medianNx.svg',plot = p,dpi = 1200)
ggplot(comb_df[comb_df$Type == "Contig", ], 
       aes(x = Nx, y = log10(median), color = Group, linetype = Project, alpha = Project)) +
  geom_line() +
  scale_color_manual(values = my_colors) +
  scale_linetype_manual(values = c("VGP" = "solid", "Other" = "dotted")) +
  scale_alpha_manual(values = c("VGP" = 1, "Other" = 1)) +  # reduce transparency for "Other"
  labs(
    x = "Genome fraction",
    y = "log10(Contig length)",
    title = "Median Contig Nx curves"
  ) +
  scale_y_continuous(
    breaks = log10(c(1e2,1e3,1e4, 1e5, 1e6, 1e7, 1e8,1e9)),   # tick positions in log space
    labels = c("100 bp","1 kb","10 kb", "100 kb", "1 Mb", "10 Mb", "100 Mb","1 Gb")
  ) +
  theme_minimal()
#ggsave(filename = 'contig_medianNx.svg',plot = p,dpi = 1200)
ggplot(comb_df[comb_df$Type == "Contig", ], 
       aes(x = Nx, y = log10(mean), color = Group, linetype = Project, alpha = Project)) +
  geom_line() +
  scale_color_manual(values = my_colors) +
  scale_linetype_manual(values = c("VGP" = "solid", "Other" = "dotted")) +
  scale_alpha_manual(values = c("VGP" = 1, "Other" = 1)) +  # reduce transparency for "Other"
  labs(
    x = "Genome fraction",
    y = "log10(Contig length)",
    title = "Mean Contig Nx curves"
  ) +
  scale_y_continuous(
    breaks = log10(c(1e2,1e3,1e4, 1e5, 1e6, 1e7, 1e8,1e9)),   # tick positions in log space
    labels = c("100 bp","1 kb","10 kb", "100 kb", "1 Mb", "10 Mb", "100 Mb","1 Gb")
  ) +
  theme_minimal()
#ggsave(filename = 'contig_meanNx.svg',plot = p,dpi = 1200)
ggplot(comb_df[comb_df$Type == "Scaffold", ], 
       aes(x = Nx, y = log10(mean), color = Group, linetype = Project, alpha = Project)) +
  geom_line() +
  scale_color_manual(values = my_colors) +
  scale_linetype_manual(values = c("VGP" = "solid", "Other" = "dotted")) +
  scale_alpha_manual(values = c("VGP" = 1, "Other" = 1)) +  # reduce transparency for "Other"
  labs(
    x = "Genome fraction",
    y = "log10(Scaffold length)",
    title = "Mean Scaffold Nx curves"
  ) +
  scale_y_continuous(
    breaks = log10(c(1e2,1e3,1e4, 1e5, 1e6, 1e7, 1e8,1e9)),   # tick positions in log space
    labels = c("100 bp","1 kb","10 kb", "100 kb", "1 Mb", "10 Mb", "100 Mb","1 Gb")
  ) +
  theme_minimal()
#ggsave(filename = 'scaffold_meanNx.svg',plot = p,dpi = 1200)

ggplot(comb_df[comb_df$Type=="Scaffold" & comb_df$Group %in% c("Birds","Mammals","Ray-finned Fishes","Amphibians"),], aes(x=Nx, color=Group, fill=Group,linetype=Project)) +
  # Shaded area between q5 and q95
  geom_ribbon(aes(ymin=log10(q5), ymax=log10(q95)), alpha=0.2, color=NA) +
  # Bold median line
  geom_line(aes(y=log10(median)), size=1) +
  scale_color_manual(values=my_colors) +
  scale_fill_manual(values=my_colors) +
  labs(
    x = "Genome fraction",
    y = "log10(Scaffold length)",
    title = "Median Nx curves with 5–95% quantile shading"
  ) +
  ylim(c(1,9.5)) +
  theme_minimal()

ggplot(comb_df[comb_df$Type=="Contig" & comb_df$Group %in% c("Birds","Mammals","Ray-finned Fishes","Amphibians"),], aes(x=Nx, color=Group, fill=Group,linetype=Project)) +
  # Shaded area between q5 and q95
  geom_ribbon(aes(ymin=log10(q5), ymax=log10(q95)), alpha=0.2, color=NA) +
  # Bold median line
  geom_line(aes(y=log10(median)), size=1) +
  scale_color_manual(values=my_colors) +
  scale_fill_manual(values=my_colors) +
  labs(
    x = "Genome fraction",
    y = "log10(Contig length)",
    title = "Median Nx curves with 5–95% quantile shading"
  ) +
  #ylim(c(1,9.5)) +
  theme_minimal()


ggplot(comb_df[comb_df$Group %in% c("Amphibians") & comb_df$Type=="Contig",], aes(x=Nx, color=Group, fill=Group,linetype=Project)) +
  # Shaded area between q5 and q95
  geom_ribbon(aes(ymin=log10(q5), ymax=log10(q95)), alpha=0.2, color=NA) +
  # Bold median line
  geom_line(aes(y=log10(median)), size=1) +
  scale_color_manual(values=my_colors) +
  scale_fill_manual(values=my_colors) +
  scale_linetype_manual(values = c("VGP" = "solid", "Other" = "dotted")) +
  labs(
    x = "Genome fraction",
    y = "log10(Contig length)",
    title = "Median Nx curves with 5–95% quantile shading"
  ) +
  ylim(c(1,8.5)) +
  theme_minimal()

ggplot(comb_df[comb_df$Group %in% c("Amphibians") & comb_df$Type=="Scaffold",], aes(x=Nx, color=Group, fill=Group,linetype=Project)) +
  # Shaded area between q5 and q95
  geom_ribbon(aes(ymin=log10(q5), ymax=log10(q95)), alpha=0.2, color=NA) +
  # Bold median line
  geom_line(aes(y=log10(median)), size=1) +
  scale_color_manual(values=my_colors) +
  scale_fill_manual(values=my_colors) +
  scale_linetype_manual(values = c("VGP" = "solid", "Other" = "dotted")) +
  labs(
    x = "Genome fraction",
    y = "log10(Scaffold length)",
    title = "Median Nx curves with 5–95% quantile shading"
  ) +
  theme_minimal()

  
