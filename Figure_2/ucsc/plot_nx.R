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
scaffolds <- bind_rows(
  compute_summary(list.files('./reptiles/',       pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Lepidosauria", "Scaffold"),
  compute_summary(list.files('./cyclostomes/',    pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Cyclostomes", "Scaffold"),
  compute_summary(list.files('./lobe-finned_fishes/', pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Lobe-finned Fishes", "Scaffold"),
  compute_summary(list.files('./crocodiles/',     pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Crocodiles", "Scaffold"),
  compute_summary(list.files('./turtles/',        pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Turtles", "Scaffold"),
  compute_summary(list.files('./birds/',          pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Birds", "Scaffold"),
  compute_summary(list.files('./mammals/',        pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Mammals", "Scaffold"),
  compute_summary(list.files('./sharks/',         pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Cartilaginous Fishes", "Scaffold"),
  compute_summary(list.files('./fish/',           pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Ray-finned Fishes", "Scaffold"),
  compute_summary(list.files('./amphibians/',           pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE), "Amphibians", "Scaffold")
)

# ---- Compute contig summaries ----
contigs <- bind_rows(
  compute_summary(list.files('./reptiles/',       pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Lepidosauria", "Contig"),
  compute_summary(list.files('./cyclostomes/',    pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Cyclostomes", "Contig"),
  compute_summary(list.files('./lobe-finned_fishes/', pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Lobe-finned Fishes", "Contig"),
  compute_summary(list.files('./crocodiles/',     pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Crocodiles", "Contig"),
  compute_summary(list.files('./turtles/',        pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Turtles", "Contig"),
  compute_summary(list.files('./birds/',          pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Birds", "Contig"),
  compute_summary(list.files('./mammals/',        pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Mammals", "Contig"),
  compute_summary(list.files('./sharks/',         pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Cartilaginous Fishes", "Contig"),
  compute_summary(list.files('./fish/',           pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Ray-finned Fishes", "Contig"),
  compute_summary(list.files('./amphibians/',           pattern="*.contigs_Nx_table.tsv", full.names=TRUE), "Amphibians", "Contig")
)

# ---- Combine everything ----
plot_df <- bind_rows(scaffolds, contigs)

# ---- Plot ----
ggplot(plot_df, aes(x=Nx, y=log10(mean), color=Group, linetype=Type)) +
  geom_line() +
  scale_color_manual(values=my_colors) +
  scale_linetype_manual(values=c("Scaffold"="solid", "Contig"="dotted")) +
  labs(
    x = "Genome fraction",
    y = "log10(Sequence length)",
    linetype = "Assembly type"
  ) +
  ylim(c(1,9.5)) +
  theme_minimal()

ggplot(plot_df[plot_df$Type=="Scaffold",], aes(x=Nx, y=log10(mean), color=Group)) +
  geom_line() +
  scale_color_manual(values=my_colors) +
  #scale_linetype_manual(values=c("Scaffold"="solid", "Contig"="dotted")) +
  labs(
    x = "Genome fraction",
    y = "log10(Scaffold length)"
  ) +
  ylim(c(1,9.5)) +
  theme_minimal()

ggplot(plot_df[plot_df$Type=="Contig",], aes(x=Nx, y=log10(mean), color=Group)) +
  geom_line() +
  scale_color_manual(values=my_colors) +
  #scale_linetype_manual(values=c("Scaffold"="solid", "Contig"="dotted")) +
  labs(
    x = "Genome fraction",
    y = "log10(Contig length)"
  ) +
  ylim(c(1,9.5)) +
  theme_minimal()

plot_df_birds=plot_df[plot_df$Group=="Birds" & plot_df$Type=="Scaffold",]
ggplot(plot_df_birds, aes(x=Nx)) +
  geom_ribbon(aes(ymin=log10(q5), ymax=log10(q95)),
              fill="#00796B", alpha=0.3) +
  geom_line(aes(y=log10(mean)),
            color="#00796B", size=1) +
  labs(
    x = "Genome fraction",
    y = "log10(Scaffold length)",
    title = "Bird genomes: Mean Nx curve with 5–95% shaded"
  ) +
  theme_minimal()

ggplot(plot_df[plot_df$Type=="Scaffold" & plot_df$Group %in% c("Birds","Mammals","Ray-finned Fishes","Amphibians"),], aes(x=Nx, color=Group, fill=Group)) +
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

ggplot(plot_df[plot_df$Type=="Contig" & plot_df$Group %in% c("Birds","Mammals","Ray-finned Fishes","Amphibians"),], aes(x=Nx, color=Group, fill=Group)) +
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
  ylim(c(1,9.5)) +
  theme_minimal()


ggplot(plot_df[plot_df$Group %in% c("Birds"),], aes(x=Nx, color=Group, fill=Group,linetype=Type)) +
  # Shaded area between q5 and q95
  geom_ribbon(aes(ymin=log10(q5), ymax=log10(q95)), alpha=0.2, color=NA) +
  # Bold median line
  geom_line(aes(y=log10(median)), size=1) +
  scale_color_manual(values=my_colors) +
  scale_fill_manual(values=my_colors) +
  labs(
    x = "Genome fraction",
    y = "log10(Sequence length)",
    title = "Median Nx curves with 5–95% quantile shading"
  ) +
  ylim(c(3,8.5)) +
  theme_minimal()
  
