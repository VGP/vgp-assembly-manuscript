library('stringr')
library('ggplot2')

plot_df<-read.table(file = 'Percent_assignedChromosome.tsv',header = T,sep = '\t')
p<-ggplot(
  data = plot_df,
  aes(
    y = as.numeric(perc),
    x = Project,
    fill = Project,
    color = Project
  )
) +
  geom_jitter(
    width = 0.4,
    height = 0,
    alpha = 0.3,
    size = 0.5
  ) +
  geom_violin(
    alpha = 0.35,        # keep fill but make it restrained
    width = 0.7,
    scale = "area",
    linewidth = 0.5      # slightly stronger outline
  ) +
  geom_boxplot(
    fill = NA,           # box stays unfilled
    width = 0.2,
    linewidth = 1.1,     # clearly heavier box + whiskers
    outlier.shape = NA
  ) +
  geom_hline(
    yintercept = 90,
    linetype = "dashed",
    color = "darkgrey"
  ) +
  scale_y_continuous(limits = c(0.000001, 100)) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  ylab("Genome Assigned to Chromosomes (%)") +
  xlab("") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 14, face = "bold")
  )
p
ggsave(filename = 'Fig2d_chromosomeAssigned.svg',plot=p,dpi=1200)
