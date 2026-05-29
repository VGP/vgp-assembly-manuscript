library('ggplot2')
library('matrixStats')
files <- list.files(pattern="*.scaffolds_Nx_table.tsv", full.names=TRUE, recursive=FALSE)

stats<-c()

for (f in files){
  stats_in<-read.delim(f,sep='\t',header = T)
  stats<-cbind(stats,stats_in[,2])
}
quantiles <- rowQuantiles(stats, probs=c(0.05, 0.95))
q5 <- quantiles[,1]
q95 <- quantiles[,2]


plot_df<-data.frame(Nx=c(1:1000)/1000,mean=rowMeans(stats),median=rowMedians(stats),
                      sd_min=rowMeans(stats)-rowSds(stats),sd_max=rowMeans(stats)+rowSds(stats),
                    q5 = q5,
                    q95 = q95)

ggplot(plot_df, aes(x=Nx)) +
  geom_line(aes(y=log10(median)), color="blue") +
  #geom_line(aes(y=log10(sd_min)), color="blue", linetype="dashed") +
  #geom_line(aes(y=log10(sd_max)), color="blue", linetype="dashed") +
  geom_line(aes(y=log10(q5)), color="blue", linetype="dashed") +
  geom_line(aes(y=log10(q95)), color="blue", linetype="dashed") +
  labs(x="Genome fraction", y="log10(Scaffold length)", title="Nx curve with 5%, 50% & 95% quantiles") +
  theme_minimal()


files <- list.files(pattern="*.contigs_Nx_table.tsv", full.names=TRUE, recursive=FALSE)

stats<-c()

for (f in files){
  stats_in<-read.delim(f,sep='\t',header = T)
  stats<-cbind(stats,stats_in[,2])
}
quantiles <- rowQuantiles(stats, probs=c(0.05, 0.95))
q5 <- quantiles[,1]
q95 <- quantiles[,2]

plot_df<-data.frame(Nx=c(1:1000)/1000,mean=rowMeans(stats),median=rowMedians(stats),
                    sd_min=rowMeans(stats)-rowSds(stats),sd_max=rowMeans(stats)+rowSds(stats),
                    q5 = q5,
                    q95 = q95)

ggplot(plot_df, aes(x=Nx)) +
  geom_line(aes(y=log10(median)), color="blue") +
  geom_line(aes(y=log10(q5)), color="blue", linetype="dashed") +
  geom_line(aes(y=log10(q95)), color="blue", linetype="dashed") +
  labs(x="Genome fraction", y="log10(Contig length)", title="Nx curve with 5%, 50% & 95% quantiles") +
  theme_minimal()
