setwd(dirname(rstudioapi::getSourceEditorContext()$path))

requiredpackages <- c("readxl", "data.table", "tidyverse", "ape", "ggplot2", "TreeTools", "ape", "devtools", "stringr", "ggnewscale")

install_load <- function(packages){
  for (p in packages) {
    if (p %in% rownames(installed.packages())) {
      library(p, character.only=TRUE)
    } else {
      install.packages(p, repos = "http://cran.us.r-project.org")
      library(p,character.only = TRUE)
    }
  }
}

install_load(requiredpackages)
lapply(requiredpackages, packageVersion)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("treeio")
BiocManager::install("ggtree")
BiocManager::install("ggtreeExtra")
BiocManager::install("tidytree")
library(treeio)
library(ggtree)
library(ggtreeExtra)
library(tidytree)

library(remotes)
install_github("ropensci/bold") # needed by taxize
install_github("ropensci/taxize") # needed by datelife
remotes::install_github("ropensci/rotl") # default version is buggy
library(rotl)

devtools::install_github("phylotastic/datelife")
devtools::install_github("phylotastic/datelifeplot")

#load VGP ordinal list (https://docs.google.com/spreadsheets/d/1Jwjv6Kwc6VIn1UMMhnG6kvFCxjwGdC5b7p_HtbDOMOs)
ordinal_list <- read_excel("VGP Ordinal List.xlsx", sheet = 1)
setDT(ordinal_list)
setnames(ordinal_list, c("Orders Scientific Name (inferred >50 MYA divergence times)"), c("order_50MYA"))
ordinal_list <- ordinal_list %>% mutate(order_50MYA = gsub(")", "", order_50MYA))
ordinal_list[, c("order", "suborder1", "suborder2", "suborder3", "suborder4") := tstrsplit(order_50MYA, "\\ \\(|\\ > |\\|", fixed=FALSE)]

#List all superorders:
as.data.frame(table(ordinal_list$Superorder))
#List all orders:
as.data.frame(table(ordinal_list$order))

# representation of ordinal species
ordinal_species <- ordinal_list %>% 
  mutate(complete_status = Status == 4) %>%
  group_by(Order, complete_status) %>% 
  slice_head(n = 1) %>%
  ungroup %>%
  group_by(Order) %>% 
  mutate(count_complete_status = n()) %>%
  ungroup %>%
  filter((count_complete_status > 1 & complete_status) | count_complete_status == 1)

conditions <- c("Ambystoma mexicanum x Ambystoma tigrinum", "Saccoglossus sp. HW-2024a", "Aspidoscelis tigris stejnegeri", "Polymixia hollisterae")
replacement_values <- c("Ambystoma mexicanum", "Saccoglossus kowalevskii", "Aspidoscelis tigris", "Polymixia_berndti")
inds <- match(ordinal_species$`Scientific Name`, conditions)
ordinal_species$`Scientific Name`[!is.na(inds)] <- replacement_values[na.omit(inds)]

VGP_ordinal_resolved_names <- rotl::tnrs_match_names(names = ordinal_species$`Scientific Name`)

VGP_ordinal_resolved_names_with_status <- cbind(VGP_ordinal_resolved_names, ordinal_species$complete_status, ordinal_species$Lineage)
in_tree <- rotl::is_in_tree(VGP_ordinal_resolved_names_with_status$ott_id)
VGP_ordinal_resolved_names_with_status[!in_tree,] # species missing in tree
VGP_ordinal_subtree <- rotl::tol_induced_subtree(VGP_ordinal_resolved_names_with_status$ott_id[in_tree])

frequency <- as.data.frame(table(ordinal_species$complete_status))
frequency$Percent=frequency$Freq/sum(frequency$Freq)*100
frequency

completed_status <- VGP_ordinal_resolved_names_with_status$`ordinal_species$complete_status`
ott_id_completed <- VGP_ordinal_resolved_names_with_status$ott_id[completed_status]

VGP_ordinal_resolved_names_with_status <- VGP_ordinal_resolved_names_with_status%>%
  mutate(tip=sapply(as.character(ott_id), function(x) {
    str_subset(VGP_ordinal_subtree$tip.label, paste(".*_ott",x,"$", sep=""))}))

completeness_grp <- list(missing   = sub("_", " ", strip_ott_ids(unlist(VGP_ordinal_resolved_names_with_status[VGP_ordinal_resolved_names_with_status$`ordinal_species$complete_status` == FALSE,]$tip, use.names = FALSE))),
            completed = sub("_", " ", strip_ott_ids(unlist(VGP_ordinal_resolved_names_with_status[VGP_ordinal_resolved_names_with_status$`ordinal_species$complete_status` == TRUE,]$tip, use.names = FALSE))))

VGP_ordinal_subtree$tip.label <- sub("_", " ", strip_ott_ids(VGP_ordinal_subtree$tip.label))

l <- VGP_ordinal_resolved_names_with_status[lengths(VGP_ordinal_resolved_names_with_status$tip)>0,] %>% select(`ordinal_species$Lineage`)
row.names(l) <- sub("_", " ", strip_ott_ids(VGP_ordinal_resolved_names_with_status[lengths(VGP_ordinal_resolved_names_with_status$tip)>0,]$tip))

metadata <- data.frame (
  label = row.names(l),
  lineage = l
)

lineage_colors <- metadata %>%
  select('ordinal_species.Lineage') %>%
  distinct()

metadata$lineage <- factor(metadata$ordinal_species.Lineage, levels=lineage_colors$ordinal_species.Lineage)

p<-ggtree(VGP_ordinal_subtree, layout="fan") + 
  geom_tiplab(size = 2, offset=2)

p <- p %<+% metadata

p1 <- groupOTU(p, completeness_grp, 'status') + aes(color=status) +
  theme(legend.position="right",
        legend.margin=margin(0,0,0,40),
        legend.box.spacing = margin(4)) + 
  scale_color_manual(values = c("black", "red"))

p2 <- p1 + new_scale_fill() +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=ordinal_species.Lineage),
    width=2,
    offset=0
  ) +
  scale_fill_manual(
    name="Lineage",
    values=c(1:6),
    guide=guide_legend(keywidth=0.3, keyheight=0.3, ncol=2, order=2)
  )+
  theme(plot.margin = margin(40,0,40,0))
ggsave(dpi=600, filename='ordinal_tree.png', width = 12, height = 8)

#import vertebrate tree
id <- 801601
vertebrate_tree <- rotl::tol_subtree(ott_id = id)
ggtree(vertebrate_tree, layout="fan") + 
  theme_tree2()
