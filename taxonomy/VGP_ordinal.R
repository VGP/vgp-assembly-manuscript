setwd(dirname(rstudioapi::getSourceEditorContext()$path))

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

required_packages <- c("readxl", "data.table", "tidyverse", "ape", "ggplot2", "TreeTools", "ape", "devtools", "stringr", "ggnewscale", "BiocManager")

install_load(required_packages)
lapply(required_packages, packageVersion)

BiocManager_packages <- c("treeio", "ggtree", "ggtreeExtra", "tidytree")
for (p in BiocManager_packages) {
  if(!require(p, quietly=TRUE)){
    BiocManager::install(p)
  }
}

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

# completeness stats
frequency <- as.data.frame(table(ordinal_species$complete_status))
frequency$Percent=frequency$Freq/sum(frequency$Freq)*100
frequency

# replace missing names in taxonomy with closely related species
conditions <- c("Ambystoma mexicanum x Ambystoma tigrinum", "Saccoglossus sp. HW-2024a", "Aspidoscelis tigris stejnegeri", "Polymixia hollisterae")
replacement_values <- c("Ambystoma mexicanum", "Saccoglossus kowalevskii", "Aspidoscelis tigris", "Polymixia berndti")
inds <- match(ordinal_species$`Scientific Name`, conditions)
ordinal_species$`Scientific Name`[!is.na(inds)] <- replacement_values[na.omit(inds)]

# match names
VGP_ordinal_resolved_names <- rotl::tnrs_match_names(names = ordinal_species$`Scientific Name`)

# add status and lineage
VGP_ordinal_resolved_names_with_status <- cbind(VGP_ordinal_resolved_names, ordinal_species$complete_status, ordinal_species$Lineage)

#check presence in tree
in_tree <- rotl::is_in_tree(VGP_ordinal_resolved_names_with_status$ott_id)
VGP_ordinal_resolved_names_with_status[!in_tree,] # species missing in tree

# induce tree from species list
VGP_ordinal_subtree <- rotl::tol_induced_subtree(VGP_ordinal_resolved_names_with_status$ott_id[in_tree])

# open tree taxonomy identifiers (ott_id) for completed species
completed_status <- VGP_ordinal_resolved_names_with_status$`ordinal_species$complete_status`
ott_id_completed <- VGP_ordinal_resolved_names_with_status$ott_id[completed_status]

# generate full tip names
VGP_ordinal_resolved_names_with_status <- VGP_ordinal_resolved_names_with_status%>%
  mutate(tip=sapply(as.character(ott_id), function(x) {
    str_subset(VGP_ordinal_subtree$tip.label, paste(".*_ott",x,"$", sep=""))}))

# group Operational Taxonomic Units (OTUs) by completeness
completeness_grp <- list(missing   = sub("_", " ", rotl::strip_ott_ids(unlist(VGP_ordinal_resolved_names_with_status[VGP_ordinal_resolved_names_with_status$`ordinal_species$complete_status` == FALSE,]$tip, use.names = FALSE))),
            completed = sub("_", " ", rotl::strip_ott_ids(unlist(VGP_ordinal_resolved_names_with_status[VGP_ordinal_resolved_names_with_status$`ordinal_species$complete_status` == TRUE,]$tip, use.names = FALSE))))

# remove underscore from tip labels
VGP_ordinal_subtree$tip.label <- sub("_", " ", rotl::strip_ott_ids(VGP_ordinal_subtree$tip.label))

# build metadata table
l <- VGP_ordinal_resolved_names_with_status[lengths(VGP_ordinal_resolved_names_with_status$tip)>0,] %>% select(`ordinal_species$Lineage`)
row.names(l) <- sub("_", " ", rotl::strip_ott_ids(VGP_ordinal_resolved_names_with_status[lengths(VGP_ordinal_resolved_names_with_status$tip)>0,]$tip))

metadata <- data.frame (
  label = row.names(l),
  lineage = l$`ordinal_species$Lineage`
)

row.names(metadata) <- row.names(l)

lineage_colors <- metadata %>%
  select('lineage') %>%
  distinct()

metadata$lineage <- factor(metadata$lineage, levels=lineage_colors$lineage)

p_ordinal <- ggtree(VGP_ordinal_subtree, layout="fan") + 
  geom_tiplab(size = 2, offset=2)

p_ordinal <- p_ordinal %<+% metadata

p_ordinal1 <- groupOTU(p_ordinal, completeness_grp, 'status') + aes(color=status) +
  theme(legend.position="right",
        legend.margin=margin(0,0,0,40),
        legend.box.spacing = margin(4)) + 
  scale_color_manual(values = c("black", "red"))

p_ordinal2 <- p_ordinal1 + new_scale_fill() +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=lineage),
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

###### Full vertebrate tree

#import vertebrate tree
if (!file.exists("vertebrate_tree.rds")) {
  id <- 801601
  vertebrate_tree <- rotl::tol_subtree(ott_id = id)
  saveRDS(vertebrate_tree, file = "vertebrate_tree.rds")
} else {
  vertebrate_tree <- readRDS("vertebrate_tree.rds")
}
# remove underscore from tip labels
vertebrate_tree$tip.label <- sub("_", " ", rotl::strip_ott_ids(vertebrate_tree$tip.label))

if (!file.exists("taxa.rds")) {
  orders <- datelife::get_ott_children(ott_ids = 801601, ott_rank = "order")
  orders <- orders[orders$rank == 'order',]
  families <- datelife::get_ott_children(ott_ids = 801601, ott_rank = "family")
  #species <- datelife::get_ott_children(ott_ids = 801601, ott_rank = "species")
  saveRDS(list(orders, families), file = "taxa.rds")
} else {
  taxa <- readRDS("taxa.rds")
  orders <- taxa[[1]]$Vertebrata
  families <- taxa[[2]]$Vertebrata
}

# find orders
VGP_orders_resolved_names <- rotl::tnrs_match_names(names = row.names(orders))

#check presence in tree
orders_in_tree <- rotl::is_in_tree(VGP_orders_resolved_names$ott_id)
VGP_orders_resolved_names[!orders_in_tree,] # species missing in tree

# completed species
completed_grp <- list(completed = sub("_", " ", rotl::strip_ott_ids(unlist(VGP_ordinal_resolved_names_with_status[VGP_ordinal_resolved_names_with_status$`ordinal_species$complete_status` == TRUE,]$tip, use.names = FALSE))))

# Statistics:
length(vertebrate_tree$tip.label)
length(orders$Vertebrata[orders$Vertebrata$rank == 'order',]$ott_id)
length(families$Vertebrata[families$Vertebrata$rank == 'family',]$ott_id)

# plot tree
if (!file.exists("vertebrate_tree_plot.rds")) {
  p<-ggtree(vertebrate_tree, layout="fan")
  saveRDS(vertebrate_tree, file = "vertebrate_tree_plot.rds")
}else{
  readRDS("vertebrate_tree_plot.rds")
}

# add VGP metadata
p <- p %<+% metadata

p1 <- groupOTU(p, completed_grp, 'status') + aes(color=status) +
  theme(legend.position="right",
        legend.margin=margin(0,0,0,40),
        legend.box.spacing = margin(4)) + 
  scale_color_manual(values = c("black","green"))

p2 <- p1 + new_scale_fill() +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=lineage),
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


ggtree(vertebrate_tree, layout="fan") + 
  theme_tree2()
