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

required_packages <- c("readxl", "data.table", "tidyverse", "ape", "ggplot2", 
                       "TreeTools", "ape", "devtools", "stringr", "ggnewscale", 
                       "BiocManager", "RColorBrewer", "googlesheets4")

install_load(required_packages)
lapply(required_packages, packageVersion)

BiocManager_packages <- c("treeio", "ggtree", "ggtreeExtra", "tidytree")
for (p in BiocManager_packages) {
  if(require(p, quietly=TRUE)){
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

# load VGP ordinal list
gs4_auth(token = NULL, scopes = "https://www.googleapis.com/auth/spreadsheets.readonly")
ordinal_list <- read_sheet(ss = "1zTqPvFvvsJqadnNDMvY-VO2xJSn6-92Y3La-jHwQ2fU", sheet = 1)
#ordinal_list <- read_excel("VGP Ordinal List.xlsx", sheet = 1) # ord download as Excel first
setDT(ordinal_list)
setnames(ordinal_list, c("Orders Scientific Name (inferred >50 MYA divergence times)"), c("order_50MYA"))
ordinal_list <- ordinal_list %>% mutate(order_50MYA = gsub(")", "", order_50MYA))
ordinal_list[, c("order", "suborder1", "suborder2", "suborder3") := tstrsplit(order_50MYA, "\\ \\(|\\ > |\\|", fixed=FALSE)]

#List all superorders:
write.table(as.data.frame(table(ordinal_list$Superorder[!is.na(ordinal_list$`Accession # for main haplotype`)])), ,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
#List all orders:
write.table(as.data.frame(table(ordinal_list$order[!is.na(ordinal_list$`Accession # for main haplotype`)])), ,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

# representation of ordinal species
ordinal_species <- ordinal_list %>% 
  mutate(complete_status = Status >= 4) %>%
  group_by(Order, complete_status) %>% 
  slice_head(n = 1) %>%
  ungroup %>%
  group_by(Order) %>% 
  mutate(count_complete_status = n()) %>%
  ungroup %>%
  filter((count_complete_status > 1 & complete_status) | count_complete_status == 1 | (count_complete_status > 1 & !complete_status))

# completeness stats
frequency <- as.data.frame(table(ordinal_species$complete_status))
frequency$Percent=round(frequency$Freq/sum(frequency$Freq)*100,2)
write.table(frequency, sep = "\t", row.names = FALSE, quote = FALSE)

# replace missing names in taxonomy with closely related species
conditions <- c("Ambystoma mexicanum x Ambystoma tigrinum", "Aspidoscelis tigris stejnegeri", "Aegotheles albertisi", "Osmerus mordax", "Hydrolagus colliei")
replacement_values <- c("Ambystoma mexicanum", "Aspidoscelis tigris", "Aegotheles albertisi albertisi", "Osmerus mordax mordax", "Hydrolagus bemisi")
inds <- match(ordinal_species$`Scientific Name`, conditions)
ordinal_species$`Scientific Name updated` <- ordinal_species$`Scientific Name` # Copy original column to a new one
ordinal_species$`Scientific Name updated`[!is.na(inds)] <- replacement_values[na.omit(inds)] # Replace only matched rows in the new column

# match names
VGP_ordinal_resolved_names <- rotl::tnrs_match_names(names = ordinal_species$`Scientific Name updated`)

# add status and lineage
VGP_ordinal_resolved_names_with_status <- cbind(VGP_ordinal_resolved_names, ordinal_species$complete_status, ordinal_species$Lineage, ordinal_species$`Scientific Name`)

#check presence in tree
in_tree <- rotl::is_in_tree(VGP_ordinal_resolved_names_with_status$ott_id)
VGP_ordinal_resolved_names_with_status[!in_tree,] # species missing in tree, this should be empty after the replacements above

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
completeness_grp <- list(missing   = unlist(VGP_ordinal_resolved_names_with_status[VGP_ordinal_resolved_names_with_status$`ordinal_species$complete_status` == FALSE,]$`ordinal_species$\`Scientific Name\``, use.names = FALSE),
            completed = unlist(VGP_ordinal_resolved_names_with_status[VGP_ordinal_resolved_names_with_status$`ordinal_species$complete_status` == TRUE,]$`ordinal_species$\`Scientific Name\``, use.names = FALSE))

# remove underscore from tip labels
VGP_ordinal_subtree$tip.label <- sapply(test, function(x) {
  i <- match(x, VGP_ordinal_resolved_names_with_status$tip)
  if (!is.na(i)) VGP_ordinal_resolved_names_with_status$`ordinal_species$\`Scientific Name\``[i] else x
}, USE.NAMES = FALSE)

# build metadata table for lineage
l <- VGP_ordinal_resolved_names_with_status[lengths(VGP_ordinal_resolved_names_with_status$`ordinal_species$\`Scientific Name\``)>0,] %>% select(`ordinal_species$Lineage`)
row.names(l) <- VGP_ordinal_resolved_names_with_status[lengths(VGP_ordinal_resolved_names_with_status$`ordinal_species$\`Scientific Name\``)>0,]$`ordinal_species$\`Scientific Name\``

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

class_colors <- c(
  Amphibians    = "#228B22",  # Forest green
  Birds         = "#DAA520",  # Goldenrod
  Fishes        = "#4682B4",  # Steel blue
  Invertebrates = "#DA70D6",  # Orchid
  Mammals       = "#B22222",  # Firebrick
  Reptiles      = "#556B2F"   # Dark olive green
)

p_ordinal2 <- p_ordinal1 + new_scale_fill() +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=lineage),
    width=2,
    offset=0
  ) +
  scale_fill_manual(
    name="Lineage",
    values=class_colors,
    guide=guide_legend(keywidth=0.3, keyheight=0.3, ncol=2, order=2)
  )+
  theme(
    legend.position = "right",
    legend.margin = margin(10, 10, 10, 10),
    legend.box.spacing = unit(1, "cm"),
    plot.margin = margin(40,60,40,20))
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
length(orders$ott_id)
length(families$ott_id)

# plot tree
if (!file.exists("vertebrate_tree_plot.rds")) {
  p<-ggtree(vertebrate_tree, layout="fan")
  saveRDS(p, file = "vertebrate_tree_plot.rds")
}else{
  p<-readRDS("vertebrate_tree_plot.rds")
}

# build metadata table for lineage
row_idx <- lapply(VGP_orders_resolved_names$unique_name, function(x) grep(x, vertebrate_tree$node.label))
internal_order_nodes <- as.list(vertebrate_tree$node.label[unlist(row_idx)])
subtrees <- sapply(internal_order_nodes, function(x) {extract.clade(vertebrate_tree, x)})
internal_order_nodes_grp <- as.list(subtrees['tip.label',])

# add order names to lists
names(internal_order_nodes_grp) <- internal_order_nodes

p1 <- groupOTU(p, completed_grp, 'status') + aes(color=status) +
  theme(legend.position="right",
        legend.margin=margin(0,0,0,40),
        legend.box.spacing = margin(4)) + 
  scale_color_manual(values = c("black","green"))

# Original levels (e.g. "Notacanthiformes ott925748", etc.)
lineage_levels <- internal_order_nodes

# Pretty labels for the legend only
legend_labels <- gsub(" ott\\d+$", "", lineage_levels)

# Name the colors with the full lineage values used in the data
names(order_colors) <- lineage_levels

p2 <- groupOTU(p1, internal_order_nodes_grp, 'lineage') + new_scale_fill() +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=lineage),
    width=2,
    offset=0
  ) +
  scale_fill_manual(
    name = "Lineage",
    values = order_colors,
    labels = legend_labels,
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3, ncol = 2, order = 2)
  )+
  theme(plot.margin = margin(40,0,40,0))
ggsave(dpi=600, filename='full_tree.png', width = 12, height = 8)
