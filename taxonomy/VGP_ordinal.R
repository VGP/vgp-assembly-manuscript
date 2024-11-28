setwd(dirname(rstudioapi::getSourceEditorContext()$path))

requiredpackages <- c("readxl", "data.table", "tidyverse", "ape", "ggplot2", "TreeTools", "ape", "devtools", "stringr")

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

#import vertebrate tree
ids <- 801601
vertebrate_tree <- datelife::get_ott_children(ott_ids = ids, ott_rank = "species")
vertebrate_tree %>% 
  ggtree() + 
  geom_tiplab() +
  theme_tree2()

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

VGP_ordinal_resolved_names <- rotl::tnrs_match_names(names = ordinal_species$`Scientific Name`)
VGP_ordinal_resolved_names <- cbind(VGP_ordinal_resolved_names, ordinal_species$complete_status)

VGP_ordinal_in_tree <- rotl::is_in_tree(ott_id(VGP_ordinal_resolved_names))
VGP_ordinal_subtree <- rotl::tol_induced_subtree(VGP_ordinal_resolved_names$ott_id[VGP_ordinal_in_tree])

completed_status <- VGP_ordinal_resolved_names$`ordinal_species$complete_status`
names(completed_status) <- VGP_ordinal_subtree$node.label
VGP_ordinal_subtree_grouped <- groupOTU(VGP_ordinal_subtree, completed_status, "completed")

frequency <- as.data.frame(table(completed_status))
frequency$Percent=frequency$Freq/sum(frequency$Freq)*100
frequency

p <- ggtree(VGP_ordinal_subtree_grouped, mapping=aes(color=completed), layout="fan") + 
#  geom_tiplab() + 
  theme_tree2() + 
  scale_color_manual(values = c("red", "black"))

p

