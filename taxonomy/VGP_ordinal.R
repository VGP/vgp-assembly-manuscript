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
                       "BiocManager", "RColorBrewer", "googlesheets4", "ggforce", "gtable", "ggpubr", "grid", "systemfonts", "svglite")

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
library(dplyr)
library(ggforce)

library(remotes)
install_github("ropensci/bold") # needed by taxize
install_github("ropensci/taxize") # needed by datelife
remotes::install_github("ropensci/rotl") # default version is buggy
library(rotl)

devtools::install_github("phylotastic/datelife")
devtools::install_github("phylotastic/datelifeplot")

# load VGP ordinal list
gs4_auth(token = NULL, scopes = "https://www.googleapis.com/auth/spreadsheets.readonly", email = "giulio.formenti@gmail.com")
ordinal_list <- read_sheet(ss = "17aOjpVgclwdDcDccx7Jpuy4IL6PUcVYntB-iwfj0aZ4", sheet = 1)

#ordinal_list <- read_excel("VGP Ordinal List.xlsx", sheet = 1) # ord download as Excel first
setDT(ordinal_list)
setnames(ordinal_list, c("Orders Scientific Name (inferred >50 MYA divergence times)"), c("order_50MYA"))
ordinal_list <- ordinal_list %>% mutate(order_50MYA = gsub(")", "", order_50MYA))
ordinal_list[, c("order", "suborder1", "suborder2", "suborder3") := tstrsplit(order_50MYA, "\\ \\(|\\ > |\\|", fixed=FALSE)]

# fix the hybrid
ordinal_list$`Scientific Name`[ordinal_list$`Scientific Name` == "Ambystoma mexicanum x Ambystoma tigrinum"] <- "Ambystoma mexicanum"

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

# Count sequenced species (Status >= 4) per Order from the full table
order_seq_counts <- ordinal_list %>%
  transmute(
    Order,
    species = `Scientific Name`,
    is_seq = !is.na(`Accession # for main haplotype`)
  ) %>%
  distinct(Order, species, .keep_all = TRUE) %>%   # avoid duplicate species rows, if any
  group_by(Order) %>%
  summarise(
    n_sequenced = sum(is_seq, na.rm = TRUE),
    n_species   = n(),                              # total species listed for that Order
    pct_seq     = n_sequenced / n_species,
    .groups = "drop"
  )

# Attach to existing ordinal_species dataframe
ordinal_species <- ordinal_species %>%
  left_join(order_seq_counts, by = "Order")

# completeness stats
frequency <- as.data.frame(table(ordinal_species$complete_status))
frequency$Percent=round(frequency$Freq/sum(frequency$Freq)*100,2)
write.table(frequency, sep = "\t", row.names = FALSE, quote = FALSE)

# replace missing names in taxonomy with closely related species
conditions <- c("Aspidoscelis tigris stejnegeri", "Aegotheles albertisi", "Osmerus mordax", "Hydrolagus colliei")
replacement_values <- c("Aspidoscelis tigris", "Aegotheles albertisi albertisi", "Osmerus mordax mordax", "Hydrolagus bemisi")
inds <- match(ordinal_species$`Scientific Name`, conditions)
ordinal_species$`Scientific Name updated` <- ordinal_species$`Scientific Name` # Copy original column to a new one
ordinal_species$`Scientific Name updated`[!is.na(inds)] <- replacement_values[na.omit(inds)] # Replace only matched rows in the new column

# match names
VGP_ordinal_resolved_names <- rotl::tnrs_match_names(names = ordinal_species$`Scientific Name updated`)

# add status and extended lineage
VGP_ordinal_resolved_names_with_status <- cbind(VGP_ordinal_resolved_names, ordinal_species$complete_status, ordinal_species$`Extended lineage`, ordinal_species$`Scientific Name`)

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
VGP_ordinal_subtree$tip.label <- sapply(VGP_ordinal_subtree$tip.label, function(x) {
  i <- match(x, VGP_ordinal_resolved_names_with_status$tip)
  if (!is.na(i)) VGP_ordinal_resolved_names_with_status$`ordinal_species$\`Scientific Name\``[i] else x
}, USE.NAMES = FALSE)

# build metadata table for lineage
l <- VGP_ordinal_resolved_names_with_status[lengths(VGP_ordinal_resolved_names_with_status$`ordinal_species$\`Scientific Name\``)>0,] %>% select(`ordinal_species$\`Extended lineage\``)
row.names(l) <- VGP_ordinal_resolved_names_with_status[lengths(VGP_ordinal_resolved_names_with_status$`ordinal_species$\`Scientific Name\``)>0,]$`ordinal_species$\`Scientific Name\``

metadata <- data.frame (
  label = row.names(l),
  lineage = l$`ordinal_species$\`Extended lineage\``
  )

lineage_colors <- metadata %>%
  select('lineage') %>%
  distinct()

metadata$order = ordinal_list$`Order (NCBI)`[match(metadata$label, ordinal_list$`Scientific Name`)]
metadata$assembly_size = ordinal_list$`Assembly Size`[match(metadata$label, ordinal_list$`Scientific Name`)]
metadata$total_species_order = ordinal_list$`# species/order`[match(metadata$label, ordinal_list$`Scientific Name`)]
metadata$sequenced_species_order = ordinal_species$n_sequenced[match(metadata$label, ordinal_species$`Scientific Name`)]
metadata$completed = VGP_ordinal_resolved_names_with_status$`ordinal_species$complete_status`[match(metadata$label, VGP_ordinal_resolved_names_with_status$`ordinal_species$\`Scientific Name\``)]
metadata$label_to_plot = ifelse(metadata$completed, metadata$label, NA)
metadata$lineage <- factor(metadata$lineage, levels=lineage_colors$lineage)

class_colors <- c(
  "Mammals" = "#E69F00",  # Okabe–Ito Orange
  "Birds" = "#00796B",    # Teal 700
  "Crocodilians" = "#009688",# Teal 500
  "Turtles" = "#4DB6AC",  # Teal 300
  "Lepidosauria" = "#80CBC4", # Teal 200
  "Amphibians" = "#984EA3",   # Purple
  "Lobe-finned fishes" = "#A6761D", # Brown
  "Ray-finned fishes" = "#56B4E9",  # Okabe–Ito Sky Blue
  "Cartilaginous fishes" = "#0072B2", # Okabe–Ito Blue
  "Cyclostomes" = "#CC79A7", # Okabe–Ito Magenta
  "Other Deurostomes" = "#999999"  # Neutral Gray
)

# 1) Choose which lineages go on row 1 vs row 2
top_row <- c("Mammals","Birds","Crocodilians","Turtles","Lepidosauria","Amphibians")
bottom_row <- c("Lobe-finned fishes",
                "Ray-finned fishes","Cartilaginous fishes",
                "Cyclostomes","Other Deurostomes")

desired_order <- c(top_row, bottom_row)

# Apply the order to the factor used by the plot
metadata$lineage <- factor(metadata$lineage, levels = desired_order)

# 2) Two-line legend labels (insert \n where you want the wrap and adjust wording here)
display_labels <- c(
  "Mammals" = "Mammals",
  "Birds" = "Birds",
  "Crocodilians" = "Crocodilians",
  "Turtles" = "Turtles",
  "Lepidosauria" = "Lepidosauria\n(Squamates and Tuatara)",
  "Amphibians" = "Amphibians",
  "Lobe-finned fishes" = "Lobe-finned fishes\n(coelacanth/lungfish)",
  "Ray-finned fishes" = "Ray-finned fishes",
  "Cartilaginous fishes" = "Cartilaginous fishes",
  "Cyclostomes" = "Cyclostomes (Jawless fishes)",
  "Other Deurostomes" = "Other Deurostomes\n(non-vertebrates)"
)

row.names(metadata) <- row.names(l)

# Make "Genus species subspecies" -> "G. species"
# and "Genus1 sp1 x Genus2 sp2"   -> "G. sp1 × G. sp2"
shorten_species_label <- function(x) {
  normalize <- function(s) gsub("[[:space:]]+", " ", trimws(s))
  abbr_one <- function(s) {
    s <- normalize(s)
    if (s == "" || is.na(s)) return(NA_character_)
    toks <- strsplit(s, " ", fixed = TRUE)[[1]]
    if (length(toks) == 1L) return(s)                # leave one-word names alone
    genus_abbr <- paste0(substr(toks[1], 1, 1), ".") # G.
    species    <- toks[2]                            # keep only species epithet
    paste(genus_abbr, species)
  }
  vapply(x, function(s) {
    if (is.na(s)) return(NA_character_)
    s <- normalize(s)
    parts <- strsplit(s, "(?i)\\s+(?:x|×)\\s+", perl = TRUE)[[1]]
    paste(vapply(parts, abbr_one, character(1)), collapse = " × ")
  }, character(1))
}

# Use it for your plotted labels (only for completed species, as before)
metadata$label_to_plot_short <- ifelse(
  metadata$completed,
  shorten_species_label(metadata$label),
  NA_character_
)

# Prepare inward (negative) values and shared limits
metadata$sequenced_inward <- -metadata$sequenced_species_order
xmax <- max(
  metadata$total_species_order,
  metadata$sequenced_species_order,
  na.rm = TRUE
)

# Create a scatterplot
plot(metadata$sequenced_species_order, metadata$total_species_order,
     main = "Scatterplot of X vs Y", # Plot title
     xlab = "Sequenced species per order",            # X-axis label
     ylab = "Total species per order",            # Y-axis label
     pch = 19,                       # Point style (filled circles)
     col = "blue",
     log = "yx")                   # Point color

# choose readable ticks
axis_breaks <- c(0,2000,4000)
axis_labs   <- c(0, "2K", "4K")  # show magnitudes only

#### plotting starts here ####

# Start from the circular tree
p_ordinal <- ggtree(VGP_ordinal_subtree, layout = "circular") %<+% metadata

# (1) INNER ring: lineage tiles
p_ordinal <- groupOTU(p_ordinal, completeness_grp, 'status') + aes(color=status) + new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(fill = lineage),
    width = 2,          # ring thickness
    offset = 0        # distance from tips
  )

# (2) ONE ultra-compact back-to-back ring
# Outward = total species per order (positive)
p_ordinal <- p_ordinal +
  geom_fruit(
    geom        = geom_col,
    mapping     = aes(x = total_species_order),
    orientation = "y",
    width       = 0.4,     # ring thickness (keep small)
    offset      = 0.7,     # distance from tips
    axis.params = list(
      axis       = "x",
      text.size  = 2
    ),
    fill="darkgrey"
  )

# Inward = sequenced species per order (negative)
# Same offset/width/limits, hide axis to avoid duplication
p_ordinal <- p_ordinal +
  geom_fruit(
    geom        = geom_col,
    mapping     = aes(x = sequenced_inward),
    orientation = "y",
    width       = 0.4,
    offset      = -0.2,
    axis.params = list(axis = "x", text.size = 2),
    fill="lightgrey"
  )

# (3) Labels outside both rings
p_ordinal <- p_ordinal + 
  geom_tiplab2(aes(label = label_to_plot_short, angle = angle), size = 2, offset = 12, fontface = "italic")

## orders
# --- 1) Tip -> Order data taken ONLY from the plotted tree ---
tip_order <- p_ordinal$data %>%
  filter(isTip) %>%
  transmute(
    order = !!sym("order"),                 # order must already be in metadata you %<+% attached
    label = trimws(as.character(label)),    # exact tip labels used by the plot
    angle = as.numeric(angle)
  )

tip_order <- tip_order %>%
  mutate(order = ifelse(is.na(order), paste0("U_"), order))

# --- 2) Helper: build one strip row for a single order ---
build_order_strip <- function(df_one_order) {
  df <- df_one_order %>% dplyr::filter(!is.na(angle))
  if (nrow(df) == 0) return(NULL)
  
  a <- (df$angle %% 360)
  o <- order(a)
  a <- a[o]
  labs <- trimws(as.character(df$label[o]))
  
  # Match ggtree's vertical flip rule
  make_upright <- function(theta) {
    theta <- theta %% 360
    is_left <- theta > 90 & theta < 270
    angle_text <- if (is_left) theta - 180 else theta
    hjust_val  <- if (is_left) 1 else 0
    list(angle_text = angle_text, hjust_val = hjust_val)
  }
  
  if (length(a) == 1L) {
    adj <- make_upright(a[1])
    return(tibble::tibble(
      order      = df$order[1],
      taxa1      = labs[1],
      taxa2      = labs[1],
      angle_text = adj$angle_text,
      hjust_val  = adj$hjust_val,
      theta      = a[1]                 # use the tip angle for ordering
    ))
  }
  
  # robust midpoint over the circular boundary
  span <- (a[length(a)] - a[1]) %% 360
  mid_deg <- (a[1] + span / 2) %% 360
  adj <- make_upright(mid_deg)
  
  tibble::tibble(
    order      = df$order[1],
    taxa1      = labs[1],
    taxa2      = if (identical(labs[1], labs[length(labs)]) && length(labs) >= 2) labs[2] else labs[length(labs)],
    angle_text = adj$angle_text,
    hjust_val  = adj$hjust_val,
    theta      = mid_deg               # center angle for alternating around the circle
  )
}

# --- 3) Build the strips table (one row per order) + alternating color ---
order_strips <- tip_order %>%
  dplyr::group_split(order, .keep = TRUE) %>%
  lapply(build_order_strip) %>%
  dplyr::bind_rows() %>%
  dplyr::filter(!is.na(taxa1), !is.na(taxa2)) %>%
  dplyr::arrange(theta) %>%
  dplyr::mutate(
    idx         = dplyr::row_number(),
    strip_color = ifelse(idx %% 2 == 1, "gray", "black")
  )

# Color orders red ONLY if none of their species is represented among the plotted (completed) tips
order_rep <- p_ordinal$data %>%
  dplyr::filter(isTip) %>%
  dplyr::group_by(order = !!rlang::sym("order")) %>%
  dplyr::summarise(n_rep = sum(!is.na(label_to_plot_short)), .groups = "drop") %>%  # completed & labeled
  dplyr::mutate(has_rep = n_rep > 0)

order_strips <- order_strips %>%
  dplyr::left_join(order_rep, by = "order") %>%
  dplyr::mutate(status_color = ifelse(is.na(has_rep) | !has_rep, "red", "black"))

# quick diagnostic
cat("Orders in tip_order:", n_distinct(tip_order$order),
    " | Orders with strips:", n_distinct(order_strips$order), "\n")

# --- 4) Draw: loop geom_strip so text actually appears ---
p_ordinal1 <- p_ordinal
for (i in seq_len(nrow(order_strips))) {
  p_ordinal1 <- p_ordinal1 +
    geom_striplab(
      taxa1   = order_strips$taxa1[i],     # exact tip labels from p_tree$data
      taxa2   = order_strips$taxa2[i],
      label   = order_strips$order[i],
      angle   = order_strips$angle_text[i],
      hjust   = order_strips$hjust_val[i],
      offset  = 1,                         # tweak relative to your tiplab offset (10)
      offset.text = 1,
      align   = TRUE,
      barsize = 1,                         # set >0 for a thin connector
      extend = 0,
      fill    = NA,                        # or e.g. "black", alpha = 0.1 for a band
      barcolor = "black",  # Color for the bar, order_strips$strip_color[i]
      textcolor = order_strips$status_color[i],   # Color for the label text
      fontsize = 2
    )
}

# keep the tip label size in one place
tip_font_mm <- 2     # <- this matches geom_tiplab2(size = 2)
mm_to_pt    <- function(mm) mm / 0.3527778
legend_pt   <- mm_to_pt(tip_font_mm)   # ~ 5.67 pt

# legends
p_ordinal2 <- p_ordinal1 +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.margin = margin(t = 4, r = 0, b = 0, l = 0),
    legend.box.spacing = unit(4, "mm"),
    legend.text  = element_text(size = legend_pt),
    legend.title = element_text(size = legend_pt)
  ) +
  scale_color_manual(
    values = c("black", "red"),
    guide = guide_legend(title = "Status", nrow = 1, byrow = TRUE, order = 1, title.theme = element_text(size = legend_pt),
                         label.theme = element_text(size = legend_pt))
  ) +
  scale_fill_manual(
    name   = "Lineage",
    values = class_colors[desired_order],
    breaks = desired_order,                      # order in legend
    labels = display_labels[desired_order],      # two-line text
    guide  = guide_legend(nrow = 2, ncol = length(top_row), byrow = TRUE, keywidth = 0.5, keyheight = 0.5,
                          title.theme = element_text(size = legend_pt),
                          label.theme = element_text(size = legend_pt))
  )

ggsave(filename='ordinal_tree.svg', width = 15, height = 12, device = "svg")

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
