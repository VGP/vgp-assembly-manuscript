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

bioc_packages <- c(
  "treeio",
  "ggtree",
  "ggtreeExtra",
  "tidytree"
)

for (p in bioc_packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
  library(p, character.only = TRUE)
}

library(dplyr)
library(ggforce)

# ---------------------------------------------------------------------------
# GitHub-only dependencies required for OpenTree descendant retrieval
# ---------------------------------------------------------------------------

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("bold", quietly = TRUE)) {
  remotes::install_github(
    "ropensci/bold",
    upgrade = "never",
    dependencies = TRUE
  )
}

if (!requireNamespace("datelife", quietly = TRUE)) {
  remotes::install_github(
    "phylotastic/datelife",
    upgrade = "never",
    dependencies = TRUE
  )
}

library(datelife)

github_packages <- c(
  bold         = "ropensci/bold",
  taxize       = "ropensci/taxize",
  rotl         = "ropensci/rotl",
  datelife     = "phylotastic/datelife",
  datelifeplot = "phylotastic/datelifeplot"
)

for (p in names(github_packages)) {
  if (!requireNamespace(p, quietly = TRUE)) {
    remotes::install_github(github_packages[[p]], upgrade = "never")
  }
  library(p, character.only = TRUE)
}

# load VGP ordinal list
gs4_auth(token = NULL, scopes = "https://www.googleapis.com/auth/spreadsheets.readonly", email = "giulio.formenti@gmail.com", cache = TRUE)
ordinal_list <- read_sheet(ss = "118_VJMfdvHzZfcJPIN4jJwSZ157TjtoAwfMVIAIYSmI", sheet = 1)

#ordinal_list <- read_excel("VGP Ordinal List.xlsx", sheet = 1) # ord download as Excel first
setDT(ordinal_list)

# April 2026 spreadsheet: updated header for the inferred >50 MYA order field
setnames(
  ordinal_list,
  "Inferred Orders (>50 MYA divergence times)",
  "order_50MYA"
)

ordinal_list <- ordinal_list %>%
  mutate(order_50MYA = gsub(")", "", order_50MYA, fixed = TRUE))

ordinal_list[, c("order", "suborder1", "suborder2", "suborder3") :=
               tstrsplit(order_50MYA, "\\ \\(|\\ > |\\|", fixed = FALSE)]

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

# representation of ordinal species (VGP Orders)
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

# representation of official NCBI orders by completed assemblies
ncbi_order_completion <- ordinal_list %>%
  filter(!is.na(`Order (NCBI)`)) %>%
  group_by(`Order (NCBI)`) %>%
  summarise(
    complete_status = any(Status >= 4, na.rm = TRUE),
    .groups = "drop"
  )

frequency_ncbi <- as.data.frame(table(ncbi_order_completion$complete_status))
frequency_ncbi$Percent <- round(
  frequency_ncbi$Freq / sum(frequency_ncbi$Freq) * 100,
  2
)

write.table(
  frequency_ncbi,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# replace missing names in taxonomy with closely related species
conditions <- c("Aspidoscelis tigris stejnegeri", "Aegotheles albertisi", "Osmerus mordax", "Hydrolagus colliei", "Bassozetus sp. 2 HX-2024")
replacement_values <- c("Aspidoscelis tigris", "Aegotheles albertisi albertisi", "Osmerus mordax mordax", "Hydrolagus bemisi", "Bassozetus zenkevitchi")
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

# ---------------------------------------------------------------------------
# Match OpenTree tips to VGP species, then use clean VGP names for plotting
# ---------------------------------------------------------------------------

tip_map_ordinal <- tibble::tibble(
  raw_tip = VGP_ordinal_subtree$tip.label,
  ott_id = as.numeric(
    stringr::str_extract(raw_tip, "(?<=ott)[0-9]+$")
  )
)

if (any(is.na(tip_map_ordinal$ott_id))) {
  print(tip_map_ordinal %>% dplyr::filter(is.na(ott_id)))
  stop("Could not recover OTT identifiers from one or more VGP_ordinal_subtree tips.")
}

# Add raw OpenTree tip identifier to the species table
VGP_ordinal_resolved_names_with_status_in_tree <-
  VGP_ordinal_resolved_names_with_status %>%
  dplyr::inner_join(tip_map_ordinal, by = "ott_id") %>%
  dplyr::mutate(
    plot_tip = `ordinal_species$\`Scientific Name\``
  )

# Map current raw tree tips to the original VGP scientific names
tip_label_map <- VGP_ordinal_resolved_names_with_status_in_tree %>%
  dplyr::select(raw_tip, plot_tip)

new_tip_labels <- tip_label_map$plot_tip[
  match(VGP_ordinal_subtree$tip.label, tip_label_map$raw_tip)
]

if (any(is.na(new_tip_labels))) {
  stop("Some OpenTree tip labels could not be matched back to VGP scientific names.")
}

# From this point onward, the plotted tree uses clean scientific names
VGP_ordinal_subtree$tip.label <- new_tip_labels

# groupOTU now uses the clean labels that are actually in the plotted tree
completeness_grp <- list(
  missing = VGP_ordinal_resolved_names_with_status_in_tree$plot_tip[
    !VGP_ordinal_resolved_names_with_status_in_tree$`ordinal_species$complete_status`
  ],
  completed = VGP_ordinal_resolved_names_with_status_in_tree$plot_tip[
    VGP_ordinal_resolved_names_with_status_in_tree$`ordinal_species$complete_status`
  ]
)

stopifnot(
  all(unlist(completeness_grp) %in% VGP_ordinal_subtree$tip.label)
)

# ---------------------------------------------------------------------------
# Build metadata keyed by the clean plotted tip labels
# ---------------------------------------------------------------------------

metadata <- VGP_ordinal_resolved_names_with_status_in_tree %>%
  dplyr::transmute(
    label = plot_tip,
    lineage = `ordinal_species$\`Extended lineage\``,
    completed = `ordinal_species$complete_status`
  ) %>%
  dplyr::mutate(
    order = ordinal_list$`Order (NCBI)`[
      match(label, ordinal_list$`Scientific Name`)
    ],
    assembly_size = ordinal_list$`Assembly Size`[
      match(label, ordinal_list$`Scientific Name`)
    ],
    total_species_order = ordinal_list$`# species/order`[
      match(label, ordinal_list$`Scientific Name`)
    ],
    sequenced_species_order = ordinal_species$n_sequenced[
      match(label, ordinal_species$`Scientific Name`)
    ],
    label_to_plot = ifelse(completed, label, NA_character_)
  )

row.names(metadata) <- metadata$label

lineage_colors <- metadata %>%
  dplyr::select(lineage) %>%
  dplyr::distinct()

metadata$lineage <- factor(metadata$lineage, levels = lineage_colors$lineage)

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

ggsave(
  filename = "ordinal_tree.svg",
  plot = p_ordinal2,
  width = 15,
  height = 12,
  device = "svg"
)

###### Full vertebrate tree ######################################################################################################################################################

# ---------------------------------------------------------------------------
# Read updated VGP table
# ---------------------------------------------------------------------------

vgp_tbl <- googlesheets4::read_sheet(
  ss = "1Jwjv6Kwc6VIn1UMMhnG6kvFCxjwGdC5b7p_HtbDOMOs",
  sheet = 1
) %>%
  tibble::as_tibble()

required_columns <- c(
  "Scientific Name",
  "Extended Lineage",
  "Accession # for main haplotype",
  "# species/order"
)

missing_columns <- setdiff(required_columns, names(vgp_tbl))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns in the updated VGP table: ",
    paste(missing_columns, collapse = ", ")
  )
}

# ---------------------------------------------------------------------------
# Extended Lineages represented in the vertebrate figure
# ---------------------------------------------------------------------------

lineage_components <- tibble::tribble(
  ~lineage,                ~ott_query,             ~fallback_tip,
  "Mammals",               "Mammalia",             NA_character_,
  "Birds",                 "Aves",                 NA_character_,
  "Crocodilians",          "Alligatoridae",        NA_character_,
  "Crocodilians",          "Longirostres",         NA_character_,
  "Turtles",               "Testudines",           NA_character_,
  "Lepidosauria",          "Lepidosauria",         NA_character_,
  "Amphibians",            "Amphibia",             NA_character_,
  "Lobe-finned fishes",    "Dipnoi",               "Protopterus annectens",
  "Lobe-finned fishes",    "Coelacanthimorpha",    "Latimeria chalumnae",
  "Ray-finned fishes",     "Actinopterygii",       NA_character_,
  "Cartilaginous fishes",  "Chondrichthyes",       NA_character_,
  "Cyclostomes",           "Cyclostomata",         NA_character_
)

extended_lineage_order <- unique(lineage_components$lineage)

class_colors <- c(
  "Mammals"              = "#E69F00",
  "Birds"                = "#00796B",
  "Crocodilians"         = "#009688",
  "Turtles"              = "#4DB6AC",
  "Lepidosauria"         = "#80CBC4",
  "Amphibians"           = "#984EA3",
  "Lobe-finned fishes"   = "#A6761D",
  "Ray-finned fishes"    = "#56B4E9",
  "Cartilaginous fishes" = "#0072B2",
  "Cyclostomes"          = "#CC79A7"
)

vgp_green <- "#00A651"

# ---------------------------------------------------------------------------
# Select VGP vertebrate species with a genome in the data freeze
# ---------------------------------------------------------------------------

has_accession <- function(x) {
  !is.na(x) & trimws(as.character(x)) != ""
}

represented_metadata <- vgp_tbl %>%
  dplyr::filter(
    !is.na(`Scientific Name`),
    !is.na(`Extended Lineage`),
    has_accession(`Accession # for main haplotype`),
    `Extended Lineage` %in% extended_lineage_order
  ) %>%
  dplyr::distinct(`Scientific Name`, `Extended Lineage`, .keep_all = TRUE)

excluded_outgroups <- vgp_tbl %>%
  dplyr::filter(
    !is.na(`Scientific Name`),
    has_accession(`Accession # for main haplotype`),
    is.na(`Extended Lineage`) |
      !`Extended Lineage` %in% extended_lineage_order
  ) %>%
  dplyr::distinct(`Scientific Name`, `Extended Lineage`, .keep_all = TRUE) %>%
  dplyr::select(
    `Scientific Name`,
    `Extended Lineage`,
    `Accession # for main haplotype`
  )

cat(
  "\nVGP input for vertebrate figure\n",
  "-------------------------------\n",
  "All accession-bearing rows:                    ",
  sum(has_accession(vgp_tbl$`Accession # for main haplotype`)), "\n",
  "Vertebrate species retained for the figure:    ",
  nrow(represented_metadata), "\n",
  "Outgroup/non-vertebrate species excluded:      ",
  nrow(excluded_outgroups), "\n",
  sep = ""
)

if (nrow(excluded_outgroups) > 0) {
  cat("\nExcluded outgroup/non-vertebrate species:\n")
  print(tibble::as_tibble(excluded_outgroups), n = Inf)
}

write.table(
  excluded_outgroups,
  file = "full_tree_excluded_outgroups.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ---------------------------------------------------------------------------
# Resolve Extended Lineage clades in OpenTree
# ---------------------------------------------------------------------------

lineage_resolved <- rotl::tnrs_match_names(
  names = lineage_components$ott_query
)

lineage_components <- dplyr::bind_cols(
  lineage_components,
  lineage_resolved %>%
    dplyr::select(
      unique_name,
      ott_id,
      approximate_match,
      score,
      is_synonym
    )
)

if (any(is.na(lineage_components$ott_id))) {
  cat("\nExtended Lineage components not resolved by OpenTree:\n")
  print(
    tibble::as_tibble(
      lineage_components %>%
        dplyr::filter(is.na(ott_id))
    ),
    n = Inf
  )
  stop("One or more Extended Lineage components could not be resolved.")
}

# ---------------------------------------------------------------------------
# Retrieve OpenTree species descendants for each Extended Lineage component
# ---------------------------------------------------------------------------

extract_descendant_table <- function(x) {
  
  if (inherits(x, "data.frame")) {
    return(tibble::as_tibble(x, rownames = "taxon_name"))
  }
  
  if (is.list(x)) {
    
    is_df <- vapply(
      x,
      function(z) inherits(z, "data.frame"),
      logical(1)
    )
    
    if (sum(is_df) == 0) {
      stop("No data.frame was returned by datelife::get_ott_children().")
    }
    
    return(
      tibble::as_tibble(
        x[[which(is_df)[1]]],
        rownames = "taxon_name"
      )
    )
  }
  
  stop("Unexpected output from datelife::get_ott_children().")
}

get_component_species <- function(lineage, ott_query, ott_id, fallback_tip) {
  
  descendant_result <- datelife::get_ott_children(
    ott_ids = setNames(ott_id, ott_query),
    ott_rank = "species"
  )
  
  descendants <- extract_descendant_table(descendant_result)
  
  if ("rank" %in% names(descendants)) {
    descendants <- descendants %>%
      dplyr::filter(rank == "species")
  }
  
  descendant_ids <- unique(stats::na.omit(descendants$ott_id))
  
  if (length(descendant_ids) == 0 && !is.na(fallback_tip)) {
    
    fallback_resolved <- rotl::tnrs_match_names(
      names = fallback_tip
    )
    
    descendant_ids <- unique(
      stats::na.omit(fallback_resolved$ott_id)
    )
  }
  
  if (length(descendant_ids) == 0) {
    stop("No species descendants retrieved for lineage component: ", ott_query)
  }
  
  tibble::tibble(
    lineage = lineage,
    component = ott_query,
    ott_id = descendant_ids
  )
}

lineage_species_cache <- "full_tree_extended_lineage_species_ott_ids.rds"

if (!file.exists(lineage_species_cache)) {
  
  component_species_ott <- purrr::pmap_dfr(
    lineage_components %>%
      dplyr::select(lineage, ott_query, ott_id, fallback_tip),
    get_component_species
  ) %>%
    dplyr::distinct(lineage, ott_id, .keep_all = TRUE)
  
  saveRDS(
    component_species_ott,
    file = lineage_species_cache
  )
  
} else {
  
  component_species_ott <- readRDS(lineage_species_cache)
}

# The plotted universe is the union of the displayed vertebrate lineages.
# Non-vertebrate outgroups are excluded from the plotted tree.
vertebrate_species_ott_ids <- unique(component_species_ott$ott_id)

cat(
  "\nOpenTree species retrieved across displayed vertebrate lineages: ",
  length(vertebrate_species_ott_ids), "\n",
  sep = ""
)

# ---------------------------------------------------------------------------
# Build the vertebrate-only OpenTree species tree
# ---------------------------------------------------------------------------

tree_cache <- "full_tree_vertebrate_extended_lineages_species_only.rds"

if (!file.exists(tree_cache)) {
  
  vertebrate_tree_species_only <- rotl::tol_induced_subtree(
    ott_ids = vertebrate_species_ott_ids
  )
  
  saveRDS(
    vertebrate_tree_species_only,
    file = tree_cache
  )
  
} else {
  
  vertebrate_tree_species_only <- readRDS(tree_cache)
}

make_tip_map <- function(tree) {
  
  tip_map <- tibble::tibble(
    raw_tip_label = tree$tip.label,
    ott_id = as.numeric(
      stringr::str_extract(raw_tip_label, "(?<=ott)[0-9]+$")
    ),
    clean_tip_label = gsub(
      "_",
      " ",
      rotl::strip_ott_ids(raw_tip_label),
      fixed = TRUE
    )
  )
  
  if (any(is.na(tip_map$ott_id))) {
    cat("\nTree tips without recoverable OTT identifiers:\n")
    print(
      tibble::as_tibble(
        tip_map %>%
          dplyr::filter(is.na(ott_id))
      ),
      n = Inf
    )
    stop("Could not recover OTT identifiers from one or more tree tips.")
  }
  
  tip_map
}

vertebrate_tip_map <- make_tip_map(vertebrate_tree_species_only)

cat(
  "\nPlotted OpenTree vertebrate tree\n",
  "-------------------------------\n",
  "Terminal species tips: ",
  length(vertebrate_tree_species_only$tip.label), "\n",
  sep = ""
)

# ---------------------------------------------------------------------------
# Build Extended Lineage gear metadata
# ---------------------------------------------------------------------------

lineage_tip_metadata <- component_species_ott %>%
  dplyr::inner_join(
    vertebrate_tip_map %>%
      dplyr::select(ott_id, label = raw_tip_label),
    by = "ott_id"
  ) %>%
  dplyr::select(lineage, label) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    lineage = factor(lineage, levels = extended_lineage_order)
  )

# ---------------------------------------------------------------------------
# Extended Lineage coverage table for manual placement beside the legend
# ---------------------------------------------------------------------------
# Numerator: accession-bearing VGP species from the updated full table.
# Denominator: sum of '# species/order' within each Extended Lineage.
# OpenTree tips are retained only as a diagnostic for the rendered tree.
# ---------------------------------------------------------------------------

parse_species_count <- function(x) {
  suppressWarnings(
    as.numeric(
      gsub("[^0-9.]", "", as.character(x))
    )
  )
}

lineage_vgp_counts <- represented_metadata %>%
  dplyr::count(`Extended Lineage`, name = "n_vgp_species") %>%
  dplyr::rename(lineage = `Extended Lineage`)

lineage_estimated_species_counts <- vgp_tbl %>%
  dplyr::filter(
    !is.na(`Extended Lineage`),
    `Extended Lineage` %in% extended_lineage_order,
    !is.na(`# species/order`)
  ) %>%
  dplyr::mutate(
    species_count = parse_species_count(`# species/order`)
  ) %>%
  dplyr::group_by(`Extended Lineage`) %>%
  dplyr::summarise(
    estimated_total_species = sum(species_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::rename(lineage = `Extended Lineage`)

lineage_opentree_counts <- lineage_tip_metadata %>%
  dplyr::mutate(
    lineage = as.character(lineage)
  ) %>%
  dplyr::count(lineage, name = "n_opentree_species")

lineage_coverage <- tibble::tibble(
  lineage = extended_lineage_order
) %>%
  dplyr::left_join(lineage_vgp_counts, by = "lineage") %>%
  dplyr::left_join(lineage_estimated_species_counts, by = "lineage") %>%
  dplyr::left_join(lineage_opentree_counts, by = "lineage") %>%
  dplyr::mutate(
    n_vgp_species = tidyr::replace_na(n_vgp_species, 0L),
    coverage_estimated_pct = round(
      100 * n_vgp_species / estimated_total_species,
      2
    ),
    coverage_opentree_pct = round(
      100 * n_vgp_species / n_opentree_species,
      2
    )
  )

lineage_coverage_for_figure <- lineage_coverage %>%
  dplyr::transmute(
    `Extended Lineage` = lineage,
    `VGP species with genome` = n_vgp_species,
    `Estimated extant species` = estimated_total_species,
    `Coverage (%)` = coverage_estimated_pct
  )

lineage_coverage_with_opentree_diagnostic <- lineage_coverage %>%
  dplyr::transmute(
    `Extended Lineage` = lineage,
    `VGP species with genome` = n_vgp_species,
    `Estimated extant species` = estimated_total_species,
    `Coverage of estimated diversity (%)` = coverage_estimated_pct,
    `OpenTree tips in plotted tree` = n_opentree_species,
    `Coverage of OpenTree tips (%)` = coverage_opentree_pct
  )

cat("\nExtended Lineage coverage table for figure annotation:\n")
print(tibble::as_tibble(lineage_coverage_for_figure), n = Inf)

cat(
  "\nTotals across displayed vertebrate Extended Lineages\n",
  "---------------------------------------------------\n",
  "VGP species with genome:  ",
  sum(lineage_coverage_for_figure$`VGP species with genome`, na.rm = TRUE), "\n",
  "Estimated extant species: ",
  format(
    sum(lineage_coverage_for_figure$`Estimated extant species`, na.rm = TRUE),
    big.mark = ","
  ), "\n",
  "OpenTree tips in tree:    ",
  format(
    sum(
      lineage_coverage_with_opentree_diagnostic$`OpenTree tips in plotted tree`,
      na.rm = TRUE
    ),
    big.mark = ","
  ), "\n",
  sep = ""
)

write.table(
  lineage_coverage_for_figure,
  file = "full_tree_extended_lineage_coverage_for_figure.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  lineage_coverage_with_opentree_diagnostic,
  file = "full_tree_extended_lineage_coverage_with_opentree_diagnostic.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ---------------------------------------------------------------------------
# Resolve VGP species against OpenTree
# ---------------------------------------------------------------------------

species_reconciliation <- represented_metadata %>%
  dplyr::transmute(
    spreadsheet_name = `Scientific Name`,
    lineage = `Extended Lineage`,
    opentree_query = `Scientific Name`
  )

# Use only for strict naming normalization, not proxy substitution.
name_query_corrections <- c(
  "Amia calva ocellicauda" = "Amia calva"
)

correction_index <- match(
  species_reconciliation$opentree_query,
  names(name_query_corrections)
)

species_reconciliation$opentree_query[!is.na(correction_index)] <-
  unname(name_query_corrections[correction_index[!is.na(correction_index)]])

species_reconciliation <- species_reconciliation %>%
  dplyr::mutate(
    query_corrected = spreadsheet_name != opentree_query
  )

vgp_resolved <- rotl::tnrs_match_names(
  names = species_reconciliation$opentree_query
)

species_reconciliation <- dplyr::bind_cols(
  species_reconciliation,
  vgp_resolved %>%
    dplyr::select(
      unique_name,
      ott_id,
      approximate_match,
      score,
      is_synonym
    )
) %>%
  dplyr::mutate(
    resolved_by_opentree = !is.na(ott_id),
    actual_tip_present = resolved_by_opentree &
      ott_id %in% vertebrate_tip_map$ott_id,
    actual_tip_label = vertebrate_tip_map$raw_tip_label[
      match(ott_id, vertebrate_tip_map$ott_id)
    ],
    in_synthetic_tree = FALSE
  )

resolved_rows <- which(species_reconciliation$resolved_by_opentree)

if (length(resolved_rows) > 0) {
  species_reconciliation$in_synthetic_tree[resolved_rows] <-
    rotl::is_in_tree(species_reconciliation$ott_id[resolved_rows])
}

# ---------------------------------------------------------------------------
# Automatically find plotted-tip proxies for resolved but unplotted taxa
# ---------------------------------------------------------------------------

proxy_candidates <- species_reconciliation %>%
  dplyr::filter(
    resolved_by_opentree,
    in_synthetic_tree,
    !actual_tip_present
  )

cat(
  "\nInitial VGP recovery against species-only tree\n",
  "----------------------------------------------\n",
  "VGP vertebrate species with genome:           ",
  nrow(species_reconciliation), "\n",
  "Resolved by OpenTree:                         ",
  sum(species_reconciliation$resolved_by_opentree), "\n",
  "Present as plotted tips:                      ",
  sum(species_reconciliation$actual_tip_present), "\n",
  "Eligible for automatic proxy search:          ",
  nrow(proxy_candidates), "\n",
  "Unresolved or absent from synthetic tree:     ",
  sum(
    !species_reconciliation$resolved_by_opentree |
      !species_reconciliation$in_synthetic_tree
  ), "\n",
  sep = ""
)

descendant_tip_nodes <- function(tree, node) {
  
  children <- tree$edge[tree$edge[, 1] == node, 2]
  
  tip_children <- children[
    children <= ape::Ntip(tree)
  ]
  
  internal_children <- children[
    children > ape::Ntip(tree)
  ]
  
  if (length(internal_children) > 0) {
    tip_children <- c(
      tip_children,
      unlist(
        lapply(
          internal_children,
          function(x) descendant_tip_nodes(tree, x)
        )
      )
    )
  }
  
  unique(tip_children)
}

auto_proxy_results <- tibble::tibble(
  ott_id = numeric(),
  proxy_tip_label = character(),
  proxy_tip_clean_name = character(),
  proxy_ott_id = numeric(),
  proxy_candidate_count = integer(),
  proxy_ancestor_steps = integer()
)

if (nrow(proxy_candidates) > 0) {
  
  augmented_ids <- unique(
    c(
      vertebrate_tip_map$ott_id,
      proxy_candidates$ott_id
    )
  )
  
  signature_file <- tempfile()
  writeLines(as.character(sort(augmented_ids)), signature_file)
  augmented_signature <- unname(tools::md5sum(signature_file))
  unlink(signature_file)
  
  augmented_cache <- paste0(
    "full_tree_augmented_for_proxy_assignment_",
    substr(augmented_signature, 1, 12),
    ".rds"
  )
  
  if (!file.exists(augmented_cache)) {
    
    augmented_tree <- rotl::tol_induced_subtree(
      ott_ids = augmented_ids
    )
    
    saveRDS(
      augmented_tree,
      file = augmented_cache
    )
    
  } else {
    
    augmented_tree <- readRDS(augmented_cache)
  }
  
  augmented_tip_map <- make_tip_map(augmented_tree) %>%
    dplyr::mutate(
      tree_node = dplyr::row_number()
    )
  
  base_plotted_ids <- vertebrate_tip_map$ott_id
  
  find_proxy_tip <- function(target_ott_id) {
    
    target_node <- augmented_tip_map$tree_node[
      match(target_ott_id, augmented_tip_map$ott_id)
    ]
    
    if (is.na(target_node)) {
      return(
        tibble::tibble(
          ott_id = target_ott_id,
          proxy_tip_label = NA_character_,
          proxy_tip_clean_name = NA_character_,
          proxy_ott_id = NA_real_,
          proxy_candidate_count = NA_integer_,
          proxy_ancestor_steps = NA_integer_
        )
      )
    }
    
    current_node <- target_node
    ancestor_steps <- 0L
    
    repeat {
      
      parent_row <- match(current_node, augmented_tree$edge[, 2])
      
      if (is.na(parent_row)) {
        return(
          tibble::tibble(
            ott_id = target_ott_id,
            proxy_tip_label = NA_character_,
            proxy_tip_clean_name = NA_character_,
            proxy_ott_id = NA_real_,
            proxy_candidate_count = NA_integer_,
            proxy_ancestor_steps = ancestor_steps
          )
        )
      }
      
      parent_node <- augmented_tree$edge[parent_row, 1]
      ancestor_steps <- ancestor_steps + 1L
      
      descendant_nodes <- descendant_tip_nodes(
        augmented_tree,
        parent_node
      )
      
      candidates <- augmented_tip_map %>%
        dplyr::filter(
          tree_node %in% descendant_nodes,
          ott_id %in% base_plotted_ids,
          ott_id != target_ott_id
        ) %>%
        dplyr::arrange(clean_tip_label, ott_id)
      
      if (nrow(candidates) > 0) {
        
        chosen <- candidates[1, ]
        
        return(
          tibble::tibble(
            ott_id = target_ott_id,
            proxy_tip_label = chosen$raw_tip_label,
            proxy_tip_clean_name = chosen$clean_tip_label,
            proxy_ott_id = chosen$ott_id,
            proxy_candidate_count = nrow(candidates),
            proxy_ancestor_steps = ancestor_steps
          )
        )
      }
      
      current_node <- parent_node
    }
  }
  
  auto_proxy_results <- dplyr::bind_rows(
    lapply(
      unique(proxy_candidates$ott_id),
      find_proxy_tip
    )
  )
}

species_reconciliation <- species_reconciliation %>%
  dplyr::left_join(auto_proxy_results, by = "ott_id") %>%
  dplyr::mutate(
    plotted_tip_to_highlight = dplyr::case_when(
      actual_tip_present ~ actual_tip_label,
      !is.na(proxy_tip_label) ~ proxy_tip_label,
      TRUE ~ NA_character_
    ),
    highlight_method = dplyr::case_when(
      actual_tip_present ~ "Actual plotted tip",
      !is.na(proxy_tip_label) ~ "Automatic nearest plotted proxy",
      TRUE ~ "Not highlighted"
    )
  )

# ---------------------------------------------------------------------------
# Proxy audit and proxy-use statistics
# ---------------------------------------------------------------------------

proxy_audit <- species_reconciliation %>%
  dplyr::select(
    spreadsheet_name,
    lineage,
    opentree_query,
    query_corrected,
    unique_name,
    ott_id,
    resolved_by_opentree,
    in_synthetic_tree,
    actual_tip_present,
    actual_tip_label,
    proxy_tip_label,
    proxy_tip_clean_name,
    proxy_ott_id,
    proxy_candidate_count,
    proxy_ancestor_steps,
    plotted_tip_to_highlight,
    highlight_method
  ) %>%
  dplyr::arrange(highlight_method, lineage, spreadsheet_name)

proxy_used <- proxy_audit %>%
  dplyr::filter(
    highlight_method == "Automatic nearest plotted proxy"
  )

actual_used <- proxy_audit %>%
  dplyr::filter(
    highlight_method == "Actual plotted tip"
  )

not_highlighted <- proxy_audit %>%
  dplyr::filter(
    highlight_method == "Not highlighted"
  )

proxy_summary <- tibble::tibble(
  metric = c(
    "Total VGP vertebrate species with genome",
    "Represented by actual OpenTree tip",
    "Represented by automatic plotted proxy",
    "Not represented after proxy search",
    "Unique proxy tips used",
    "Proxy tips representing >1 VGP species"
  ),
  value = c(
    nrow(proxy_audit),
    nrow(actual_used),
    nrow(proxy_used),
    nrow(not_highlighted),
    dplyr::n_distinct(proxy_used$proxy_tip_label),
    proxy_used %>%
      dplyr::count(proxy_tip_label) %>%
      dplyr::filter(n > 1) %>%
      nrow()
  )
) %>%
  dplyr::mutate(
    percent_of_vgp_species = round(
      100 * value / nrow(proxy_audit),
      2
    )
  )

proxy_by_lineage <- proxy_audit %>%
  dplyr::group_by(lineage) %>%
  dplyr::summarise(
    n_vgp_species = dplyr::n(),
    n_actual_tips = sum(highlight_method == "Actual plotted tip"),
    n_proxy_tips = sum(highlight_method == "Automatic nearest plotted proxy"),
    n_not_highlighted = sum(highlight_method == "Not highlighted"),
    percent_using_proxy = round(
      100 * n_proxy_tips / n_vgp_species,
      2
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    factor(lineage, levels = extended_lineage_order)
  )

proxy_substitution_table <- proxy_used %>%
  dplyr::transmute(
    `VGP species` = spreadsheet_name,
    `Extended Lineage` = lineage,
    `Resolved OpenTree taxon` = unique_name,
    `Highlighted proxy tip` = proxy_tip_clean_name,
    `Proxy candidate tips in nearest clade` = proxy_candidate_count,
    `Ancestor steps to proxy clade` = proxy_ancestor_steps
  ) %>%
  dplyr::arrange(
    factor(`Extended Lineage`, levels = extended_lineage_order),
    `VGP species`
  )

cat("\nProxy usage summary:\n")
print(tibble::as_tibble(proxy_summary), n = Inf)

cat("\nProxy use by Extended Lineage:\n")
print(tibble::as_tibble(proxy_by_lineage), n = Inf)

cat("\nAutomatic proxy substitutions:\n")
print(tibble::as_tibble(proxy_substitution_table), n = Inf)

if (nrow(not_highlighted) > 0) {
  cat("\nSpecies remaining without a plotted representation:\n")
  print(tibble::as_tibble(not_highlighted), n = Inf)
}

write.table(
  proxy_audit,
  file = "full_tree_species_highlight_and_proxy_audit.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  proxy_summary,
  file = "full_tree_proxy_usage_summary.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  proxy_by_lineage,
  file = "full_tree_proxy_usage_by_extended_lineage.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  proxy_substitution_table,
  file = "full_tree_proxy_substitutions.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  not_highlighted,
  file = "full_tree_species_not_highlighted_after_proxy.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ---------------------------------------------------------------------------
# Select displayed tips and report tip collisions
# ---------------------------------------------------------------------------

represented_tips <- unique(
  stats::na.omit(species_reconciliation$plotted_tip_to_highlight)
)

tip_collisions <- species_reconciliation %>%
  dplyr::filter(!is.na(plotted_tip_to_highlight)) %>%
  dplyr::count(plotted_tip_to_highlight, name = "n_vgp_species") %>%
  dplyr::filter(n_vgp_species > 1)

if (nrow(tip_collisions) > 0) {
  cat("\nMultiple VGP species represented by the same displayed tip:\n")
  print(tibble::as_tibble(tip_collisions), n = Inf)
}

cat(
  "\nFinal highlighted-tip summary\n",
  "-----------------------------\n",
  "VGP vertebrate species with genome: ",
  nrow(species_reconciliation), "\n",
  "VGP species assigned a display tip: ",
  sum(!is.na(species_reconciliation$plotted_tip_to_highlight)), "\n",
  "Unique tree tips highlighted:       ",
  length(represented_tips), "\n",
  "VGP species not represented:        ",
  nrow(not_highlighted), "\n",
  sep = ""
)

if (length(represented_tips) == 0) {
  stop("No VGP tips were assigned for plotting.")
}

represented_grp <- list("VGP Phase I" = represented_tips)

# ---------------------------------------------------------------------------
# Plot: black tree, green represented/proxy branches, Extended Lineage gear
# ---------------------------------------------------------------------------

tree_annotated <- ggtree::groupOTU(
  vertebrate_tree_species_only,
  represented_grp,
  group_name = "status"
)

p <- ggtree::ggtree(
  tree_annotated,
  layout = "fan",
  color = "black",
  linewidth = 0.10
)

tree_overlay_data <- p$data %>%
  dplyr::mutate(
    is_vgp = !is.na(status) & as.character(status) != "0",
    vgp_overlay = ifelse(is_vgp, vgp_green, NA_character_)
  )

green_tips <- tree_overlay_data %>%
  dplyr::filter(isTip, is_vgp) %>%
  dplyr::distinct(label)

cat(
  "\nGreen-overlay diagnostic\n",
  "------------------------\n",
  "Assigned display tips:               ",
  length(represented_tips), "\n",
  "Terminal tips highlighted green:     ",
  nrow(green_tips), "\n",
  "Green branches/nodes with internals: ",
  sum(tree_overlay_data$is_vgp), "\n",
  sep = ""
)

# Green VGP / proxy branches above the black tree
p1 <- p +
  ggtree::geom_tree(
    data = tree_overlay_data,
    mapping = aes(color = vgp_overlay),
    linewidth = 0.22,
    na.rm = TRUE
  ) +
  scale_color_identity(
    name = "Species represented",
    breaks = vgp_green,
    labels = "VGP Phase I or nearest plotted proxy",
    na.value = "transparent",
    guide = guide_legend(order = 2)
  )

# Extended Lineage gear added last so it remains above branch layers
p2 <- p1 +
  ggnewscale::new_scale_fill() +
  geom_fruit(
    data = lineage_tip_metadata,
    geom = geom_tile,
    mapping = aes(y = label, fill = lineage),
    width = 0.85,
    offset = 0
  ) +
  scale_fill_manual(
    name = "Extended Lineage",
    values = class_colors[extended_lineage_order],
    breaks = extended_lineage_order,
    guide = guide_legend(
      keywidth = 0.35,
      keyheight = 0.35,
      ncol = 2,
      order = 1
    )
  ) +
  theme(
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 20),
    legend.box.spacing = unit(3, "mm"),
    plot.margin = margin(40, 80, 40, 20)
  )

ggsave(
  filename = "full_tree_extended_lineage_gear_with_proxies.svg",
  plot = p2,
  width = 14,
  height = 12,
  device = "svg"
)