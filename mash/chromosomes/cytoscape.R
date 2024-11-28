############# SETUP ##############

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

requiredpackages <- c("purrrlyr", 'dplyr')

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

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RCy3")
library(RCy3)

############# FUNCTIONS ##############
#move nodes
place <- function(row, centroid_x, centroid_y, x, y){
  
  setNodePropertyBypass(
    row.names(row),
    row$x_location-centroid_x+x,
    "NODE_X_LOCATION",
    bypass = TRUE
  )
  setNodePropertyBypass(
    row.names(row),
    row$y_location-centroid_y+y,
    "NODE_Y_LOCATION",
    bypass = TRUE
  )
}

##########################

openSession("chromosomes.cys")

setVisualStyle('default')
setEdgeLineWidthDefault(0.5)
setEdgeColorDefault('#f0f0f0')

network <- getNetworkSuid()
#getVisualPropertyNames() # handy to find node properties

commandsRun('table list type=node')
getNetworkList()
setCurrentNetwork('louvain_(none)_outliers.csv')
communities <- getTableColumns(columns = 'CD_MemberList') 
communities <- sapply(strsplit(as.character(communities$CD_MemberList), " "), function(x) x)

communities <- Filter(function(x){length(x)<1000},communities) # remove overlapping large community


len = lengths(communities)
community_assignments <- data.frame(name = unlist(communities), community_leuvain = rep(seq_along(len), len), id = sequence(len))

setCurrentNetwork('outliers.csv')

loadTableData(
  community_assignments,
  data.key.column = "name",
  table.key.column = "name",
)

# clock coordinates
large_communities <- community_assignments %>% group_by(community_leuvain) %>% count(community_leuvain, .drop = FALSE) %>% filter(n > 10) %>% sort_by(~n, decreasing=TRUE) # large communities
community_list <- large_communities$community_leuvain
community_count <- length(community_list)
degrees_per_iter <- 360 / community_count
start_angle_deg <- 270

for (i in 1:community_count) {
  
  createGroupByColumn(as.character(community_list[[i]]), column = 'community_leuvain', value = community_list[[i]])
  group <- collapseGroup(groups = as.character(community_list[[i]]))

  angle_deg <- (start_angle_deg + ((i-1) * degrees_per_iter)) %% 360
  x <- 5000 * cos(angle_deg * pi/180)
  y <- 5000 * sin(angle_deg * pi/180)
  
  setNodePropertyBypass(community_list[[i]], x, "NODE_X_LOCATION", bypass = TRUE)
  setNodePropertyBypass(community_list[[i]], y, "NODE_Y_LOCATION", bypass = TRUE)
  
  expandGroup(groups = community_list[[i]])
}

#small communities
small_communities <- community_assignments %>% group_by(community_leuvain) %>% count(community_leuvain, .drop = FALSE) %>% filter(n <= 10) %>% sort_by(~n, decreasing=TRUE) # small communities

community_list <- small_communities$community_leuvain
community_count <- length(community_list)
degrees_per_iter <- 360 / community_count
start_angle_deg <- 270
e <- 0
y <- 6000

for (i in 1:community_count) {
  
  createGroupByColumn(as.character(community_list[[i]]), column = 'community_leuvain', value = community_list[[i]])
  group <- collapseGroup(groups = as.character(community_list[[i]]))
  print(group)
  
  x <- -5000 + e * 200

  if (i > community_count/2) {
    e <- 0
    y = 6500
  }
  e <- e + 1
  
  setNodePropertyBypass(community_list[[i]], x, "NODE_X_LOCATION", bypass = TRUE)
  setNodePropertyBypass(community_list[[i]], y, "NODE_Y_LOCATION", bypass = TRUE)
  
  expandGroup(groups = community_list[[i]])
}
