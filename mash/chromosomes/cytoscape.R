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
setNodeColorDefault('#D8D8D8')

network <- getNetworkSuid()
#getVisualPropertyNames() # handy to find node properties

commandsRun('table list type=node')
getNetworkList()
setCurrentNetwork('louvain_(none)_outliers.csv')
communities <- getTableColumns(columns = 'CD_MemberList') 
communities <- sapply(strsplit(as.character(communities$CD_MemberList), " "), function(x) x)

communities <- Filter(function(x){length(x)<1000},communities) # remove very large communities

len = lengths(communities)
community_assignments <- data.frame(name = unlist(communities), community_leuvain = rep(seq_along(len), len), id = sequence(len))

setCurrentNetwork('outliers.csv')

loadTableData(
  community_assignments,
  data.key.column = "name",
  table.key.column = "name",
)

# clock coordinates
community_list <- unique(community_assignments$community_leuvain)
community_count <- length(community_list)
degrees_per_iter <- 360 / community_count
start_angle_deg <- 90

for (i in 1:community_count) {
  
  filter <- createColumnFilter(filter.name=as.character(community_list[[i]]), column='community_leuvain', community_list[[i]], 'IS')
  #getSelectedNodes() # handy handle to selected nodes
  positions <- getNodePosition(filter$nodes)
  
  angle_deg <- (start_angle_deg + ((i-1) * degrees_per_iter)) %% 360
  
  x <- 100 * cos(angle_deg * pi/180)
  y <- 100 * sin(angle_deg * pi/180)
  
  positions %>% by_row(place, mean(positions$x_location), mean(positions$y_location), x, y)
}
