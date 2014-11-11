# ---------------------------------------------------------------------
# R script takes name of a graphml file as input, opens that file extracts
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# load igraph
# --------------------------------------------------------------------
library(igraph)

# ---------------------------------------------------------------------
# get command line arguement with graph file name
# ---------------------------------------------------------------------
comArgs = commandArgs(T)
graphFile = comArgs[1]
outputFile = comArgs[2]

# ---------------------------------------------------------------------
# a function to extract the largest connected component in the graph
# ---------------------------------------------------------------------
giantComponent <- function(graph, ...) {
      cl <- clusters(graph, ...)
        induced.subgraph(graph, which(cl$membership == which.max(cl$csize)))
  }


# ---------------------------------------------------------------------
# read in graph file and extract attributes 
# ---------------------------------------------------------------------
brainGraph = read.graph(graphFile, format='graphml')
brainGraph = giantComponent(brainGraph)
atlasLabels = V(brainGraph)$atlas_1_region_num
write.table(atlasLabels, file=outputFile)
