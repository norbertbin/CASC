# ---------------------------------------------------------------------
# load reqired libraries
# ---------------------------------------------------------------------
if (!require(gplots)) {
    install.packages('gplots', dependencies = T)
    require(gplots)
}

# ---------------------------------------------------------------------
# some cluster ploting helper functions
# ---------------------------------------------------------------------
# get the layer assignments for nodes along a given coordinate
getLayer = function(nLayers, covVec, nNodes) {
    cRange = range(covVec)
    cutoffs = seq(cRange[1], cRange[2], length.out = nLayers + 1)
    layer = rep(0, nNodes)
    #get layer assignments for each node
    for(i in 1:nNodes) {
        layer[i] = sum(covVec[i] >= cutoffs[1:nLayers])
    }

    return(layer)
}

# plot the nodes in a given layer
plotLayer = function(cluster, layerData, orthCoord, colSeq,
    outFileName) {

    # determine the correct plot formula
    if(orthCoord == 'z'){
        plotFormula = y~x
    } else if (orthCoord == 'x') {
        plotFormula = z~y
    } else {
        plotFormula = z~x
    }

    png(outFileName)
    print(
        xyplot(plotFormula, group = cluster, data = layerData,
               pch = 20, col = colSeq[sort(unique(cluster))],
               main = paste(orthCoord, ' = ', round(min(layerData$z), 3),
                   ' to ', round(max(layerData$z), 3)))                         
        )
    dev.off()
}

# ---------------------------------------------------------------------
# functions for estimating and plotting block matrix
# ---------------------------------------------------------------------
# estimates block matrix based on given adjacency matrix and block memberships 
# includeBlocks specifies which block should be included in the block matrix
estBlockMat = function(adjMat, nodeBlocks) {
  
  nodeCounts = tabulate(nodeBlocks)
  nodeCounts = nodeCounts[nodeCounts > 0]
  uniqueBlocks = sort(unique(nodeBlocks))
  nBlocks = length(uniqueBlocks)
  bMat = matrix(0, nrow = nBlocks, ncol = nBlocks)
  
  for(i in 1:nBlocks) {
    for(j in i:nBlocks) {
      if(i == j){
        bMat[i, i] = nnzero(adjMat[nodeBlocks == uniqueBlocks[i],
                nodeBlocks == uniqueBlocks[i]]) /
          (nodeCounts[i]^2 - nodeCounts[i])
      }
      else {
        bMat[i, j] = nnzero(adjMat[nodeBlocks == uniqueBlocks[i],
                nodeBlocks == uniqueBlocks[j]]) /
          (nodeCounts[i]*nodeCounts[j])
        bMat[j, i] = bMat[i, j]
      }
    }
  }

  return(bMat)
} 

#plot heat map of estimated block matrix
plotB = function(bMat, outFileName) {

    pdf(outFileName, width = 10, height = 10)
    heatmap.2(bMat, Colv = "Rowv", dendrogram = "none", key = F,
              trace = "none",
              col = colorRampPalette(c("white", "black"))(256),
              breaks = getVectorBins(as.vector(bMat), 257),
              lwid = c(.09, .99), lhei = c(.01, .99),
              margins = c(5, 5))
    dev.off()
    
    return(0)
}

# generate bin break points so that each bin has equal frequency
getVectorBins = function(vec, nBins) {
    vec = unique(vec)
    nVec = length(vec)
    return(sort(vec)[round(nVec/nBins*(1:nBins))])
}
