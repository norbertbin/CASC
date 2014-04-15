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
