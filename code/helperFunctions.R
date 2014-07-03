# ---------------------------------------------------------------------
# load reqired libraries
# ---------------------------------------------------------------------
if (!require(gplots)) {
    install.packages('gplots', dependencies = T)
    require(gplots)
}
if (!require(clue)) {
    install.packages('clue', dependencies = T)
    require(clue)
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

estBlockMatSparse = function(adjMat, nodeBlocks) {
  
   nodeCounts = tabulate(nodeBlocks)
   nodeCounts = nodeCounts[nodeCounts > 0]
   uniqueBlocks = sort(unique(nodeBlocks))
   nBlocks = length(uniqueBlocks)
   bMat = matrix(0, nrow = nBlocks, ncol = nBlocks)  
	adjMat = as.matrix(summary(adjMat))

   for(i in 1:dim(adjMat)[1]) {
		bi = nodeBlocks[adjMat[i,1]]
		bj = nodeBlocks[adjMat[i,2]]
		bMat[bi, bj] = bMat[bi, bj] + 1 
	}

	bMat[upper.tri(bMat)] = (bMat + t(bMat))[upper.tri(bMat)]
	bMat[lower.tri(bMat)] = 0

	for(i in 1:nBlocks) {
    for(j in i:nBlocks) {
      if(i == j){
        bMat[i, i] = bMat[i, i] /
          (nodeCounts[i]^2 - nodeCounts[i])
      }
      else {
        bMat[i, j] = bMat[i, j] /
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
    heatmap.2(bMat, Colv = "Rowv", dendrogram = "none", key = T,
              trace = "none", symm = T,
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

# compute column sd
colSd = function(x, na.rm=TRUE) {
  if (na.rm) {
    n = colSums(!is.na(x))
  } else {
    n = nrow(x)
  }
  colVar =  colMeans(x*x, na.rm=na.rm) - (colMeans(x, na.rm=na.rm))^2
  return(sqrt(colVar * n/(n-1)))
}

# ---------------------------------------------------------------------
# functions for comparison analysis of graph clusters
# ---------------------------------------------------------------------

# get node and edge counts
getNodeEdgeCount = function(preVec, rawDataDir, procDataDir) {

    nGraphs = length(preVec)
    nNodes = vector(length = nGraphs)
    nEdges = vector(length = nGraphs)

    for(i in 1:nGraphs) {
        graphInputFile = paste(rawDataDir, preVec[i], '_big_graph.mat', sep='')
        lccCovInputFile = paste(procDataDir, preVec[i], '_big_lcc.txt', sep='')

        # read file with largest connected component and coordinates
        coorData = read.table(lccCovInputFile)
        
        # read in connectome data
        fiberGraph = readMat(graphInputFile)$fibergraph

        # only keep the largest connected component
        fiberGraph = fiberGraph[coorData$V1 + 1, coorData$V1 + 1]
        nNodes[i] = length(coorData$V1)

        # symmetrize the matrix
        fiberGraph = forceSymmetric(fiberGraph)
        
        # compute number of edges
        rSums = rowSums(fiberGraph)
        nEdges[i] = sum(rSums)/2
    }
    
    return(list(nNodes, nEdges))
}

# get block prob and mean/sd estimates
writeEstBX = function(preVec, procDataDir, rawDataDir, outDir) {

    nGraphs = length(preVec)
    bList = list()
    meanList = list()
    sdList = list()

    for(i in 1:nGraphs) {
        if(!file.exists(paste(outDir, preVec[i], '_estimatedBX.bin', sep=''))) {
            graphInputFile = paste(rawDataDir, preVec[i], '_big_graph.mat',
                sep='')
            lccCovInputFile = paste(procDataDir, preVec[i], '_big_lcc.txt',
                sep='')
            procCoordFile = paste(procDataDir, preVec[i], '_proc_coord.txt',
                sep='')

            # read file with largest connected component and coordinates
            coorData = read.table(lccCovInputFile)
            
            # read in connectome data
            fiberGraph = readMat(graphInputFile)$fibergraph

            # only keep the largest connected component
            fiberGraph = fiberGraph[coorData$V1 + 1, coorData$V1 + 1]

            # symmetrize the matrix
            fiberGraph = forceSymmetric(fiberGraph)
            
            # load cluster assignments
            cascCluster = as.vector(loadMatrix(paste(outDir, preVec[i], "_CASC",
                sep=""), 1))

            # estimate block matrix
            bList[[i]] = estBlockMatSparse(fiberGraph, cascCluster)

            # estimate block location mean and sd
            nClust = length(unique(cascCluster))
            clusterMeans = matrix(0, ncol = 3, nrow = nClust)
            clusterSd = matrix(0, ncol = 3, nrow = nClust)
            coorMat = as.matrix(coorData)[,2:4] + 1

            for(j in 1:nClust) {
                clusterMeans[i,] = colMeans(coorMat[cascCluster[[i]] == j, ])
                clusterSd[i,] = colSd(coorMat[cascCluster[[i]] == j, ])
            }

            meanList[[i]] = clusterMeans
            sdList[[i]] = clusterSd

            # write results
            saveMatrixList(paste(outDir, preVec[i], '_estimatedBX', sep=''),
                           list(bList[[i]], meanList[[i]], sdList[[i]] ))
        }
    }

    return(true)
}

# compute the cluster alignment matrix for brain graph pair
hDist = function(mean1, mean2, sd1, sd2) {

	nClust = dim(mean1)[1]
	dMat = matrix(0, nrow = nClust, ncol = nClust)

	for(i in 1:nClust) {
		for(j in i:nClust) {
			dMat[i,j] = sqrt(sum((mean1[i,] - mean2[j,])^2
                    / (sd1[i,]^2+sd2[j,]^2)))
		}
	}
    dMat[lower.tri(dMat)] = dMat[upper.tri(dMat)]

	return(dMat)
}

# compute the cluster alignment matrix for all brain graphs
writeHDistAll = function(preVec, outDir) {

    nGraphs = length(preVec)
    
    for(i in 1:nGraphs) {

        if(!file.exists(paste(outDir, preVec[i], '_hDistMat.bin', sep=''))) {
            mean1 = loadMatrix(paste(outDir, preVec[i], '_estimatedBX',
                sep=''), 2)
            sd1 = loadMatrix(paste(outDir, preVec[i], '_estimatedBX',
                sep=''), 3)
            hDistList = list()

            for(j in 1:nGraphs) {
                mean2 = loadMatrix(paste(outDir, preVec[j], '_estimatedBX',
                    sep=''), 2)
                sd2 = loadMatrix(paste(outDir, preVec[j], '_estimatedBX',
                    sep=''), 3)

                hDistList[[j]] = hDist(mean1, mean2, sd1, sd2)
            }

            # write results
            saveMatrixList(paste(outDir, preVec[i], '_hDistMat', sep=''),
                       hDistList)
        }

    }

    return(hDistList)
}

# match clusters for a pair of brain graphs using the Hungarian algo
matchPairClusters = function(hDist, id1, id2) {
    return( as.vector(solve_LSAP(hDist)) )
}

# get Frobenius norm of a pair of brain graph block prob.
calcFrobNormB = function(bMat1, bMat2, permVec) {
    return( sqrt(sum((bMat1 - bMat2[permVec, permVec])^2)) )
}

# get Frobenius norm of all pairs of brain graphs
getFrobNormBAll = function(preVec, outDir) {
    nGraphs = length(preVec)
    fMat = matrix(0, nrow = nGraphs, ncol = nGraphs)

    for(i in 1:(nGraphs-1)) {
        bMati = loadMatrix(paste(outDir, preVec[i], '_estimatedBX',
            sep=''), 1)
        
        for(j in (i+1):nGraphs) {
            hDist = loadMatrix(paste(outDir, preVec[i], '_hDistMat', sep=''), j)
            permVec = matchPairClusters(hDist, i, j)

            bMatj = loadMatrix(paste(outDir, preVec[j], '_estimatedBX',
                sep=''), 1)
       
            fMat[i, j] = calcForbNormB(bMati, bMatj, permVec)
        }
    }

    return(fMat)
}
