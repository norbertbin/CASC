# ---------------------------------------------------------------------
# load packages and source files
# ---------------------------------------------------------------------
if (!require(rARPACK)) {
    install.packages('rARPACK', dependencies = T)
    require(rARPACK)
}
if (!require(foreach)) {
    install.packages('foreach', dependencies = T)
    require(foreach)
}
if (!require(doMC)) {
    install.packages('doMC', dependencies = T)
    require(doMC)
}
if (!require(biganalytics)) {
    install.packages('biganalytics', dependencies = T)
    require(biganalytics)
}

initialDir = getwd()
setwd('../code/')
source('simulateGraphData.R')
source('readWriteMatrix.R')
source('extEvalMethods.R')
source('spectralClusteringMethods.R')
setwd(initialDir)
# ---------------------------------------------------------------------
# load block matrix and covariate means 
# ---------------------------------------------------------------------
comArgs = commandArgs(T)
filePre = comArgs[1]
nCores = as.numeric(comArgs[2])
nH = as.numeric(comArgs[3])

# set number of cores
registerDoMC(nCores)

outDir = 'cache/'

# block matrix file
blockMatFile = paste(outDir, filePre, '_estimatedB')

bMat = loadMatrix(blockMatFile, 1)

# covariate parameters file
covParamFile = paste(outDir, filePre, '_covParamEst')

clusterMeans = loadMatrix(covParamFile, 1)
clusterSd = loadMatrix(covParamFile, 2)
nodeCounts = loadMatrix(covParamFile, 3)

# ---------------------------------------------------------------------
# run simulations
# ---------------------------------------------------------------------
nIter = 20
nBlocks = dim(bMat)[1]
nCov = dim(clusterMeans)[2]
nNodesSeq = c(3*10^4, 6*10^4, 9*10^4)

#adjust block matrix to ensure all nodes connected
#also increase signal in covariates by reducing clusterSd
bMat = bMat*(sum(nodeCounts)/nNodesSeq[1]) /
    (log(sum(nodeCounts))/log(nNodesSeq[1]))
clusterSd = clusterSd/3

misRateSc = matrix(0, nrow = length(nNodesSeq), ncol = nIter)
misRateCasc = misRateSc
misRateCca = misRateSc
misRateScx = misRateSc

for(i in 1:length(nNodesSeq)) {
    nNodes = nNodesSeq[i]
    nMembers = round(nNodes * nodeCounts/sum(nodeCounts))
    nMembers[which(nMembers == max(nMembers))] = nMembers[which(nMembers ==
                max(nMembers))] + nNodes - sum(nMembers)#hack to keep count
    
    for(j in 1:nIter) {
        adjMat = simSparseAdjMat(bMat, nMembers)
        coordMat = simCoordMat(clusterMeans, clusterSd, nMembers)

        # add normal noise
        coordMat = coordMat + cbind(rnorm(nNodes, 0, .25),
            rnorm(nNodes, 0, .25), rnorm(nNodes, 0, .25))
            
        # demean the columns
        coordMat = coordMat - vector(length = nNodes) %*%
            t(colMeans(coordMat)) 

        # rescale
        coordMat = coordMat/max(coordMat)

        # compute regularized graph Laplacian
        rSums = rowSums(adjMat)
        tau = sum(rSums)/length(rSums)
        adjMat = Diagonal(nNodes, 1/sqrt(rSums + tau)) %*% adjMat %*%
            Diagonal(nNodes, 1/sqrt(rSums + tau))
        
        # compute svd's of L and X
        lapSvd = eigs(adjMat, nBlocks + 1)
        lapSvd$vectors = lapSvd$vectors[,1:nBlocks]
        covSvd = svd(coordMat)

        # compute upper and lower bounds for h
        hMin = (lapSvd$values[nBlocks] - lapSvd$values[nBlocks+1])/covSvd$d[1]^2
        hMax = lapSvd$values[1]/covSvd$d[min(nBlocks, nCov)]

        # ---------------------------------------------------------------------
        # for comparison compute RSC, CCA and spectral clustering on X
        # ---------------------------------------------------------------------
        # do regularized spectral clustering for comparison
        scSv = lapSvd$vectors/sqrt(rowSums(lapSvd$vectors^2))
        scKM = bigkmeans(scSv, nBlocks, iter.max = 200, nstart = 10)
        scCluster = scKM$cluster

        # do CCA
        ccaSvd = svd(adjMat%*%coordMat)
        ccaKM = bigkmeans(ccaSvd$u, nBlocks, iter.max = 200, nstart = 10)
        ccaCluster = ccaKM$cluster

        # do SC on X
        scxKM = bigkmeans(covSvd$u, nBlocks, iter.max = 200, nstart = 10)
        scxCluster = scxKM$cluster

        # ---------------------------------------------------------------------
        # compute SVD and clusters for set of tuning parameters of CASC
        # ---------------------------------------------------------------------
        hSet = seq(hMin, hMax, length.out = nH)

        wcssVec = rep(0, nH)
        clusterMat = matrix(rep(0, nH*nNodes), nrow = nH)

        kmList = foreach(h = hSet) %dopar% {
            cascSvd = getCascSvd(adjMat, coordMat, h, nBlocks)$singVec

            # project eigenvectors onto sphere
            cascSvd = cascSvd/sqrt(rowSums(cascSvd^2))

            registerDoMC(1) #dont reparallelize
            cascKM = bigkmeans(cascSvd, nBlocks, iter.max = 100, nstart = 10)

        }

        for(h in hSet) {
            wcssVec[match(h, hSet)] = mean(kmList[[match(h, hSet)]]$withinss)
            clusterMat[match(h, hSet),] = kmList[[match(h, hSet)]]$cluster
        }

        # get clusters with min wcss
        iMin = match(min(wcssVec), wcssVec)
        cascCluster = clusterMat[iMin,]

        # compute and store the misclustering rate for each method
        misRateSc[i, j] = misClustRateClust(scCluster, nMembers)
        misRateCasc[i, j] = misClustRateClust(cascCluster, nMembers)
        misRateCca[i, j] = misClustRateClust(ccaCluster, nMembers)
        misRateScx[i, j] = misClustRateClust(scxCluster, nMembers)
    }
}

# save misclustering results
misRates = data.frame(n = rep(nNodesSeq, 4), misClustRate =
    c(rowSums(misRateCasc)/nIter, rowSums(misRateCca)/nIter,
    rowSums(misRateSc)/nIter, rowSums(misRateScx)/nIter),
    group = rep(1:4, each = nPoints))

write.table(misRates, paste(outDir, 'simAtlasBlocks.txt')) 

# plot misclustering results
legLab = c("CASC", "CCA", "RSC", "SC on X")
pdf(paste(outDir, "simAtlasBlocks.pdf"), width = 7, height = 7)
print(
    xyplot(misClustRate ~ n, group = group, type = "b", pch = 1:4,
           cex = 1.2, data = misRates, ylab = "Average mis-clustering rate",
           xlab = "Number of nodes (N)", lwd = 2, key = list(
                               text = list(legLab),
                               lines = list(col = 
			trellis.par.get()$superpose.symbol$col[1:4], lwd = 2),
                               points = list(pch = 1:4, cex = 1.2,
			col = trellis.par.get()$superpose.symbol$col[1:4]),
                                        corner = c(.95, .95)) )
    )
dev.off()
