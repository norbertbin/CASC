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

outDir = 'cache/'

# block matrix file
blockMatFile = paste(outDir, filePre, '_estimatedB', sep='')

bMat = loadMatrix(blockMatFile, 1)

# covariate parameters file with node counts
covParamFile = paste(outDir, filePre, '_covParamEst', sep='')

clusterMeans = loadMatrix(covParamFile, 1)
clusterSd = loadMatrix(covParamFile, 2)
nodeCounts = loadMatrix(covParamFile, 3)

# ---------------------------------------------------------------------
# run simulations
# ---------------------------------------------------------------------
nIter = 20
nBlocks = dim(bMat)[1]
nCov = dim(clusterMeans)[2]

misRateSc = vector(length = nIter)
misRateCasc = vector(length = nIter)
misRateCca = vector(length = nIter)
misRateScx = vector(length = nIter)


nNodes = sum(nodeCounts)
nMembers = nodeCounts

for(i in 1:nIter) {
        # set number of cores
        registerDoMC(nCores)

        adjMat = simSparseAdjMat(bMat, nMembers)
        coordMat = simCoordMat(clusterMeans, clusterSd, nMembers)

        # add normal noise
        coordMat = coordMat + cbind(rnorm(nNodes, 0, .25),
            rnorm(nNodes, 0, .25), rnorm(nNodes, 0, .25))
            
        # demean the columns
        coordMat = coordMat - rep(1, nNodes) %*%
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
        misRateSc[i] = misClustRateClust(scCluster, nMembers)
        misRateCasc[i] = misClustRateClust(cascCluster, nMembers)
        misRateCca[i] = misClustRateClust(ccaCluster, nMembers)
        misRateScx[i] = misClustRateClust(scxCluster, nMembers)

        #track progress
        write.table(c(misRateSc[i], misRateCasc[i], misRateCca[i],
                      misRateScx[i]), append = T, row.names = F,
                    col.names = F, paste(outDir, 'logSimAtlasBlocks.txt',
                        sep=''))

}

# save misclustering results
misRates = data.frame(casc = rowSums(misRateCasc)/nIter,
    cca = rowSums(misRateCca)/nIter,
    sc = rowSums(misRateSc)/nIter,
    scx = rowSums(misRateScx)/nIter)

write.table(misRates, paste(outDir, 'simAtlasBlocks.txt', sep='')) 

