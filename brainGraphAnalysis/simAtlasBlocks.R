# ---------------------------------------------------------------------
# load packages and source files
# ---------------------------------------------------------------------
#if (!require(rARPACK)) {
#    install.packages('rARPACK', dependencies = T)
#    require(rARPACK)
#}
if (!require(irlba)) {
    install.packages('irlba', dependencies = T)
    require(irlba)
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
if (!require(mclust)) {
    install.packages('mclust', dependencies = T)
    require(mclust)
}

initialDir = getwd()
setwd('../code/')
source('simulateGraphData.R')
source('readWriteMatrix.R')
source('extEvalMethods.R')
source('spectralClusteringMethods.R')
source('helperFunctions.R')
setwd(initialDir)

# ---------------------------------------------------------------------
# load block matrix and covariate means 
# ---------------------------------------------------------------------
comArgs = commandArgs(T)
filePre = comArgs[1]
nCores = as.numeric(comArgs[2])
nH = as.numeric(comArgs[3])

outDir = 'cache/'
procDataDir = '../procData/'

# block matrix file
blockMatFile = paste(outDir, filePre, '_estimatedB', sep='')

bMat = loadMatrix(blockMatFile, 1)

# covariate parameters file with node counts
covParamFile = paste(outDir, filePre, '_covParamEst', sep='')

clusterMeans = loadMatrix(covParamFile, 1)
clusterSd = loadMatrix(covParamFile, 2)
nodeClusters = loadMatrix(covParamFile, 3)

# ---------------------------------------------------------------------
# run simulations
# ---------------------------------------------------------------------
nIter = 10
nBlocks = dim(bMat)[1]
nCov = dim(clusterMeans)[2]

misRateSc = vector(length = nIter)
misRateCasc = vector(length = nIter)
misRateCca = vector(length = nIter)
misRateScx = vector(length = nIter)

ariSc = vector(length = nIter)
ariCasc = vector(length = nIter)
ariCca = vector(length = nIter)
ariScx = vector(length = nIter)

nNodes = length(nodeClusters)
nMembers = table(nodeClusters)
nodeMembership = rep(1:nBlocks, times = nMembers)

# set number of cores
registerDoMC(nCores)

foreach(i = 1:nIter) %dopar% {

        adjMat = simSparseAdjMat(bMat, nMembers)

        # ensure graph is connnected
        while( min(rowSums(adjMat)) == 0) {
            adjMat = simSparseAdjMat(bMat, nMembers)
        }

        # compute regularized graph Laplacian
        rSums = rowSums(adjMat)
        tau = sum(rSums)/length(rSums)
        adjMat = Diagonal(nNodes, 1/sqrt(rSums + tau)) %*% adjMat %*%
            Diagonal(nNodes, 1/sqrt(rSums + tau))
        
        # compute svd's of L and X
        # lapSvd = eigs(adjMat, nBlocks + 1, opts = list(maxitr = 10000))
        # compute svd's using irlba since we are using this for CASC
        lapSvd = irlba(adjMat, nu = nBlocks + 1, m_b = 2*(nBlocks + 1)) 
        lapSvd$u = lapSvd$u[ ,1:nBlocks]

        #save the result matrix
        saveMatrixList(paste(outDir, 'simLapSvd_', i, sep=''),
                       list(lapSvd$u, as.matrix(lapSvd$d)))
        saveSparseMatrix(paste(outDir, 'simAdjMat_', i, sep=''), adjMat)
}

for(i in 1:nIter) {
        #simulate location from normal dist
        #coordMat = simCoordMat(clusterMeans, clusterSd, nMembers)

        #use fixed location
        coordMat = loadMatrix(paste(procDataDir, filePre, '_big_lcc_sim',
            sep=''), 1)

        #permute fixed location for correct block membership
        coordMat = coordMat[sort(nodeClusters, index.return=T)$ix,]
        
        # add normal noise
        #coordMat = coordMat + cbind(rnorm(nNodes, 0, .25),
        #    rnorm(nNodes, 0, .25), rnorm(nNodes, 0, .25))
            
        # demean the columns
        coordMat = coordMat - rep(1, nNodes) %*%
            t(colMeans(coordMat)) 

        # rescale
        coordMat = coordMat/max(coordMat)
        covSvd = svd(coordMat)

        # load matrix
        lapSvd = list()
        lapSvd$u = loadMatrix(paste(outDir, 'simLapSvd_', i, sep=''), 1)
        lapSvd$d = as.vector(loadMatrix(
            paste(outDir, 'simLapSvd_', i, sep=''), 2))

        adjMat = loadSparseMatrix(paste(outDir, 'simAdjMat_', i, sep=''))
        
        # compute upper and lower bounds for h
        hMin = (lapSvd$d[nBlocks] - lapSvd$d[nBlocks+1])/covSvd$d[1]^2
        hMax = lapSvd$d[1]/covSvd$d[min(nBlocks, nCov)]

        # ---------------------------------------------------------------------
        # for comparison compute RSC, CCA and spectral clustering on X
        # ---------------------------------------------------------------------
        # do regularized spectral clustering for comparison
        scSv = lapSvd$u/sqrt(rowSums(lapSvd$u^2))
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

        # compute and store the ARI for each method
        ariSc[i] = adjustedRandIndex(scCluster, nodeMembership)
        ariCasc[i] = adjustedRandIndex(cascCluster, nodeMembership)
        ariCca[i] = adjustedRandIndex(ccaCluster, nodeMembership)
        ariScx[i] = adjustedRandIndex(scxCluster, nodeMembership)
        
        #track progress
        write.table(t(c(misRateSc[i], misRateCasc[i], misRateCca[i],
                      misRateScx[i])), append = T, row.names = F,
                    col.names = F, paste(outDir, 'logSimAtlasMisRate.txt',
                        sep=''))
        
        write.table(t(c(ariSc[i], ariCasc[i], ariCca[i],
                        ariScx[i])), append = T, row.names = F,
                    col.names = F, paste(outDir, 'logSimAtlasAri.txt',
                        sep=''))
}

