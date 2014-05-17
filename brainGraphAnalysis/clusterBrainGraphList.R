# ---------------------------------------------------------------------
# load packages and source files
# ---------------------------------------------------------------------
if (!require(Matrix)) {
    install.packages('Matrix', dependencies = T)
    require(Matrix)
}
if (!require(R.matlab)) {
    install.packages('R.matlab', dependencies = T)
    require(R.matlab)
}
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

cwd = getwd()
setwd("../code/")
source("readWriteMatrix.R")
source("irlbaMod.R")
source("spectralClusteringMethods.R")
setwd(cwd)

# ---------------------------------------------------------------------
# get command line arguements with data directories
# ---------------------------------------------------------------------
comArgs = commandArgs(T)
rawDataDir = comArgs[1]
procDataDir = comArgs[2]
preListFile = comArgs[3]
nBlocks = as.numeric(comArgs[4]) #num of blocks/clusters in model
nCores = as.numeric(comArgs[5]) #max num of cores that can be used
nH = as.numeric(comArgs[6]) #num of h values to use in CASC

# define other directories
outDir = "cache/"
figDir = "figs/"

# read in file prefixes
preVec = readLines(preListFile)

# number of graphs
nGraphs = length(preVec)

# set constant number of coordinate covariates
nCov = 3

# number of cores 
registerDoMC(nCores)

# ---------------------------------------------------------------------
# run RSC, SCX, CCA in parallel
# ---------------------------------------------------------------------
hList = foreach(i = 1:nGraphs) %dopar% {
    # input files
    graphInputFile = paste(rawDataDir, preVec[i], '_big_graph.mat', sep='')
    lccCovInputFile = paste(procDataDir, preVec[i], '_big_lcc.txt', sep='')
    procCoordFile = paste(procDataDir, preVec[i], '_proc_coord.txt', sep='')

    if(!file.exists(procCoordFile)) {
        # read file with largest connected component and coordinates
        coorData = read.table(lccCovInputFile)

        nNodes = length(coorData$V1)
        
        # create covariate matrix with centered xyz and rescaled to 1 
        # add white noise with sigma = 1/4
        covData = cbind(coorData$V2, coorData$V3, coorData$V4) +
            cbind(rnorm(nNodes, 0, 1/4), rnorm(nNodes, 0, 1/4),
                  rnorm(nNodes, 0, 1/4))
        covData = covData - rep(1, nNodes) %*% t(colMeans(covData))
        covData = covData/max(covData)
        
        # store coordinate data with added noise
        write.table(covData, procCoordFile)
    } else {
        covData = as.matrix(read.table(procCoordFile, header = T))
    }

    # ---------------------------------------------------------------------
    # compute SVD to find tuning range
    # ---------------------------------------------------------------------

    if(!file.exists(paste(outDir, preVec[i], "_LSVD.bin", sep="")) |
       !file.exists(paste(outDir, preVec[i], "_covSVD.bin", sep=""))) {

        # read file with largest connected component and coordinates
        coorData = read.table(lccCovInputFile)
        
        # read in connectome data
        fiberGraph = readMat(graphInputFile)$fibergraph

        # only keep the largest connected component
        fiberGraph = fiberGraph[coorData$V1 + 1, coorData$V1 + 1]
        nNodes = length(coorData$V1)

        # symmetrize the matrix
        fiberGraph = forceSymmetric(fiberGraph)
        
        # compute the regularized Laplacian
        rSums = rowSums(fiberGraph)
        tau = mean(rSums)
        fiberGraph = Diagonal(nNodes, 1/sqrt(rSums + tau)) %*% fiberGraph %*%
            Diagonal(nNodes, 1/sqrt(rSums + tau))

        # compute svd's of L and X
        lapSvd = irlba(fiberGraph, nu = nBlocks + 1, m_b = 2*nBlocks)
        covSvd = svd(covData)

        # save these to outDir
        saveMatrixList(paste(outDir, preVec[i], "_LSVD", sep=""), list(
            matrix(lapSvd$d), lapSvd$u[,1:nBlocks]))
        saveMatrixList(paste(outDir, preVec[i], "_covSVD", sep=""), list(
            matrix(covSvd$d), covSvd$u))

    } else {
        # load svd from outDir
        lapSvd = list()
        lapSvd$d = loadMatrix(paste(outDir, preVec[i], "_LSVD", sep=""), 1)
        lapSvd$u = loadMatrix(paste(outDir, preVec[i], "_LSVD", sep=""), 2)
        covSvd = list()
        covSvd$d = loadMatrix(paste(outDir, preVec[i], "_covSVD", sep=""), 1)
        covSvd$u = loadMatrix(paste(outDir, preVec[i], "_covSVD", sep=""), 2)
    }
    
    # ---------------------------------------------------------------------
    # for comparison compute RSC, CCA and spectral clustering on X
    # ---------------------------------------------------------------------
    # do regularized spectral clustering for comparison
    if(!file.exists(paste(outDir, preVec[i], "_SC.bin", sep=""))) {
        scSv = lapSvd$u/sqrt(rowSums(lapSvd$u^2))

        # don't reparallelize 
        registerDoMC(1)

        scKM = bigkmeans(scSv, nBlocks, iter.max = 200, nstart = 10)
        scCluster = scKM$cluster
        saveMatrixList(paste(outDir, preVec[i], "_SC", sep=""), list(
            matrix(scCluster) ))
    } 

    #do CCA
    if(!file.exists(paste(outDir, preVec[i], "_CCA.bin", sep=""))) {
        ccaSvd = svd(fiberGraph%*%covData)

        # number of cores 
        registerDoMC(1)

        ccaKM = bigkmeans(ccaSvd$u, nBlocks, iter.max = 200, nstart = 10)
        ccaCluster = ccaKM$cluster
        saveMatrixList(paste(outDir, preVec[i], "_CCA", sep=""),
               list(matrix(ccaCluster)))
    } 

    #do SC on X
    if(!file.exists(paste(outDir, preVec[i], "_SCX.bin", sep=""))) {

        # number of cores 
        registerDoMC(1)

        scxKM = bigkmeans(covSvd$u, nBlocks, iter.max = 200, nstart = 10)
        scxCluster = scxKM$cluster
        saveMatrixList(paste(outDir, preVec[i], "_SCX", sep=""),
           list(matrix(scxCluster)))
    }

    # compute upper and lower bounds for h
    hMin = (lapSvd$d[nBlocks] - lapSvd$d[nBlocks+1])/covSvd$d[1]^2
    hMax = lapSvd$d[1]/covSvd$d[min(nBlocks, nCov)]

    # foreach output values
    c(hMin, hMax)
}

# do memory clean up
rm(ccaSvd)
rm(covSvd)
rm(lapSvd)
rm(scSv)
rm(coorData)
gc()

for(i in 1:nGraphs) {
    hMin = hList[[i]][1]
    hMax = hList[[i]][2]

    # ---------------------------------------------------------------------
    # compute SVD and clusters for set of tuning parameters of CASC
    # ---------------------------------------------------------------------
    hSet = seq(hMin, hMax, length.out = nH)

    wcssVec = rep(0, nH)
    clusterMat = matrix(rep(0, nH*nNodes), nrow = nH)

    # number of cores 
    registerDoMC(nCores)
    
    kmList = foreach(h = hSet) %dopar% {
        cascSvd = getCascSvd(fiberGraph, covData, h, nBlocks)

        eigengap = cascSvd$d[nBlocks] - cascSvd$d[nBlocks+1]

        #save eigengap and h for reference
        write.table(t(c(h, eigengap)), append = T, row.names = F,
                col.names = F, paste(outDir, preVec[i],
                    '_h_eigenGap.txt', sep=''))
    
        #project eigenvectors onto sphere
        cascSvd$singVec = cascSvd$singVec/sqrt(rowSums(cascSvd$singVec^2))

        registerDoMC(1) #dont reparallelize
        cascKM = bigkmeans(cascSvd$singVec, nBlocks, iter.max = 100,
            nstart = 10)
    }

    for(h in hSet) {
        wcssVec[match(h, hSet)] = mean(kmList[[match(h, hSet)]]$withinss)
        clusterMat[match(h, hSet),] = kmList[[match(h, hSet)]]$cluster
    }

    # ---------------------------------------------------------------------
    # get clusters with min wcss and save to file
    # ---------------------------------------------------------------------
    iMin = match(min(wcssVec), wcssVec)

    # save final h value
    write.table(c(hSet[iMin]), append = T, row.names = F,
            col.names = F, paste(outDir, preVec[i], '_h_eigenGap.txt', sep=''))

    cascCluster = clusterMat[iMin,]
    saveMatrixList(paste(outDir, preVec[i], "_CASC", sep=""), list(
        matrix(cascCluster) ))
}

