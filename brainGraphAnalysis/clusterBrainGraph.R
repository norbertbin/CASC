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
filePre = comArgs[3]
nBlocks = as.numeric(comArgs[4]) #num of blocks/clusters in model
nCores = as.numeric(comArgs[5]) #max num of cores that can be used
nH = as.numeric(comArgs[6]) #num of h values to use in CASC

# input files
graphInputFile = paste(rawDataDir, filePre, '_big_graph.mat', sep='')
lccCovInputFile = paste(procDataDir, filePre, '_big_lcc.txt', sep='')
print(graphInputFile)
print(lccCovInputFile)

# input/output files
procCoordFile = paste(procDataDir, filePre, '_proc_coord.txt', sep='')

# ---------------------------------------------------------------------
# define directories
# load brain graph data
# ---------------------------------------------------------------------
outDir = "cache/"
figDir = "figs/"

# set constant number of coordinate covariates
nCov = 3

#set number of cores to be used for parallel methods
registerDoMC(nCores)

# read in connectome data
fiberGraph = readMat(graphInputFile)$fibergraph

# read file with largest connected component and coordinates
coorData = read.table(lccCovInputFile)

# ---------------------------------------------------------------------
# preprocess data 
# ---------------------------------------------------------------------

# only keep the largest connected component
fiberGraph = fiberGraph[coorData$V1 + 1, coorData$V1 + 1]
nNodes = length(coorData$V1)

# symmetrize the matrix
fiberGraph = fiberGraph + t(fiberGraph)

if(!file.exists(procCoordFile)) {
# create covariate matrix with centered xyz and rescaled to 1 
# add white noise with sigma = 1/4
    coorData$V2 = coorData$V2 + rnorm(nNodes, 0, 1/4)
    coorData$V3 = coorData$V3 + rnorm(nNodes, 0, 1/4)
    coorData$V4 = coorData$V4 + rnorm(nNodes, 0, 1/4)
    maxCoor = max(c(coorData$V2, coorData$V3, coorData$V4))
    covData = matrix(c(coorData$V2 - mean(coorData$V2),
        coorData$V3 - mean(coorData$V3),
        coorData$V4 - mean(coorData$V4)), ncol = 3)/maxCoor

# store coordinate data with added noise
    write.table(covData, procCoordFile)
} else {
    covData = as.matrix(read.table(procCoordFile, header = T))
}

# compute the regularized Laplacian
rSums = rowSums(fiberGraph)
tau = mean(rSums)
fiberGraph = Diagonal(nNodes, 1/sqrt(rSums + tau)) %*% fiberGraph %*%
    Diagonal(nNodes, 1/sqrt(rSums + tau))

# ---------------------------------------------------------------------
# compute SVD to find tuning range
# ---------------------------------------------------------------------

if(!file.exists(paste(outDir, filePre, "_LSVD.bin", sep="")) |
   !file.exists(paste(outDir, filePre, "_covSVD.bin", sep=""))) {
    # compute svd's of L and X
    lapSvd = irlba(fiberGraph, nu = nBlocks + 1, m_b = 2*nBlocks)
    covSvd = svd(covData)

    # save these to outDir
    saveMatrixList(paste(outDir, filePre, "_LSVD", sep=""), list(
        matrix(lapSvd$d), lapSvd$u[,1:nBlocks]))
    saveMatrixList(paste(outDir, filePre, "_covSVD", sep=""), list(
        matrix(covSvd$d), covSvd$u))

} else {
    # load svd from outDir
    lapSvd = list()
    lapSvd$d = loadMatrix(paste(outDir, filePre, "_LSVD", sep=""), 1)
    lapSvd$u = loadMatrix(paste(outDir, filePre, "_LSVD", sep=""), 2)
    covSvd = list()
    covSvd$d = loadMatrix(paste(outDir, filePre, "_covSVD", sep=""), 1)
    covSvd$u = loadMatrix(paste(outDir, filePre, "_covSVD", sep=""), 2)
}

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
saveMatrixList(paste(outDir, filePre, "_SC", sep=""), list(
    matrix(scCluster) ))

#do CCA
ccaSvd = svd(fiberGraph%*%covData)
ccaKM = bigkmeans(ccaSvd$u, nBlocks, iter.max = 200, nstart = 10)
ccaCluster = ccaKM$cluster
saveMatrixList(paste(outDir, filePre, "_CCA", sep=""),
               list(matrix(ccaCluster)))

#do SC on X
scxKM = bigkmeans(covSvd$u, nBlocks, iter.max = 200, nstart = 10)
scxCluster = scxKM$cluster
saveMatrixList(paste(outDir, filePre, "_SCX", sep=""),
           list(matrix(scxCluster)))

# ---------------------------------------------------------------------
# compute SVD and clusters for set of tuning parameters of CASC
# ---------------------------------------------------------------------
hSet = seq(hMin, hMax, length.out = nH)

wcssVec = rep(0, nH)
clusterMat = matrix(rep(0, nH*nNodes), nrow = nH)

kmList = foreach(h = hSet) %dopar% {
    cascSvd = getCascSvd(fiberGraph, covData, h, nBlocks)$singVec

    #project eigenvectors onto sphere
    cascSvd = cascSvd/sqrt(rowSums(cascSvd^2))

    registerDoMC(1) #dont reparallelize
    cascKM = bigkmeans(cascSvd, nBlocks, iter.max = 100, nstart = 10)

}

for(h in hSet) {
    wcssVec[match(h, hSet)] = mean(kmList[[match(h, hSet)]]$withinss)
    clusterMat[match(h, hSet),] = kmList[[match(h, hSet)]]$cluster
}

# ---------------------------------------------------------------------
# get clusters with min wcss and save to file
# ---------------------------------------------------------------------
iMin = match(min(wcssVec), wcssVec)
cascCluster = clusterMat[iMin,]
saveMatrixList(paste(outDir, filePre, "_CASC", sep=""), list(
    matrix(cascCluster) ))


