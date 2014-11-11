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

# ---------------------------------------------------------------------
# define directories
# load brain graph data
# ---------------------------------------------------------------------
outDir = "cache/"
figDir = "figs/"

# input/output file
covFile = paste(outDir, filePre, '_big_graph_w_inv_attr.txt', sep='')

#set number of cores to be used for parallel methods
registerDoMC(nCores)

# read in connectome data
fiberGraph = readMat(graphInputFile)$fibergraph

# read file with largest connected component and coordinates
coorData = read.table(lccCovInputFile)

# read file with atlas labels
covData = as.matrix(read.table(covFile))

# construct the sparse model matrix using the atlas labels
covData = Matrix(model.matrix(~., data = data.frame(as.factor(covData))))

# set constant number of covariates
nCov = dim(covData)[2]

# ---------------------------------------------------------------------
# preprocess data 
# ---------------------------------------------------------------------

# only keep the largest connected component
fiberGraph = forceSymmetric(fiberGraph[coorData$V1 + 1, coorData$V1 + 1])
nNodes = length(coorData$V1)

# center and scale covariate data
#covData = scale(covData, scale = F)
#covData = covData/sqrt(sum(covData^2))

# compute the regularized Laplacian
rSums = rowSums(fiberGraph)
tau = mean(rSums)
fiberGraph = forceSymmetric(Diagonal(nNodes, 1/sqrt(rSums + tau))
    %*% fiberGraph %*% Diagonal(nNodes, 1/sqrt(rSums + tau)))

# ---------------------------------------------------------------------
# compute SVD to find tuning range
# ---------------------------------------------------------------------

if(!file.exists(paste(outDir, filePre, "_LSVD.bin", sep=""))) {
    # compute svd's of L and X
    lapSvd = irlba(fiberGraph, nu = nBlocks + 1, m_b = 2*nBlocks)

    # save these to outDir
    saveMatrixList(paste(outDir, filePre, "_LSVD", sep=""), list(
        matrix(lapSvd$d), lapSvd$u[,1:nBlocks]))

} else {
    # load svd from outDir
    lapSvd = list()
    lapSvd$d = loadMatrix(paste(outDir, filePre, "_LSVD", sep=""), 1)
    lapSvd$u = loadMatrix(paste(outDir, filePre, "_LSVD", sep=""), 2)
}

# compute svd of sparse model matrix
covSvd = svd(covData)

# compute upper and lower bounds for h
hMin = (lapSvd$d[nBlocks] - lapSvd$d[nBlocks+1])/covSvd$d[1]^2
hMax = lapSvd$d[1]/covSvd$d[min(nBlocks, nCov)]

print(paste("Range:", hMin, "to", hMax))

# do memory clean up
rm(lapSvd)
rm(coorData)
gc()

# ---------------------------------------------------------------------
# compute SVD and clusters for set of tuning parameters of CASC
# ---------------------------------------------------------------------
hSet = seq(hMin, hMax, length.out = nH)

wcssVec = rep(0, nH)
clusterMat = matrix(rep(0, nH*nNodes), nrow = nH)

kmList = foreach(h = hSet) %dopar% {
    cascSvd = getCascSvd(fiberGraph, covData, h, nBlocks)

    eigengap = cascSvd$singVal[nBlocks] - cascSvd$singVal[nBlocks+1]
    
    #project eigenvectors onto sphere
    cascSvd$singVec = cascSvd$singVec/sqrt(rowSums(cascSvd$singVec^2))

    registerDoMC(1) #dont reparallelize
    cascKM = bigkmeans(cascSvd$singVec, nBlocks, iter.max = 100, nstart = 10)
    returnVal = list(cluster = cascKM$cluster, withinss = cascKM$withinss,
        eigengap = eigengap, h = h)
}

#save list of lists with cluster information as an R object
save(kmList, file = paste(outDir, filePre, "_CASC_AtlasClust.RData", sep=""))

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
            col.names = F, paste(outDir, filePre, '_h_eigenGap_Atlas.txt', sep=''))

cascCluster = clusterMat[iMin,]
saveMatrixList(paste(outDir, filePre, "_CASC_Atlas", sep=""), list(
    matrix(cascCluster) ))


