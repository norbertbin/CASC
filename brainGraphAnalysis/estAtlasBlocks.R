# ---------------------------------------------------------------------
# load packages and source files
# ---------------------------------------------------------------------
if (!require(oro.nifti)) {
    install.packages('oro.nifti', dependencies = T)
    require(oro.nifti)
}

if (!require(reshape2)) {
    install.packages('reshape2', dependencies = T)
    require(reshape2)
}

if (!require(mclust)) {
    install.packages('mclust', dependencies = T)
    require(mclust)
}

if (!require(R.matlab)) {
    install.packages('R.matlab', dependencies = T)
    require(R.matlab)
}

if (!require(Matrix)) {
    install.packages('Matrix', dependencies = T)
    require(Matrix)
}

source('../code/helperFunctions.R')
source('../code/readWriteMatrix.R')

# ---------------------------------------------------------------------
# load atlas and other data
# ---------------------------------------------------------------------
comArgs = commandArgs(T)
atlasFileName = comArgs[1]
rawDataDir = comArgs[2]
procDataDir = comArgs[3]
filePre = comArgs[4]

#define other directories
outDir = "cache/"
figDir = "figs/"

#load graph
graphInputFile = paste(rawDataDir, filePre, '_big_graph.mat', sep='')
fiberGraph = readMat(graphInputFile)$fibergraph

#load atlas
atlasDat = readNIfTI(atlasFileName)

#load clusters from cache file
cascInputFile = paste(outDir, filePre, '_CASC', sep='')
cascCluster = loadMatrix(cascInputFile, 1)
scInputFile = paste(outDir, filePre, '_SC', sep='')
scCluster = loadMatrix(scInputFile, 1)

#load graph coordinates
coordInputFile = paste(procDataDir, filePre, '_big_lcc.txt', sep='')
coorData = read.table(coordInputFile)

# only keep the largest connected component
fiberGraph = fiberGraph[coorData$V1 + 1, coorData$V1 + 1]

#create coordinate matrix starting at 1 not 0
coorMat = as.matrix(coorData)[,2:4] + 1

# ---------------------------------------------------------------------
# convert data tensor to coordinate label matrix
# match atlas to graph
# ---------------------------------------------------------------------
atlasMat = melt(atlasDat@.Data)
atlasMat = atlasMat[atlasMat$value > 0,]

# match atlas with data using location coordinates
indexAM = cbind(atlasMat$Var1, atlasMat$Var2, atlasMat$Var3) %*%
    c(1, 10^3, 10^6)
indexCM = coorMat %*% c(1, 10^3, 10^6)
matchIndex = match(indexCM, indexAM)
sum(is.na(matchIndex))

# atlas assignments for graph
atlasCluster = atlasMat$value[matchIndex]

# compute ari for atlas vs casc and sc
adjustedRandIndex(atlasCluster, cascCluster)
adjustedRandIndex(atlasCluster, scCluster)

# ---------------------------------------------------------------------
# estimate block matrix using atlas and graph
# ---------------------------------------------------------------------
# exclude very small clusters
nodeCounts = tabulate(atlasCluster)
smallClusters = which(nodeCounts < 500)

# set na's to new cluster
idRemove = c( which(is.na(atlasCluster)), which(atlasCluster %in%
    smallClusters))

# estimate block matrix
bMat = estBlockMat(fiberGraph[-idRemove, -idRemove], atlasCluster[-idRemove])

# save block matrix in cache
saveMatrixList(paste(outDir, filePre, '_estimatedB', sep=''), list(bMat))

# plot estimate of block matrix (gives margin error)
# plotB(bMat, paste(figDir, filePre, '_estimatedB.pdf', sep=''))

# ---------------------------------------------------------------------
# estimate covariate means and sd
# ---------------------------------------------------------------------
retainedClusters = unique(atlasCluster[-idRemove])
clusterMeans = matrix(0, ncol = 3, nrow = length(retainedClusters))
clusterSd = matrix(0, ncol = 3, nrow = length(retainedClusters))

# change na's to 0 to avoid errors
atlasCluster[is.na(atlasCluster)] = 0

for(i in retainedClusters) {
    clusterMeans[which(retainedClusters == i),] =
        colMeans(coorMat[atlasCluster == i, ])
    clusterSd[which(retainedClusters == i),] =
        colSd(coorMat[atlasCluster == i, ])
}

# save matrix and vector of block node counts           
saveMatrixList(paste(outDir, filePre, '_covParamEst', sep=''),
               list(clusterMeans, clusterSd,
                    matrix(table(atlasCluster[-idRemove]))))
