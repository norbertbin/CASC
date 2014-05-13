# compares brain clusters from two graphs of the same brain

# ---------------------------------------------------------------------
# load packages and source files
# ---------------------------------------------------------------------
if (!require(mclust)) {
    install.packages('mclust', dependencies = T)
    require(mclust)
}
if (!require(foreach)) {
    install.packages('foreach', dependencies = T)
    require(foreach)
}
if (!require(doMC)) {
    install.packages('doMC', dependencies = T)
    require(doMC)
}

source('../code/readWriteMatrix.R')

# ---------------------------------------------------------------------
# helper function to compute ari with aligned nodes
# ---------------------------------------------------------------------
ariAligned = function(clusters1, clusters2, matchInd) {
    return( adjustedRandIndex(clusters1[!is.na(matchInd)],
        na.omit(clusters2[matchInd])) )
}

# ---------------------------------------------------------------------
# load 
# ---------------------------------------------------------------------
comArgs = commandArgs(T)
preListFile = comArgs[1]
procDataDir = comArgs[2]
nCores = as.numeric(comArgs[3])

# define other directories
outDir = "cache/"
figDir = "figs/"

# read in file prefixes
preVec = readLines(preListFile)

# number of cores 
registerDoMC(nCores)

# ---------------------------------------------------------------------
# load clusters
# ---------------------------------------------------------------------
scCluster = list()
scxCluster = list()
ccaCluster = list()
cascCluster = list()

nGraphs = length(preVec)

for(i in 1:nGraphs) {
    scCluster[[i]] = as.vector(loadMatrix(paste(outDir, preVec[i], "_SC",
                 sep=""), 1))
    scxCluster[[i]] = as.vector(loadMatrix(paste(outDir, preVec[i], "_SCX",
                 sep=""), 1))
    ccaCluster[[i]] = as.vector(loadMatrix(paste(outDir, preVec[i], "_CCA",
                 sep=""), 1))
    cascCluster[[i]] = as.vector(loadMatrix(paste(outDir, preVec[i], "_CASC",
                 sep=""), 1))
}

# ---------------------------------------------------------------------
# get ari matrix for rsc, cca, casc, scx
# ---------------------------------------------------------------------
scAri = matrix(0, nrow = nGraphs, ncol = nGraphs)
scxAri = matrix(0, nrow = nGraphs, ncol = nGraphs)
ccaAri = matrix(0, nrow = nGraphs, ncol = nGraphs)
cascAri = matrix(0, nrow = nGraphs, ncol = nGraphs)

for(i in 1:(nGraphs-1)) {
    for(j in (i+1):nGraphs) {

        # load node locations for matching
        coordIF1 = paste(procDataDir, preVec[i], '_big_lcc.txt', sep='')
        coordIF2 = paste(procDataDir, preVec[j], '_big_lcc.txt', sep='')

        coorMat1 = as.matrix(read.table(coordIF1))[,2:4]
        coorMat2 = as.matrix(read.table(coordIF2))[,2:4]

        coorInd1 = coorMat1 %*% c(1, 10^3, 10^6)
        coorInd2 = coorMat2 %*% c(1, 10^3, 10^6)
    
        # match nodes
        matchInd = match(coorInd1, coorInd2)

        scAri[i,j] = ariAligned(scCluster[[i]],
                 scCluster[[j]], matchInd)
        scxAri[i,j] = ariAligned(scxCluster[[i]],
                 scxCluster[[j]], matchInd)
        ccaAri[i,j] = ariAligned(ccaCluster[[i]],
                 ccaCluster[[j]], matchInd)
        cascAri[i,j] = ariAligned(cascCluster[[i]],
                 cascCluster[[j]], matchInd)
    }
}

# write ari matricies
saveMatrixList(paste(outDir, "compareBrainClusters", sep=""),
               list(cascAri, ccaAri, scAri, scxAri))
