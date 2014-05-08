# compares brain clusters from two graphs of the same brain

# ---------------------------------------------------------------------
# load packages and source files
# ---------------------------------------------------------------------
if (!require(mclust)) {
    install.packages('mclust', dependencies = T)
    require(mclust)
}

source('../code/readWriteMatrix.R')

# ---------------------------------------------------------------------
# load 
# ---------------------------------------------------------------------
comArgs = commandArgs(T)
preListFile = comArgs[1]
procDataDir = comArgs[2]

# define other directories
outDir = "cache/"
figDir = "figs/"

# read in file prefixes
preVec = readLines(preListFile)

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

# load pairing information as vector with person id
pairVec = rep(1:12, each = 2)

# ---------------------------------------------------------------------
# get ari matrix for rsc, cca, casc, scx
# ---------------------------------------------------------------------
scAri = rep(0, nGraphs/2)
scxAri = rep(0, nGraphs/2)
ccaAri = rep(0, nGraphs/2)
cascAri = rep(0, nGraphs/2)

for(i in unique(pairVec)) {

    # get ith pair of graphs
    pairInd = which(i == pairVec) 

    # load node locations for matching
    coordIF1 = paste(procDataDir, preVec[pairInd[1]], '_big_lcc.txt', sep='')
    coordIF2 = paste(procDataDir, preVec[pairInd[2]], '_big_lcc.txt', sep='')

    coorMat1 = as.matrix(read.table(coordIF1))[,2:4]
    coorMat2 = as.matrix(read.table(coordIF2))[,2:4]

    coorInd1 = coorMat1 %*% c(1, 10^3, 10^6)
    coorInd2 = coorMat2 %*% c(1, 10^3, 10^6)
    
    # match nodes
    matchInd = match(coorInd1, coorInd2)

    scAri[i] = ariAligned(scCluster[[pairInd[1]]], scCluster[[pairInd[2]]],
        matchInd)
    scxAri[i] = ariAligned(scxCluster[[pairInd[1]]],
              scxCluster[[pairInd[2]]], matchInd)
    ccaAri[i] = ariAligned(ccaCluster[[pairInd[1]]],
              ccaCluster[[pairInd[2]]], matchInd)
    cascAri[i] = ariAligned(cascCluster[[pairInd[1]]],
               cascCluster[[pairInd[2]]], matchInd)
}

# ---------------------------------------------------------------------
# helper function to compute ari with aligned nodes
# ---------------------------------------------------------------------
ariAligned = function(clusters1, clusters2, matchInd) {
    return( adjustedRandIndex(clusters1[!is.na(matchInd)],
        na.omit(clusters2[matchInd])) )
}
