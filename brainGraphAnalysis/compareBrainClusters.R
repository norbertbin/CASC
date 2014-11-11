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
if (!require(gplots)) {
    install.packages('gplots', dependencies = T)
    require(gplots)
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

if(!file.exists(paste(outDir, "compareBrainClusters.bin", sep=""))) {
# ---------------------------------------------------------------------
# load clusters
# ---------------------------------------------------------------------
scCluster = list()
scxCluster = list()
ccaCluster = list()
cascCluster = list()
baCluster = list()

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
    baCluster[[i]] = read.table(paste(outDir, preVec[i],
                 "_big_graph_w_inv_attr.txt", sep=""))$x
}

# ---------------------------------------------------------------------
# get ari matrix for rsc, cca, casc, scx
# ---------------------------------------------------------------------
scAri = matrix(0, nrow = nGraphs, ncol = nGraphs)
scxAri = matrix(0, nrow = nGraphs, ncol = nGraphs)
ccaAri = matrix(0, nrow = nGraphs, ncol = nGraphs)
cascAri = matrix(0, nrow = nGraphs, ncol = nGraphs)
baAri = matrix(0, nrow = nGraphs, ncol = nGraphs)

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
        baAri[i,j] = ariAligned(baCluster[[i]],
                 baCluster[[j]], matchInd)
    }
}

# write ari matricies
saveMatrixList(paste(outDir, "compareBrainClusters", sep=""),
               list(cascAri, ccaAri, scAri, scxAri, baAri))
} else {
	cascAri = loadMatrix(paste(outDir, "compareBrainClusters", sep=""), 1)
	ccaAri = loadMatrix(paste(outDir, "compareBrainClusters", sep=""), 2)
	scAri = loadMatrix(paste(outDir, "compareBrainClusters", sep=""), 3)
	scxAri = loadMatrix(paste(outDir, "compareBrainClusters", sep=""), 4)
    baAri = loadMatrix(paste(outDir, "compareBrainClusters", sep=""), 5)
}

ordLabels = c(1, 25, 2, 37, 3, 22, 4, 11, 5, 31, 6, 20, 7, 34, 8, 29, 9, 42, 10, 21, 12, 19, 13, 24, 14, 17, 15, 26, 16, 35, 18, 38, 23, 27, 28, 40, 30, 33, 32, 36, 39, 41)

# draw ari matrix
pdf(paste(figDir, "compareScAri.pdf", sep=""), width = 10, height = 10)
heatmap.2(scAri, Rowv = NULL, Colv = NULL, dendrogram = "none", key = F, trace = "none", notecol="dodgerblue2", col = colorRampPalette(c("white", "black"))(256), lwid = c(.09, .99), lhei = c(.01, .99), margins = c(5, 5), cellnote = round(scAri, 2), labRow = ordLabels, labCol = ordLabels)
dev.off()

pdf(paste(figDir, "compareCcaAri.pdf", sep=""), width = 10, height = 10)
heatmap.2(ccaAri, Rowv = NULL, Colv = NULL, dendrogram = "none", key = F, trace = "none", notecol="dodgerblue2", col = colorRampPalette(c("white", "black"))(256), lwid = c(.09, .99), lhei = c(.01, .99), margins = c(5, 5), cellnote = round(ccaAri, 2), labRow = ordLabels, labCol = ordLabels)
dev.off()

pdf(paste(figDir, "compareCascAri.pdf", sep=""), width = 10, height = 10)
heatmap.2(cascAri, Rowv = NULL, Colv = NULL, dendrogram = "none", key = F, trace = "none", notecol="dodgerblue2", col = colorRampPalette(c("white", "black"))(256), lwid = c(.09, .99), lhei = c(.01, .99), margins = c(5, 5), cellnote = round(cascAri, 2), labRow = ordLabels, labCol = ordLabels)
dev.off()

# do hierarchical clustering using ari
pdf(paste(figDir, "cascAriHier.pdf", sep=""), width = 10, heigh = 10)
plot(hclust(as.dist(t(1-cascAri))), xlab = "", main = "Cluster Dendrogram Using
ARI for CASC", labels = ordLabels)
dev.off()

pdf(paste(figDir, "rscAriHier.pdf", sep=""), width = 10, heigh = 10)
plot(hclust(as.dist(t(1-scAri))), xlab = "", main = "Cluster Dendrogram Using
ARI for RSC", labels = ordLabels)
dev.off()

pdf(paste(figDir, "ccaAriHier.pdf", sep=""), width = 10, heigh = 10)
plot(hclust(as.dist(t(1-ccaAri))), xlab = "", main = "Cluster Dendrogram Using
ARI for CCA", labels = ordLabels)
dev.off()
