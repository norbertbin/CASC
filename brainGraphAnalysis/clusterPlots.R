# ---------------------------------------------------------------------
# load packages and source files
# ---------------------------------------------------------------------
if (!require(colorspace)) {
    install.packages('colorspace', dependencies = T)
    require(colorspace)
}

if (!require(lattice)) {
    install.packages('lattice', dependencies = T)
    require(lattice)
}

source("../code/readWriteMatrix.R")
source("../code/helperFunctions.R")

# ---------------------------------------------------------------------
# define directories
# ---------------------------------------------------------------------
comArgs = commandArgs(T)
procDataDir = comArgs[1]
filePre = comArgs[2]

# input file with coordinate location of nodes
lccCovInputFile = paste(procDataDir, filePre, '_big_lcc.txt', sep='')
print(lccCovInputFile)

# input file with cluster assignments
scInputFile = paste(outDir, filePre, "_SC", sep="")
cascInputFile = paste(outDir, filePre, "_CASC", sep="")
ccaInputFile = paste(outDir, filePre, "_CCA", sep="")
scxInputFile = paste(outDir, filePre, "_SCX", sep="")

# output directory
figDir = "figs/"

# ---------------------------------------------------------------------
# load data and clusters
# ---------------------------------------------------------------------
# read file with largest connected component and coordinates
covData = read.table(paste(procDataDir, filePre, '_proc_coord.txt', sep=''),
    header = T)

nNodes = length(covData$V1)
nLayers = 20

#load the cluster assignments
scCluster = loadMatrix(scInputFile, 1)
cascCluster = loadMatrix(cascInputFile, 1)
ccaCluster = loadMatrix(ccaInputFile, 1)
scxCluster = loadMatrix(scxInputFile, 1)

# ---------------------------------------------------------------------
# make z,y,z cross sections
# ---------------------------------------------------------------------
# get cutoff ranges along z axis
zLayer = getLayer(nLayers, covData[,3], nNodes)

# get cutoff ranges along x axis
xLayer = getLayer(nLayers, covData[,1], nNodes)

# get cutoff ranges along y axis
yLayer = getLayer(nLayers, covData[,2], nNodes)

# ---------------------------------------------------------------------
# plot sc and casc clusters
# ---------------------------------------------------------------------
# create data frame for ploting
nodeData = data.frame(x = covData[,1], y = covData[,2], z = covData[,3],
    cascCluster, scCluster, ccaCluster, scxCluster, xLayer, yLayer, zLayer)

# generate good set of colors in good order
firstCol = c(70, 30, 90, 10, 60, 40, 1, 100, 80, 20)
colOrder = c(firstCol, (1:100)[-firstCol])
sortedClusterIndex = sort.int(table(cascCluster), decreasing = T,
    index.return = T)$ix
colOrder[sortedClusterIndex] = colOrder
colSeq = rainbow(100)[colOrder]

# generate a plot for each layer
for(i in 2:(nLayers-1)) {

    #CASC results
    outFileName = paste(figDir, 'cascZLayer_', i, '.png', sep='')
    plotLayer(cascCluster[zLayer == i], nodeData[zLayer == i,],
              'z', colSeq, outFileName)

    outFileName = paste(figDir, 'cascYLayer_', i, '.png', sep='')
    plotLayer(cascCluster[yLayer == i], nodeData[yLayer == i,],
              'y', colSeq, outFileName)

    outFileName = paste(figDir, 'cascXLayer_', i, '.png', sep='')
    plotLayer(cascCluster[xLayer == i], nodeData[xLayer == i,],
              'x', colSeq, outFileName)

    #spectral clustering results
    outFileName = paste(figDir, 'scZLayer_', i, '.png', sep='')
    plotLayer(scCluster[zLayer == i], nodeData[zLayer == i,],
              'z', colSeq, outFileName)
    
    outFileName = paste(figDir, 'scYLayer_', i, '.png', sep='')
    plotLayer(scCluster[yLayer == i], nodeData[yLayer == i,],
              'y', colSeq, outFileName)

    outFileName = paste(figDir, 'scXLayer_', i, '.png', sep='')
    plotLayer(scCluster[xLayer == i], nodeData[xLayer == i,],
              'x', colSeq, outFileName)

    #cca clustering results
    outFileName = paste(figDir, 'ccaZLayer_', i, '.png', sep='')
    plotLayer(ccaCluster[zLayer == i], nodeData[zLayer == i,],
              'z', colSeq, outFileName)
    
    outFileName = paste(figDir, 'ccaYLayer_', i, '.png', sep='')
    plotLayer(ccaCluster[yLayer == i], nodeData[yLayer == i,],
              'y', colSeq, outFileName)

    outFileName = paste(figDir, 'ccaXLayer_', i, '.png', sep='')
    plotLayer(ccaCluster[xLayer == i], nodeData[xLayer == i,],
              'x', colSeq, outFileName)

    #sc on x clustering results
    outFileName = paste(figDir, 'scxZLayer_', i, '.png', sep='')
    plotLayer(scxCluster[zLayer == i], nodeData[zLayer == i,],
              'z', colSeq, outFileName)
    
    outFileName = paste(figDir, 'scxYLayer_', i, '.png', sep='')
    plotLayer(scxCluster[yLayer == i], nodeData[yLayer == i,],
              'y', colSeq, outFileName)

    outFileName = paste(figDir, 'scxXLayer_', i, '.png', sep='')
    plotLayer(scxCluster[xLayer == i], nodeData[xLayer == i,],
              'x', colSeq, outFileName)
}

# ---------------------------------------------------------------------
# plot large clusters splitting up
# ---------------------------------------------------------------------
# identify the largest sc cluster
maxScClust = which.max(table(scCluster))
maxClustNodes = which(scCluster == maxScClust)
numCascClust = length(table(cascCluster[maxClustNodes]))

# generate color seq
colSeq = rainbow(numCascClust)

# get cluster dimensions
xRange = range(covData[maxClustNodes,1])
yRange = range(covData[maxClustNodes,2])
zRange = range(covData[maxClustNodes,3])

# since x and y have the largest range project onto them
png(paste(figDir, "scMaxCluster", ".png", sep=""))
print(
    xyplot(y ~ x, group = scCluster, data = nodeData[maxClustNodes,],
           pch = 20, col = 'blue',
           main = paste("z = ", round(zRange[1], 3), " to ",
               round(xRange[2], 3)))
    )
dev.off()

png(paste(figDir, "cascMaxCluster", ".png", sep=""))
print(
    xyplot(y ~ x, group = cascCluster, data = nodeData[maxClustNodes,],
           pch = 20, col = colSeq,
           main = paste("z = ", round(zRange[1], 3), " to ",
               round(zRange[2], 3)))
    )
dev.off()

# ---------------------------------------------------------------------
# plot small clusters combining
# ---------------------------------------------------------------------
# identify the small sc clusters
smallScClust = which(table(scCluster) < 200)
smallClustNodes = which(scCluster %in% smallScClust)
numSmallScClust = length(smallScClust)

# get the casc clusters to which these nodes belong
smallCascClust = unique(cascCluster[smallClustNodes])
# get all the nodes in these clusters
cascSmallClustNodes = which(cascCluster %in% smallCascClust)
numSmallCascClust = length(smallCascClust)

# generate color seq
colSeqSc = rainbow(numSmallScClust)
colSeqCasc = rainbow(numSmallCascClust)
    
# get cluster dimensions
xRange = range(covData[smallClustNodes,1])
yRange = range(covData[smallClustNodes,2])
zRange = range(covData[smallClustNodes,3])

# since x and y have the largest range project onto them
png(paste(figDir, "scSmallCluster1", ".png", sep=""))
print(
    xyplot(y ~ x, group = scCluster, data = nodeData[smallClustNodes,],
           pch = 20, col = colSeqSc,
           main = paste("z = ", round(zRange[1], 3), " to ",
               round(xRange[2], 3)))
    )
dev.off()

png(paste(figDir, "scSmallCluster2", ".png", sep=""))
print(
    xyplot(y ~ x, group = cascCluster, data = nodeData[smallClustNodes,],
           pch = 20, col = colSeqCasc,
           main = paste("z = ", round(zRange[1], 3), " to ",
               round(xRange[2], 3)))
    )
dev.off()

png(paste(figDir, "cascSmallCluster", ".png", sep=""))
print(
    xyplot(y ~ x, group = cascCluster, data = nodeData[cascSmallClustNodes,],
           pch = 20, col = colSeqCasc,
           main = paste("z = ", round(zRange[1], 3), " to ",
               round(zRange[2], 3)))
    )
dev.off()
