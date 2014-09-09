# ---------------------------------------------------------------------
# simulations to compare CCA, CASC, SC on graph or covariates only
# versus changes in the signal strength in the graph
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# load required scripts and packages
# ---------------------------------------------------------------------
initialDir = getwd()
setwd("../code/")
source("spectralClusteringMethods.R")
source("simulateGraphData.R")
source("cascTuningMethods.R")
source("extEvalMethods.R")
setwd(initialDir)

if (!require(lattice)) {
    install.packages('lattice', dependencies = T)
    require(lattice)
}

outDir = "Figs/"

# ---------------------------------------------------------------------
# comparison of misclustering using CCA, CASC, SC on L, and SC on X
# versus p-q for fixed covProb1 - covProb2
# ---------------------------------------------------------------------
#model parameters
nBlocks = 2
nMembers = rep(500, nBlocks)
nCovs = 2
covProb1 = .5
covProb2 = .1
p = .03 #might need to adjust this
deltaPQSeq = seq(0, .02, by = .002)
nIter = 100
method = "regLaplacian"

nPoints = length(deltaPQSeq)
misClustRateCasc = rep(0, nPoints)
misClustRateCca = rep(0, nPoints)
misClustRateL = rep(0, nPoints)
misClustRateX = rep(0, nPoints)

for(iter in 1:nIter) {
    index = 0
    for(deltaPQ in deltaPQSeq) {

        #simulate the graph and covariates
        q = p - deltaPQ
        covProbMat = matrix(c(covProb1, covProb2, covProb2, covProb1)
            , nrow = nCovs)
        blockProbMat = genBlockMatPPM(p, q, nBlocks)
        simData = simAdjMatCovMat(blockProbMat, covProbMat, nMembers)

        #convert adjacency matrix to sparse matrix
        graphMat = getGraphMatrix( Matrix(simData$adjacency, sparse = T),
            method)

        #calculate the misclustering rate for this and add to total
        index = index + 1
        misClustRateCasc[index] = misClustRateCasc[index] +
            misClustRateEmp(getCascAutoSvd(graphMat, simData$covariates,
                                           nBlocks)$singVec, nMembers)
        misClustRateCca[index] = misClustRateCca[index] +
            misClustRateEmp(getCcaSvd(graphMat, simData$covariates,
                                      nBlocks)$singVec, nMembers)
        misClustRateL[index] = misClustRateL[index] +
            misClustRateEmp(getGraphScSvd(graphMat, nBlocks)$singVec,
                            nMembers)
        misClustRateX[index] = misClustRateX[index] +
            misClustRateEmp(getCovScSvd(simData$covariates,
                                        nBlocks)$singVec, nMembers)
                      
    }
}

misClustRateCasc = misClustRateCasc/nIter
misClustRateCca = misClustRateCca/nIter
misClustRateL = misClustRateL/nIter
misClustRateX = misClustRateX/nIter

#create data frame with results
misClustData = data.frame(deltaPQ = rep(deltaPQSeq, 4), misClustRate =
    c(misClustRateCasc, misClustRateCca, misClustRateL, misClustRateX),
    group = rep(1:4, each = nPoints))

#output txt file with results
write.table(misClustData, file = paste(outDir, "simCompMethodsPQ.txt",
                              sep = ""))

#create figure
pdf(paste(outDir, "simCompMethodsPQ.pdf", sep = ""),
    width = 7, height = 7)
print(
xyplot(misClustRate ~ deltaPQ, group = group, type = "b", pch = 1:4,
           data = misClustData, ylab = "Average Mis-clustering Rate",
           xlab = "p - q", lwd = 2, key = list(
                               text = list(c("CASC","CCA","SC-L","SC-X")),
                               lines = list(col = 
			trellis.par.get()$superpose.symbol$col[1:4]),
                               points = list(pch = 1:4, 
			col = trellis.par.get()$superpose.symbol$col[1:4]),
                                        corner = c(.95, .95)) )
)
dev.off()

# ---------------------------------------------------------------------
# comparison of misclustering using CASC, CCA, L, X versus n for
# increasing density
# ---------------------------------------------------------------------
#model parameters
nBlocks = 2
nCovs = 2
covProb1 = .5
covProb2 = .1
p = .03
q = .021
nNodesSeq = 500*c(1, 2, 4, 6, 8, 10)
nIter = 100
method = "regLaplacian"

nPoints = length(nNodesSeq)
misClustRateCasc = rep(0, nPoints)
misClustRateCca = rep(0, nPoints)
misClustRateL = rep(0, nPoints)
misClustRateX = rep(0, nPoints)

for(iter in 1:nIter) {
    index = 0
    for(nNodes in nNodesSeq) {

        nMembers = c(nNodes/2, nNodes/2)

        #simulate the graph and covariates
        covProbMat = matrix(c(covProb1, covProb2, covProb2, covProb1)
            , nrow = nCovs)
        blockProbMat = genBlockMatPPM(p, q, nBlocks)
        simData = simAdjMatCovMat(blockProbMat, covProbMat, nMembers)

        #convert adjacency matrix to sparse matrix
        graphMat = getGraphMatrix( Matrix(simData$adjacency, sparse = T),
            method)

        #calculate the misclustering rate for this and add to total
        index = index + 1
        misClustRateCasc[index] = misClustRateCasc[index] +
            misClustRateEmp(getCascAutoSvd(graphMat, simData$covariates,
                                           nBlocks)$singVec, nMembers)
        misClustRateCca[index] = misClustRateCca[index] +
            misClustRateEmp(getCcaSvd(graphMat, simData$covariates,
                                      nBlocks)$singVec, nMembers)
        misClustRateL[index] = misClustRateL[index] +
            misClustRateEmp(getGraphScSvd(graphMat, nBlocks)$singVec,
                            nMembers)
        misClustRateX[index] = misClustRateX[index] +
            misClustRateEmp(getCovScSvd(simData$covariates,
                                        nBlocks)$singVec, nMembers)
                      
    }
}

misClustRateCasc = misClustRateCasc/nIter
misClustRateCca = misClustRateCca/nIter
misClustRateL = misClustRateL/nIter
misClustRateX = misClustRateX/nIter

#create data frame with results
misClustData = data.frame(nNodes = rep(nNodesSeq, 4), misClustRate =
    c(misClustRateCasc, misClustRateCca, misClustRateL, misClustRateX),
    group = rep(1:4, each = nPoints))

#output txt file with results
write.table(misClustData, file = paste(outDir, "simCompMethodsDense.txt",
                              sep = ""))

#create figure
pdf(paste(outDir, "simCompMethodsDense.pdf", sep = ""), width = 7, height = 7)
print(
    xyplot(misClustRate ~ nNodes, group = group, type = "b", pch = 1:4,
           data = misClustData, ylab = "Average Mis-clustering Rate",
           xlab = "Number of Nodes", lwd = 2, key = list(
                               text = list(c("CASC","CCA","SC-L","SC-X")),
                               lines = list(col = 
			trellis.par.get()$superpose.symbol$col[1:4]),
                               points = list(pch = 1:4, 
			col = trellis.par.get()$superpose.symbol$col[1:4]),
                                        corner = c(.95, .95)) )
    )
dev.off()
