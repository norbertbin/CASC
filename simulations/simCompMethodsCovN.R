# ---------------------------------------------------------------------
# simulations to compare CCA, CASC, SC on graph or covariates only
# versus changes in the signal strength in the covariates
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
# simulations to compare CCA, CASC, SC on graph or covariates only
# versus different numbers of covariates
# ---------------------------------------------------------------------
#model parameters
nBlocks = 2
nMembers = rep(500, nBlocks)
nCovs = 2
covProb1 = .5
covProb2 = .1
p = .03 
q = .016
nCovSeq = c(2, 4, 6, 8, 10)
nIter = 20
method = "regLaplacian"

nPoints = length(nCovSeq)
misClustRateCasc = rep(0, nPoints)
misClustRateCca = rep(0, nPoints)
misClustRateL = rep(0, nPoints)
misClustRateX = rep(0, nPoints)

for(iter in 1:nIter) {
    index = 0
    for(nCov in nCovSeq) {

        #simulate the graph and covariates
        covProbMat = matrix(c( rep(c(covProb1, covProb2), nCov/2),
            rep(c(covProb2, covProb1), nCov/2) ), nrow = nBlocks)
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
misClustData = data.frame(nCov = rep(nCovSeq, 4), misClustRate =
    c(misClustRateCasc, misClustRateCca, misClustRateL, misClustRateX),
    group = rep(1:4, each = nPoints))

#output txt file with results
write.table(misClustData, file = paste(outDir, "simCompMethodsNumCov.txt",
                              sep = ""))

#create figure
pdf(paste(outDir, "simCompMethodsNumCov.pdf", sep = ""),
    width = 7, height = 7)
print(
    xyplot(misClustRate ~ nCov, group = group, type = "b", pch = 1:4,
           data = misClustData, ylab = "Average Mis-clustering Rate",
           xlab = "Number of Covariates", lwd = 2, key = list(
                               text = list(c("CASC","CCA","SC-L","SC-X")),
                               lines = list(col = 
			trellis.par.get()$superpose.symbol$col[1:4]),
                               points = list(pch = 1:4, 
			col = trellis.par.get()$superpose.symbol$col[1:4]),
                                        corner = c(.95, .95)) )
    )
dev.off()
