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



# ---------------------------------------------------------------------
# simulations to compare CCA, CASC, SC on graph or covariates only
# versus changes in the signal strength in the covariates
# ---------------------------------------------------------------------
# tuning procedure
enhancedTuning = F

# center and scale X
centerAndScale = T

outDir = paste("Figs/", "enhanced", enhancedTuning, "centered",
    centerAndScale, "/", sep="")

#model parameters
nBlocks = 2
nMembers = rep(500, nBlocks)
nCovs = 2
covProb1 = .5
p = .03 
q = .016
deltaCovProbSeq = seq(0, .4, by = .05)
nIter = 100
method = "regLaplacian"

nPoints = length(deltaCovProbSeq)
misClustRateCasc = rep(0, nPoints)
misClustRateCca = rep(0, nPoints)
misClustRateL = rep(0, nPoints)
misClustRateX = rep(0, nPoints)

for(iter in 1:nIter) {
    index = 0
    for(deltaCovProb in deltaCovProbSeq) {

        #simulate the graph and covariates
        covProb2 = covProb1 - deltaCovProb
        covProbMat = matrix(c(covProb1, covProb2, covProb2, covProb1)
            , nrow = nCovs)
        blockProbMat = genBlockMatPPM(p, q, nBlocks)
        simData = simAdjMatCovMat(blockProbMat, covProbMat, nMembers)

        if(centerAndScale == T) {
            simData$covariates = scale(simData$covariates, center = T,
                           scale = sqrt(colSums(simData$covariates^2)))
            # add a constant column
            nNodes = sum(nMembers)
            const = rep(1/sqrt(nNodes), nNodes)            
            simData$covariates = cBind(simData$covariates, const)
        }
        
        #convert adjacency matrix to sparse matrix
        graphMat = getGraphMatrix( Matrix(simData$adjacency, sparse = T),
            method)

        #calculate the misclustering rate for this and add to total
        index = index + 1
        misClustRateCasc[index] = misClustRateCasc[index] +
            misClustRateEmp(getCascAutoSvd(graphMat, simData$covariates,
               nBlocks, enhancedTuning = enhancedTuning)$singVec, nMembers)
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
misClustData = data.frame(deltaCovProb = rep(deltaCovProbSeq, 4), misClustRate =
    c(misClustRateCasc, misClustRateCca, misClustRateL, misClustRateX),
    group = rep(1:4, each = nPoints))

#output txt file with results
write.table(misClustData, file = paste(outDir, "simCompMethodsCovProb.txt",
                              sep = ""))

#create figure
pdf(paste(outDir, "simCompMethodsCovProb.pdf", sep = ""),
    width = 7, height = 7)
print(
    xyplot(misClustRate ~ deltaCovProb, group = group, type = "b", pch = 1:4,
           data = misClustData, ylab = "Average Mis-clustering Rate",
           xlab = "Difference in Covariate Probabilities", lwd = 2, key = list(
                               text = list(c("CASC","CCA","SC-L","SC-X")),
                               lines = list(col = 
			trellis.par.get()$superpose.symbol$col[1:4]),
                               points = list(pch = 1:4, 
			col = trellis.par.get()$superpose.symbol$col[1:4]),
                                        corner = c(.95, .95)) )
    )
dev.off()
