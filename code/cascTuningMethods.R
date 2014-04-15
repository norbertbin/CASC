# ---------------------------------------------------------------------
# functions to run CASC with parameter tuning
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# load packages and scripts
# ---------------------------------------------------------------------
source("spectralClusteringMethods.R")


# ---------------------------------------------------------------------
# MAIN FUNCTIONS
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# returns CASC block memberships
# ---------------------------------------------------------------------
getCascAutoClusters = function(adjacency, covariates, nBlocks,
    method = "regLaplacian", nPoints = 100) {

    graphMat = getGraphMatrix(adjacency)
    rangehTuning = getTuningRange(graphMat, covariates, nBlocks)

    hTuningSeq = seq(rangehTuning[1], rangehTuning[2], nPoints)
    wcssVec = rep(0, nPoints)
    clusterMat = matrix(rep(0, nPoints*dim(graphMat)[1]), nrow = nPoints)

    for(i in 1:nPoints) {
        cascResults = getCascResults(graphMat, covariates, hTuningSeq[i],
            nBlocks)
        wcssVec[i] = cascResults$wcss
        clusterMat[i, ] = cascResults$cluster
    }

    minWcssIndex = match(min(wcssVec), wcssVec)

    return(clusterMat[minWcssIndex, ])
}

# ---------------------------------------------------------------------
# returns CASC optimal h tuning parameter SVD
# ---------------------------------------------------------------------
getCascAutoSvd = function(graphMat, covariates, nBlocks,
    nPoints = 100) {

    rangehTuning = getTuningRange(graphMat, covariates, nBlocks)

    hTuningSeq = seq(rangehTuning$hmin, rangehTuning$hmax,
        length.out = nPoints)
    wcssVec = rep(0, nPoints)
    gapVec = rep(0, nPoints)
    
    for(i in 1:nPoints) {
        cascResults = getCascResults(graphMat, covariates, hTuningSeq[i],
            nBlocks)
        wcssVec[i] = cascResults$wcss
        gapVec[i] = cascResults$singGap
    }

    # restrict possible values to those past any phase transition
    if(min(gapVec) < .9*min(gapVec[1], gapVec[nPoints]) &
       min(gapVec) < .02) {
        starth = match(min(gapVec), gapVec) + 1
        print("transition found")
    }
    else {
        starth = 1
    }
        
    minWcssIndex = match(min(wcssVec[starth:nPoints]),
        wcssVec[starth:nPoints]) + starth - 1
print(paste("hTuning = ", hTuningSeq[minWcssIndex]))
    return( getCascSvd(graphMat, covariates, hTuningSeq[minWcssIndex],
                       nBlocks) )
}

# ---------------------------------------------------------------------
# gets a good range for the tuning parameter in CASC
# ---------------------------------------------------------------------
getTuningRange = function(graphMatrix, covariates, nBlocks) {
    
    #insure irlba internal representation is large enough
    if(nBlocks > 10) {
        internalDim = 2 * nBlocks
    }
    else {
        internalDim = 20
    }
    
    singValGraph = irlba(graphMatrix, nu = nBlocks + 1, nv = 0, m_b =
        internalDim)$d
    singValCov = svd(covariates, nu = nBlocks)$d

    hmax = singValGraph[1]/singValCov[nBlocks]^2

    hmin = (singValGraph[nBlocks] - singValGraph[nBlocks + 1])/singValCov[1]^2

    print(paste("hmax = ", hmax, " hmin = ", hmin))
    return( list( hmax = hmax, hmin = hmin ) )
}
