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
    
    singValGraph = irlba(graphMatrix, nu = n, nv = 0, m_b =
        internalDim)$d
    maxSingValCov = svd(covariates, nu = 1)$d[1]

    hmax = singValGraph[nBlocks]/maxSingValCov

    hmin = (singValGraph[nBlocks - 1] - singValGraph[nBlocks])/maxSingValCov
    
    return( list( hmax = hmax, hmin = hmin ) )
}
