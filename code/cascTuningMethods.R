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

    # value for detecting a transition
    epsilon = .01
    
    rangehTuning = getTuningRange(graphMat, covariates, nBlocks)

    hTuningSeq = seq(rangehTuning$hmin, rangehTuning$hmax,
        length.out = nPoints)
    wcssVec = vector(length = nPoints)
    gapVec = vector(length = nPoints)
    orthoX = vector(length = nPoints)
    orthoL = vector(length = nPoints)
    
    for(i in 1:nPoints) {
        cascResults = getCascResults(graphMat, covariates, hTuningSeq[i],
            nBlocks)
        orthoX[i] = sum((t(cascResults$singVecK)%*%covariates)^2)
        orthoL[i] = sum((t(cascResults$singVecK)%*%graphMat)^2)
        wcssVec[i] = cascResults$wcss
        gapVec[i] = cascResults$singGap
    }

    # restrict the range of h values
    startIndex = 1
    endIndex = nPoints
    for(i in 1:(nPoints-1)) {
        if(orthoX[i] < epsilon & orthoX[i+1] > epsilon &
           orthoL[i+1] > epsilon) {
            startIndex = i + 1
        }
        if(orthoL[i+1] < epsilon & orthoX[i] > epsilon &
           orthoL[i] > epsilon) {
            endIndex = i
        }
    }
            
    minWcssIndex = which.min(wcssVec[startIndex:endIndex]) + startIndex - 1

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

    return( list( hmax = hmax, hmin = hmin ) )
}
