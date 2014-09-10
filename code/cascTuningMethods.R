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
    epsilon = .05
    
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
        orthoX[i] = cascResults$orthoX
        orthoL[i] = cascResults$orthoL
        wcssVec[i] = cascResults$wcss
        gapVec[i] = cascResults$singGap
    }

    # restrict the range of h values using orthogonal components
    subspaces = getSubspaces(orthoX, orthoL, nPoints, epsilon)
    
    subintervalLengths = subspaces$subintervalEnd - subspaces$subintervalStart
    minCountSubspaces = which(subspaces$orthoCounts ==
        min(subspaces$orthoCounts))

# min WCSS on most overlapping set of subspaces
    startIndex = subspaces$subintervalStart[minCountSubspaces]
    endIndex = subspaces$subintervalEnd[minCountSubspaces]
    minInterval = unlist(apply(cbind(startIndex, endIndex), 1, function(x)
        {x[1]:x[2]}))
    minWcssSubindex = which.min(wcssVec[minInterval])
    hOpt = (hTuningSeq[minInterval])[minWcssSubindex]
    
#    keep longest subspace if there is a tie
#    minCountLengths = subintervalLengths[minCountSubspaces]
#    maxMinSubspace = minCountSubspaces[which.max(minCountLengths)]
#    startIndex = subspaces$subintervalStart[maxMinSubspace]
#    endIndex = subspaces$subintervalEnd[maxMinSubspace]
#    minWcssIndex = which.min(wcssVec[startIndex:endIndex]) + startIndex - 1
#    hOpt = hTuningSeq[minWcssIndex]                        
        
    return( getCascSvd(graphMat, covariates, hOpt, nBlocks) )
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

# ---------------------------------------------------------------------
# Finds leading subspace discontinuities.
# Returns the start and end of a continuous interval and
# the number of orthogonal components in the leading subspace
# on the interval.
# ---------------------------------------------------------------------
getSubspaces = function(orthoX, orthoL, nPoints, epsilon) {

    indicatorOut = vector(length = nPoints)
    indicatorIn = vector(length = nPoints)
    
    for(i in 1:(nPoints - 1)) {
        if((orthoX[i] < epsilon) & (orthoX[i+1])) {
            indicatorOut[i+1] = 1
        }
        else if((orthoL[i+1] < epsilon) & (orthoL[i] > epsilon)) {
            indicatorIn[i+1] = 1
        }
    }

    orthoCounts = cumsum(indicatorIn) - cumsum(indicatorOut) +
        max(cumsum(indicatorOut))
    subintervalStart = unique(c(which(indicatorIn == 1),
        which(indicatorOut == 1)))
    subintervalEnd = sort(c(subintervalStart-1, nPoints))
    subintervalStart = sort(c(1, subintervalStart))
    orthoCounts = orthoCounts[subintervalStart]

    return( list(orthoCounts = orthoCounts,
                 subintervalStart = subintervalStart,
                 subintervalEnd = subintervalEnd) )
}
