###--------------------------------------------------------------------
# functions to simulate graph data using SBM along with covariates
###--------------------------------------------------------------------


# ---------------------------------------------------------------------
# load any necessary libraries and scripts
# ---------------------------------------------------------------------
if (!require(Matrix)) {
    install.packages('Matrix', dependencies = T)
    require(Matrix)
}


# ---------------------------------------------------------------------
# MAIN FUNCTIONS
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# function to simulate adjacency matrix SBM and Bernoulli covariates
# ---------------------------------------------------------------------
simAdjMatCovMat = function(blockMat, covProbMat, nMembers) {

    nBlocks = dim(blockMat)[1]
    nNodes = sum(nMembers)
    nCovs = dim(covProbMat)[2]
    covariates = matrix(rep(0, nCovs*nNodes), nrow = nNodes)
    adjacency = matrix(rep(0, nNodes^2), nrow = nNodes)
    upTriMat = upper.tri(adjacency)
    startBlock = cumsum(c(0, nMembers)) + 1

    #fill the aa and xx matrices block by block
    for(i in 1:nBlocks) {
	for(j in 1:i) {
		if(i == j) {
                    upTriMatTemp = upTriMat
                    upTriMatTemp[, -(startBlock[i]:(startBlock[i+1] - 1))] = F
                    upTriMatTemp[-(startBlock[i]:(startBlock[i+1] - 1)),] = F
                    adjacency[upTriMatTemp] = rbinom(sum(upTriMatTemp), 1,
                                 blockMat[i,i])
                }
		else {
                    adjacency[startBlock[i]:(startBlock[i+1] - 1),
                              startBlock[j]:(startBlock[j+1] - 1)] = 
			rbinom(nMembers[i] * nMembers[j] , 1, blockMat[i,j])
		}
	}

        #loop over covariates to sim x
	for(k in 1:nCovs) {
		covariates[startBlock[i]:(startBlock[i+1] - 1), k] =
                    rbinom(nMembers[i], 1, covProbMat[i,k])
	}
    }

    #copy upper tri to lower tri for aa
    return( list(adjacency = adjacency + t(adjacency),
                 covariates = covariates) )
}

# ---------------------------------------------------------------------
# function to simulate adjacency matrix SBM
# ---------------------------------------------------------------------
simAdjMat = function(blockMat, nMembers) {

    nBlocks = dim(blockMat)[1]
    nNodes = sum(nMembers)
    adjacency = matrix(rep(0, nNodes^2), nrow = nNodes)
    upTriMat = upper.tri(adjacency)
    startBlock = cumsum(c(0, nMembers)) + 1

    #fill the aa and xx matrices block by block
    for(i in 1:nBlocks) {
	for(j in 1:i) {
		if(i == j) {
                    upTriMatTemp = upTriMat
                    upTriMatTemp[, -(startBlock[i]:(startBlock[i+1] - 1))] = F
                    upTriMatTemp[-(startBlock[i]:(startBlock[i+1] - 1)),] = F
                    adjacency[upTriMatTemp] = rbinom(sum(upTriMatTemp), 1,
                                 blockMat[i,i])
                }
		else {
                    adjacency[startBlock[i]:(startBlock[i+1] - 1),
                              startBlock[j]:(startBlock[j+1] - 1)] = 
			rbinom(nMembers[i] * nMembers[j] , 1, blockMat[i,j])
		}
	}
    }

    #copy upper tri to lower tri for aa
    return( adjacency + t(adjacency) )
    
}

# ---------------------------------------------------------------------
# function to simulate sparse adjacency matrix SBM
# using binding instead of subsetting
# ---------------------------------------------------------------------
simSparseAdjMat = function(blockMat, nMembers) {
    nBlocks = dim(blockMat)[1]
    nNodes = sum(nMembers)
    adjacency = NULL
    
    # fill the aa and xx matrices block by block
    for(j in 1:nBlocks) {
        adjCol = NULL
	for(i in nBlocks:1) {
		if(i == j) {
                    adjTemp = triu(simBernSparseMat(nMembers[i],
                                 nMembers[i], blockMat[i,i]))    
                }
		else if(i < j) {
                    adjTemp = simBernSparseMat(nMembers[i], nMembers[j],
                                         blockMat[i, j])
		}
                else {
                    adjTemp = Matrix(0, nrow = nMembers[i],
                        ncol = nMembers[j])
                }
                adjCol = rBind(adjTemp, adjCol)
         }
        
        adjacency = cBind(adjCol, adjacency)        
    }

    # copy upper tri to lower tri for aa
    return( adjacency + t(adjacency) )   
}

# ---------------------------------------------------------------------
# function to simulate normal covariates
# ---------------------------------------------------------------------
simCoordMat = function(blockMean, blockSd, nMembers) {

    nCov = dim(blockMean)[2]
    nBlocks = length(nMembers)
    startBlock = cumsum(c(0, nMembers)) + 1
    
    covariates = matrix(0, ncol = nCov, nrow = sum(nMembers))
    for(i in 1:nBlocks) {
        for(j in 1:nCov) {
            covariates[startBlock[i]:(startBlock[i+1] - 1),j] =
                rnorm(nMembers[i], blockMean[i,j], blockSd[i,j]) 
        }
    }

    return(covariates)
}

# ---------------------------------------------------------------------
# function to simulate adjacency matrix SBM and Bernoulli covariates with
# a proportion of incorrect group assignments
# ---------------------------------------------------------------------
simAdjMatCovMatPropInc = function(blockMat, covProbMat, nMembers,
    propIncorrect) {

    nBlocks = dim(blockMat)[1]
    nNodes = sum(nMembers)
    nCovs = dim(covProbMat)[2]
    covariates = matrix(rep(0, nCovs*nNodes), nrow = nNodes)
    adjacency = matrix(rep(0, nNodes^2), nrow = nNodes)
    upTriMat = upper.tri(adjacency)
    startBlock = cumsum(c(0, nMembers)) + 1

    #fill the aa and xx matrices block by block
    for(i in 1:nBlocks) {
	for(j in 1:i) {
		if(i == j) {
                    upTriMatTemp = upTriMat
                    upTriMatTemp[, -(startBlock[i]:(startBlock[i+1] - 1))] = F
                    upTriMatTemp[-(startBlock[i]:(startBlock[i+1] - 1)),] = F
                    adjacency[upTriMatTemp] = rbinom(sum(upTriMatTemp), 1,
                                 blockMat[i,i])
                }
		else {
                    adjacency[startBlock[i]:(startBlock[i+1] - 1),
                              startBlock[j]:(startBlock[j+1] - 1)] = 
			rbinom(nMembers[i] * nMembers[j] , 1, blockMat[i,j])
		}
	}

        #loop over covariates to sim x
	for(k in 1:nCovs) {
		covariates[startBlock[i]:(startBlock[i+1] - 1), k] =
                    rbinom(nMembers[i], 1, covProbMat[i,k])
	}
    }

    #flip the group membership of the given number of prop_inc 
    #reassign to the next group
    nFlips = round( sum(nMembers) * propIncorrect / length(nMembers) )
    covariatesTemp = covariates[startBlock[1]:(startBlock[1] + nFlips), ]

    for(i in 1:(nBlocks - 1)) {
        covariates[startBlock[i]:(startBlock[i] + nFlips), ] =
            covariates[startBlock[i+1]:(startBlock[i+1] + nFlips), ]
    }

    covariates[startBlock[i+1]:(startBlock[i+1] + nFlips), ] = covariatesTemp
    
    #copy upper tri to lower tri for aa
    return( list(adjacency = adjacency + t(adjacency),
                 covariates = covariates) )
}

# ---------------------------------------------------------------------
# HELPER FUNCTIONS
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# function to generate block matrix for planted partition model
# ---------------------------------------------------------------------
genBlockMatPPM = function(p, q, nBlocks) {
    
    return( matrix( c(rep( c(p, rep(q, nBlocks)), nBlocks - 1 ), p),
                   nrow = nBlocks )  )
}

# ---------------------------------------------------------------------
# function to generate block matrix for SBM with same off diagnols
# ---------------------------------------------------------------------
genBlockMatSBM = function(pVec, q, nBlocks) {

  bVec = rep(0, nBlocks^2)
  indx = 1

  for(i in 1:(nBlocks - 1)) {
      bVec[indx] = pVec[i]
      bVec[(indx + 1):(indx + nBlocks)] = q
      indx = indx + nBlocks + 1
  }

  bVec[indx] = pVec[nBlocks]
  
  return( matrix(bVec, nrow = nBlocks) )
}

# ---------------------------------------------------------------------
# function to generate sparse Bernoulli matrix
# ---------------------------------------------------------------------
simBernSparseMat = function(nRow, nCol, p) {

    if(p == 0) {
        return( Matrix(0, nrow = nRow, ncol = nCol) )
    }

    expNumOnes = nRow*nCol*p
    sdNumOnes = sqrt(nRow*nCol*p*(1-p))

    indOnes = rnbinom(expNumOnes + round(4*sdNumOnes), 1, p) + 1

    indCumSum = cumsum(indOnes)
    
    # in the highly unlikely case of early termination
    while (max(indCumSum) < nRow*nCol) {
        indCumSum = c(indCumSum, max(indCumSum) +
            rnbinom(1, 1, p) + 1)
    }
    
    indCumSum = indCumSum[indCumSum <= nRow*nCol]

    return( Matrix(sparseVector(1, indCumSum, nRow*nCol) ,
                   nrow = nRow, ncol = nCol) )
}
