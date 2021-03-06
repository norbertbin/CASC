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
    startBlock = cumsum(c(0, nMembers)) + 1

    #fill covariate matrices block by block
    for(i in 1:nBlocks) {
         #loop over covariates to sim x
        for(k in 1:nCovs) {
            covariates[startBlock[i]:(startBlock[i+1] - 1), k] =
                    rbinom(nMembers[i], 1, covProbMat[i,k])
        }
    }

    #copy upper tri to lower tri for aa
    return( list(adjacency = simSparseAdjMat(blockMat, nMembers),
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
simSparseAdjMat = function(bMat, nMembers) {
 nBlocks = dim(bMat)[1]
 nNodes = sum(nMembers)
 adjM = NULL
    
 for(j in nBlocks:1) {
  adjCol = NULL
  for(i in nBlocks:1) {
  	if(i <= j) {
      adjTemp = Matrix( 
       simBernSparseVec(nMembers[i],
             nMembers[j], bMat[i,j]),
        nrow = nMembers[i], ncol = nMembers[j])
    }  
    else {
     adjTemp = Matrix(0, nrow = nMembers[i],
                         ncol = nMembers[j])
    }
    adjCol = rBind(adjTemp, adjCol)
  }        
  adjM = cBind(adjCol, adjM)        
 }
 return( forceSymmetric(triu(adjM, k=1)) )
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
    startBlock = cumsum(c(0, nMembers)) + 1

    #fill the aa and xx matrices block by block
    for(i in 1:nBlocks) {
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
    return( list(adjacency = simSparseAdjMat(blockMat, nMembers),
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
simBernSparseVec = function(nRow, nCol, p) {

	#to prevent overflow problems
	nRC = round(as.numeric(nRow)*as.numeric(nCol))

    if(p == 0) {
        return( sparseVector(0, 1, length = nRC) )
    }

    expNumOnes = nRC*p
    sdNumOnes = sqrt(nRC*p*(1-p))

    indOnes = rnbinom(expNumOnes + round(3*sdNumOnes), 1, p) + 1

    indCumSum = cumsum(indOnes)
    
    while (max(indCumSum) < nRC) {
        indCumSum = c(indCumSum, max(indCumSum) +
            rnbinom(1, 1, p) + 1)
    }
    
    indCumSum = indCumSum[indCumSum <= nRC]

    return( sparseVector(1, indCumSum, nRC) )
}
