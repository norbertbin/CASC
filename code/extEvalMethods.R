# ---------------------------------------------------------------------
# functions to calculate misclustering rate
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# MAIN FUNCTIONS
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# function to calculate the empirical proportion of misclustered nodes
# ---------------------------------------------------------------------
misClustRateEmp = function(singVec, nMembers) {

  nNodes = sum(nMembers)
  nBlocks = length(nMembers)
  clusters = kmeans(singVec, nBlocks, nstart = 10)$cluster
  nMisClustNodes = 0
  uniqueClusters = unique(clusters)
  clusterCounts = matrix(rep(0, nBlocks*nBlocks),
      nrow = nBlocks)
  clusterLabels = rep(0, nBlocks)
  
  #get label counts for each cluster
  for(i in 1:nBlocks) {
    
    clustStart = sum(nMembers[1:i]) - sum(nMembers[i]) + 1
    clustEnd = sum(nMembers[1:i])
    
    for(j in uniqueClusters) {
      clusterCounts[j,i] = sum(j == clusters[clustStart:clustEnd]) 
    }    
  }
  
  #determine cluster label based on counts
  clusterCountsTemp = clusterCounts
  for(i in 1:nBlocks) {
    
    maxCoor = t(which(clusterCountsTemp ==
        max(clusterCountsTemp), arr.ind = T))
    clusterLabels[maxCoor[2]] = maxCoor[1]
    clusterCountsTemp[maxCoor[1], ] = rep(-1, nBlocks)
    clusterCountsTemp[, maxCoor[2]] = rep(-1, nBlocks)
    
  }
  
  for(i in 1:nBlocks) {
    nMisClustNodes = nMisClustNodes + sum(clusterCounts[-clusterLabels[i],i])
  }  
  
  
  return( nMisClustNodes/nNodes )
}

# ---------------------------------------------------------------------
# function to calcuale the proportion of misclustered nodes based on
# singular vectors with condition that if 95% are in one cluster all
# are misclustered
# ---------------------------------------------------------------------
misClustRate = function(singVec, nMembers) {

  clusters = getClustUsingRotMat(singVec, nMembers)
 
  #check if 95% of clusters are in one cluster
  if(sum(table(clusters)/sum(nMembers) >= .95) == 1) {
      return(1)
  }

  nBlocks = length(nMembers)
  nNodes = sum(nMembers)
  nMisClustNodes = 0
  blockMembership = rep.int(1:nBlocks, times = nMembers)
  
  for(i in 1:n) {
    if(clusters[i] != blockMembership[i]) { 
      nMisClustNodes = nMisClustNodes + 1
    }
  }
  
  return( nMisClustNodes/nNodes )
}


# ---------------------------------------------------------------------
# HELPER FUNCTIONS
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# finds the rotation matrix for the eigenvectors
# ---------------------------------------------------------------------
getRotMat = function(singVec, nMembers) {

  nBlocks = length(nMembers)
  nNodes = sum(nMembers)

  #generate assignment vectors
  z = diag(rep(1, nBlocks))
  clusters = rep(0, nNodes)
  blockMembership = rep.int(1:nBlocks, times = nMembers)
  
  #get rotation matrix
  x = matrix(rep(0,sum(nMembers)*nBlocks), nrow = nMembers)
  for(i in 1:n) {
    x[i,] = z[blockMembership[i],] 
  }
  svdRes = svd(t(x)%*%singVec)
  rotMat = svdRes$v%*%t(svdRes$u)

  return(rotMat)
}


# ---------------------------------------------------------------------
# cluster using eigenvectors and the rotation matrix
# ---------------------------------------------------------------------
getClustUsingRotMat = function(singVec, nMembers) {

  nBlocks = length(nMembers)
  nNodes = sum(nMembers)
  #generate assignment vectors
  z = diag(rep(1, nBlocks))
  clusters = rep(0, nNodes)
  blockMembership = rep.int(1:nBlocks, times = nMembers)
  
  #get rotation matrix
  x = matrix(rep(0,sum(nMembers)*nBlocks), nrow = nNodes)
  for(i in 1:nNodes) {
    x[i,] = z[blockMembership[i], ] 
  }
  svdRes = svd(t(x)%*%singVec)
  rotMat = svdRes$v%*%t(svdRes$u)
  
  
  for(i in 1:nNodes) {
    blockDist = rep(0, nBlocks)
    for( j in 1:nBlocks) {
      blockDist[j] = sum((singVec[i,]%*%rotMat-z[j,])^2)
    }
    clusters[i] = match(min(blockDist), blockDist)[1]
  }

  return(clusters)
}
