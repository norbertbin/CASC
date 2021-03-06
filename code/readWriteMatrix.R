#############################################################
#functions to load and save a matrix list into a binary file
saveMatrixList <- function(baseName, mtxList) {
    idxName <- paste(baseName, ".idx", sep="")
    idxCon <- file(idxName, 'wb')
    on.exit(close(idxCon))

    dataName <- paste(baseName, ".bin", sep="")
    con <- file(dataName, 'wb')
    on.exit(close(con))

    writeBin(0L, idxCon)

    for (m in mtxList) {
        writeBin(dim(m), con)
        writeBin(typeof(m), con)
        writeBin(c(m), con) 
        flush(con)

        offset <- as.integer(seek(con))
        cat('offset', offset)
        writeBin(offset, idxCon)
    }

    flush(idxCon)
}

loadMatrix <- function(baseName = "data", index) {
    idxName <- paste(baseName, ".idx", sep="")
    idxCon <- file(idxName, 'rb')
    on.exit(close(idxCon))

    dataName <- paste(baseName, ".bin", sep="")
    con <- file(dataName, 'rb')
    on.exit(close(con))

    seek(idxCon, (index-1)*4)
    offset <- readBin(idxCon, 'integer')

    seek(con, offset)
    d <- readBin(con, 'integer', 2)
    type <- readBin(con, 'character', 1)
    structure(readBin(con, type, prod(d)), dim=d)
}

# ---------------------------------------------------------------------
# Functions to read/write sparse matrix into a binary file
# ---------------------------------------------------------------------
# TODO: handle symmetric case explicitly 
saveSparseMatrix = function(baseName, sMat) {
    wSMat = as.matrix(summary(sMat))
    wSMat = rbind(wSMat, c(sMat@Dim, 0))

    saveMatrixList(baseName, list(wSMat))
}

loadSparseMatrix = function(baseName) {
    rSMat = loadMatrix(baseName, 1)

    nrows = dim(rSMat)[1]
    dims = rSMat[nrows, 1:2]
    rSMat = rSMat[-nrows,]

    return(sparseMatrix(i = rSMat[,1],
                        j = rSMat[,2],
                        x = rSMat[,3],
                        dims = dims))
}
