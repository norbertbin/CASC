# ---------------------------------------------------------------------
# use txt files with simulation results to regenerate plots
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# load required scripts and packages
# ---------------------------------------------------------------------
if (!require(lattice)) {
    install.packages('lattice', dependencies = T)
    require(lattice)
}

outDir = "Figs/"

# ---------------------------------------------------------------------
# read in data from txt files
# ---------------------------------------------------------------------
datPQ = read.table(paste(outDir, "simCompMethodsPQ.txt", sep=""), header=T)
datDense = read.table(paste(outDir, "simCompMethodsDense.txt", sep=""),
    header=T)
datCovProb = read.table(paste(outDir, "simCompMethodsCovProb.txt", sep=""),
    header=T)
datNumCov = read.table(paste(outDir, "simCompMethodsNumCov.txt", sep=""),
    header=T)
datMember = read.table(paste(outDir, "simDiffGroupMember.txt", sep=""),
    header=T)

# ---------------------------------------------------------------------
# generate plots
# ---------------------------------------------------------------------
#legend labels
#legLab = c("Covariate Assisted Spectral Clustering",
#    "Canonical Correlation Analysis",
#    "Regularized Spectral Clustering",
#    "Spectral Clustering on X")

legLab = c("CASC", "CCA", "RSC", "SC-X")

#plot for datPQ
pdf(paste(outDir, "simCompMethodsPQ.pdf", sep=""), width = 7, height = 7)
print(
xyplot(misClustRate ~ deltaPQ,
       group = group,
       type = "b",
       pch = 1:4,
       cex = 1.2,
       data = datPQ,
       ylab = list( label = "Average mis-clustering rate",
                    cex = 1.7),
       xlab = list( label = "Within minus between block probability (p - q)",
                    cex = 1.7),
       main = list("(a)", cex=1.7),
       lwd = 2,
       key = list( text = list(legLab,
                               cex = 1.5),
                   lines = list( col = trellis.par.get()$superpose.symbol$col[1:4],
                                 lwd = 2),
                    points = list( pch = 1:4,
                                   cex = 1.2,
                                   col = trellis.par.get()$superpose.symbol$col[1:4]),
                                   corner = c(.95, .95)),
       scales = list( cex = 1.5))
)
dev.off()

#plot for datDense
pdf(paste(outDir, "simCompMethodsDense.pdf", sep=""), width = 7, height = 7)
print(
    xyplot(misClustRate ~ nNodes,
           group = group,
           type = "b",
           pch = 1:4,
           cex = 1.2,
           data = datDense,
           ylab = list( label = "Average mis-clustering rate",
                        cex = 1.7),
           xlab = list( label = "Number of nodes (N)",
                        cex = 1.7),
           main = list( label = "(b)",
                        cex = 1.7),
           lwd = 2,
           key = list( text = list(legLab,
                                   cex = 1.5),
                       lines = list( col = trellis.par.get()$superpose.symbol$col[1:4],
                                     lwd = 2),
                       points = list( pch = 1:4,
                                      cex = 1.2,
                                      col = trellis.par.get()$superpose.symbol$col[1:4]),
                       corner = c(.95, .95)),
           scales = list( cex=1.5))
    )
dev.off()

#plot for datCovProb
pdf(paste(outDir, "simCompMethodsCovProb.pdf", sep=""), width = 7, height = 7)
print(
    xyplot(misClustRate ~ deltaCovProb,
           group = group,
           type = "b",
           pch = 1:4,
           cex = 1.2,
           data = datCovProb,
           ylab = list( label = "Average mis-clustering rate",
                        cex = 1.7),
           xlab = list( label = expression('Difference in covariate probabilities (m'[1]*' - m'[2]*')'),
                         cex = 1.7),
           main = list( label = "(a)",
                        cex = 1.7),
           lwd = 2,
           key = list( text = list(legLab,
                                   cex = 1.5),
                       lines = list( col = trellis.par.get()$superpose.symbol$col[1:4],
                                     lwd = 2),
                       points = list( pch = 1:4,
                                      cex = 1.2,
			                          col = trellis.par.get()$superpose.symbol$col[1:4]),
                       corner = c(.95, .95)),
           scale = list(cex=1.5))
    )
dev.off()

#plot for datNumCov
pdf(paste(outDir, "simCompMethodsNumCov.pdf", sep=""), width = 7, height = 7)
print(
    xyplot(misClustRate ~ nCov,
           group = group,
           type = "b",
           pch = 1:4,
           cex = 1.2,
           data = datNumCov,
           ylab = list( label = "Average mis-clustering rate",
                        cex = 1.7),
           xlab = list( label = "Number of Covariates (R)",
                        cex = 1.7),
           main = list( label = "(b)",
                        cex = 1.7),
           lwd = 2,
           key = list( text = list(legLab,
                                   cex = 1.5),
                       lines = list(col = trellis.par.get()$superpose.symbol$col[1:4],
                                    lwd = 2),
                       points = list(pch = 1:4,
                                     cex = 1.2,
			                         col = trellis.par.get()$superpose.symbol$col[1:4]),
                       corner = c(.95, .95)),
           scales = list(cex=1.5))
    )
dev.off()

#plot for datMember
pdf(paste(outDir, "simDiffGroupMember.pdf", sep=""), width = 7, height = 7)
print(
    xyplot(misClustRate ~ propIncorrect,
           group = group,
           type = "b",
           pch = 1:4,
           cex = 1.2,
           data = datMember,
           ylab = list( label = "Average mis-clustering rate",
                        cex = 1.7),
           xlab = list( label="Covariate to graph block membership agreement",
                        cex = 1.7),
           lwd = 2,
           key = list( text = list(legLab,
                                   cex = 1.5),
                       lines = list(col = trellis.par.get()$superpose.symbol$col[1:4],
                                    lwd = 2),
                       points = list(pch = 1:4,
                                     cex = 1.2,
			                         col = trellis.par.get()$superpose.symbol$col[1:4]),
                       corner = c(.95, .95)),
           scales = list(cex=1.5))
    )
dev.off()
