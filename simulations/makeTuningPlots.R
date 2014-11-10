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

outDir = "Figs"

# ---------------------------------------------------------------------
# read in data from txt files
# ---------------------------------------------------------------------
datPQ = read.table(paste(outDir, "simCompTuningPQ.txt", sep=""), header=T)
datDense = read.table(paste(outDir, "simCompTuningDense.txt", sep=""),
    header=T)

# ---------------------------------------------------------------------
# generate plots
# ---------------------------------------------------------------------
#legend labels

legLab = c("CASC Enhanced", "CASC Simple")

#plot for datPQ
pdf(paste(outDir, "simCompTuningPQ.pdf", sep=""), width = 7, height = 7)
print(
xyplot(misClustRate ~ deltaPQ,
       group = group,
       type = "b",
       pch = 1:2,
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
                   lines = list( col = trellis.par.get()$superpose.symbol$col[1:2],
                                 lwd = 2),
                    points = list( pch = 1:2,
                                   cex = 1.2,
                                   col = trellis.par.get()$superpose.symbol$col[1:2]),
                                   corner = c(.95, .95)),
       scales = list( cex = 1.5))
)
dev.off()

#plot for datDense
pdf(paste(outDir, "simCompTuningDense.pdf", sep=""), width = 7, height = 7)
print(
    xyplot(misClustRate ~ nNodes,
           group = group,
           type = "b",
           pch = 1:2,
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
                       lines = list( col = trellis.par.get()$superpose.symbol$col[1:2],
                                     lwd = 2),
                       points = list( pch = 1:2,
                                      cex = 1.2,
                                      col = trellis.par.get()$superpose.symbol$col[1:2]),
                       corner = c(.95, .95)),
           scales = list( cex=1.5))
    )
dev.off()


