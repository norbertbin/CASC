#make file to run the brain graph preprocessing and analysis

#define data directories for large files
nCores = 2
nH = 2
nBlocks = 100
filePre = NKI-TRT_0021001
rawDataDir = ../../../../Data/brainGraphs/rawData/
procDataDir = ../../../../Data/brainGraphs/procData/

#run python script to preprocess data
convertNpyToTxt: convertNpyToTxt.py
	python convertNpyToTxt.py $(rawDataDir) $(procDataDir) $(filePre)

#run R script to compute clusters
args = $(rawDataDir) $(procDataDir) $(filePre) $(nBlocks) \
	$(nCores) $(nH)
clusterBrainGraph: clusterBrainGraph.R ../code/readWriteMatrix.R \
		../code/irlbaMod.R ../code/spectralClusteringMethods.R
	R CMD BATCH '--args $(args)' clusterBrainGraph.R

#run R script to generate plots
args = $(procDataDir) $(filePre)
clusterPlots: clusterPlots.R ../code/helperFunctions.R \
		../code/helperFunctions.R
	R CMD BATCH '--args $(args)' clusterPlots.R