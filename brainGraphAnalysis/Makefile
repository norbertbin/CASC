#make file to run the brain graph preprocessing and analysis

#define data directories for large files
nCores = 2
nH = 2
nBlocks = 100
filePre = NKI-TRT_0021001
rawDataDir = ../../../../Data/brainGraphs/rawData/
procDataDir = ../../../../Data/brainGraphs/procData/
atlasFile = ../../../../Data/brainGraphs/rawData/desikan/desikan_atlas_remap.nii

#run python script to preprocess data
convertNpyToTxt: convertNpyToTxt.py
	python convertNpyToTxt.py $(rawDataDir) $(procDataDir) $(filePre)

#run R script to compute clusters
args1 = $(rawDataDir) $(procDataDir) $(filePre) $(nBlocks) \
	$(nCores) $(nH)
clusterBrainGraph: clusterBrainGraph.R ../code/readWriteMatrix.R \
		../code/irlbaMod.R ../code/spectralClusteringMethods.R
	R CMD BATCH '--args $(args1)' clusterBrainGraph.R

#run R script to generate plots
args2 = $(procDataDir) $(filePre)
clusterPlots: clusterPlots.R ../code/helperFunctions.R \
		../code/helperFunctions.R
	R CMD BATCH '--args $(args2)' clusterPlots.R

#run script to compare clusters with brain atlas and estimate block matrix
args3 = $(atlasFile) $(rawDataDir) $(procDataDir) $(filePre)
estAtlasBlocks: estAtlasBlocks.R 
	R CMD BATCH '--args $(args3)' estAtlasBlocks.R

#run script to simulate graph and coordinates from atlas estimates
args4 = $(filePre) $(nCores) $(nH)
simAtlasBlocks: simAtlasBlocks.R
	R CMD BATCH '--args $(args4)' simAtlasBlocks.R