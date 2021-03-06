#make file to run the brain graph preprocessing and analysis

#define data directories for large files
nCores = 10
nH = 10
nBlocks = 100
rawDataDir = ../../brainData/KKI-42_MIGRAINE/data/
procDataDir = ../procData/
atlasFile = ../../brainData/desikan/desikan_atlas_remap.nii

#for a single graph
filePre = KKI-21_KKI2009-12

#for running multiple graphs
fileList = ../../brainData/KKI-42_MIGRAINE/preFixes_KKI-42.txt


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

args5 = $(rawDataDir) $(procDataDir) $(filePre) $(nBlocks) \
	$(nCores) $(nH)
clusterBrainGraphWithAtlas: clusterBrainGraphWithAtlas.R \
		../code/readWriteMatrix.R ../code/irlbaMod.R \
		../code/spectralClusteringMethods.R
	R CMD BATCH '--args $(args5)' clusterBrainGraphWithAtlas.R
##########################################
# run R scripts on a set of input graphs
aFilePre = $(shell cat $(fileList))

#run python script to preprocess data
convertNpyToTxtAll:
	$(foreach fPre, $(aFilePre), python convertNpyToTxt.py $(rawDataDir) \
	$(procDataDir) $(fPre);)

#run R script to compute clusters
clusterBrainGraphAll:
	$(foreach fPre, $(aFilePre), R CMD BATCH '--args $(rawDataDir) $(procDataDir) $(fPre) $(nBlocks) $(nCores) $(nH)' clusterBrainGraph.R;)

#run R to more efficiently cluster multiple graphs in parallel
clusterBrainGraphList:
	R CMD BATCH '--args $(rawDataDir) $(procDataDir) $(fileList) $(nBlocks) $(nCores) $(nH)' clusterBrainGraphList.R

#run R script to generate plots
clusterPlotsAll:
	$(foreach fPre, $(aFilePre), R CMD BATCH '--args $(procDataDir) $(fPre)' clusterPlots.R;)	

#run R script to compare paired clusters
compareBrainClusters:
	R CMD BATCH '--args $(fileList) $(procDataDir)' compareBrainClusters.R

#run R script to plot node degree histograms
plotDegreeDist:
	R CMD BATCH '--args $(rawDataDir) $(procDataDir) $(fileList) $(nCores)' plotDegreeDist.R
