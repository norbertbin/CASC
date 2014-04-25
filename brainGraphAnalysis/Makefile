#make file to run the brain graph preprocessing and analysis

#define data directories for large files
nCores = 12
nH = 12
nBlocks = 100
rawDataDir = ../../brainData/NKI-24_MIGRAINE/data/
procDataDir = ../procData/
atlasFile = ../../brainData/desikan/desikan_atlas_remap.nii

#for a single graph
filePre = NKI-TRT_0021001

#for running multiple graphs
fileList = ../../brainData/NKI-24_MIGRAINE/preFixes_NKI-24.txt


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

#run R script to generate plots
clusterPlotsAll:
	$(foreach fPre, $(aFilePre), R CMD BATCH '--args $(procDataDir) $(fPre)' clusterPlots.R)	
