# make file to run simulation code

all: simCompMethodsGraph simCompMethodsCov simDiffGroupMember

simCompMethodsGraph: simCompMethodsGraph.R
	R CMD BATCH simCompMethodsGraph.R simCompMethodsGraph.Rout

simCompMethodsCov: simCompMethodsCov.R
	R CMD BATCH simCompMethodsCov.R simCompMethodsCov.Rout

simDiffGroupMember: simDiffGroupMember.R
	R CMD BATCH simDiffGroupMember.R simDiffGroupMember.Rout

makePlots: makePlots.R
	R CMD BATCH makePlots.R makePlots.Rout

clean:
	rm -f *.Rout figs/*
