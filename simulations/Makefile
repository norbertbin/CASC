# make file to run simulation code

all: simCompMethodsGraphDelta simCompMethodsGraphN simCompMethodsCovDelta simCompMethodsCovN simDiffGroupMember

simCompMethodsGraphDelta: simCompMethodsGraphDelta.R
	R CMD BATCH simCompMethodsGraphDelta.R simCompMethodsGraphDelta.Rout

simCompMethodsGraphN: simCompMethodsGraphN.R
	R CMD BATCH simCompMethodsGraphN.R simCompMethodsGraphN.Rout

simCompMethodsCovDelta: simCompMethodsCovDelta.R
	R CMD BATCH simCompMethodsCovDelta.R simCompMethodsCovDelta.Rout

simCompMethodsCovN: simCompMethodsCovN.R
	R CMD BATCH simCompMethodsCovN.R simCompMethodsCovN.Rout

simDiffGroupMember: simDiffGroupMember.R
	R CMD BATCH simDiffGroupMember.R simDiffGroupMember.Rout

makePlots: makePlots.R
	R CMD BATCH makePlots.R makePlots.Rout

clean:
	rm -f *.Rout figs/*
