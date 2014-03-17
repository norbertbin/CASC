# make file to run simulation code

all: simCompMethodsGraph

simCompMethodsGraph: simCompMethodsGraph.R
	R CMD BATCH simCompMethodsGraph.R simCompMethodsGraph.Rout

simCompMethodsCov: simCompMethodsCov.R
	R CMD BATCH simCompMethodsCov.R simCompMethodsCov.R

simDiffGroupMember: simDiffGroupMember.R
	R CMD BATCH simDiffGroupMember.R simDiffGroupMember.Rout

clean:
	rm *.Rout Figs/*