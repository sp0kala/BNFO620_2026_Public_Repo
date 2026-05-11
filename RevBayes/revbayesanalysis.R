#Rev Bayes Script for BNFO620:

####All of this code is based of this pipeline
### :https://revbayes.github.io/tutorials/sequential_bayes/unrooted_gene_trees.html

#reading in sequences
genes <- readDiscreteCharacterData("concatenated.nex")

#Variables of interest

##number of taxa
num_taxa <- genes.ntaxa()
##number of branches
num_branches <- 2 * num_taxa - 3
##Taxa
taxa <- genes.taxa()

##Vectors for moves and monitors
moves    = VectorMoves()
monitors = VectorMonitors()

#Substitution Model
#A: 36.117%, C: 18.322%, T: 34.723%, G: 10.836%

##stationary frequency parameters
pi_prior <- v(0.36117, 0.18322, 0.34723, 0.10836)
pi ~ dnDirichlet(pi_prior)
moves.append( mvBetaSimplex(pi, weight=2.0) )
moves.append( mvDirichletSimplex(pi, weight=1.0) )

##Exchangeability rate parameters

#AC, AG, AT, CG, CT, GT
#transitions (AG, CT) favored
er_prior <- v(1, 3, 1, 1, 3, 1) 
er ~ dnDirichlet(er_prior)
moves.append( mvBetaSimplex(er, weight=3.0) )
moves.append( mvDirichletSimplex(er, weight=1.5) )

##deterministic rate matrix variable
Q := fnGTR(er,pi)

##gamma distribution
alpha ~ dnUniform( 0, 1.0 )
alpha.setValue(1.0)

sr := fnDiscretizeGamma( alpha, alpha, 4, false )
moves.append( mvScale(alpha, weight=2.0) )


#Tree Model

#outgroup
outgroup = clade("G_borealis")

#topology prior
topology ~ dnUniformTopology(taxa)

#Branch length prior
for (i in 1:num_branches){
  branch_len[i] ~ dnExponential(10.0)
  moves.append ( mvScale(branch_len[i], weight=1.0) )
}

#Tree Assembly
psi := treeAssembly(topology, branch_len)

#Tree Length
TL :=psi.treeLength()

##Continuous-Time Markov Chain (CTMC) Model 
# the sequence evolution model
seq ~ dnPhyloCTMC(tree=psi, Q=Q, siteRates=sr, type="DNA", rootFrequencies=pi)

#attaching data to variables
seq.clamp(genes)

#looking at analysis

#model with all variables
mymodel = model(psi)

#tree moves for topology
moves.append( mvNNI(topology, weight=num_taxa) )
moves.append( mvSPR(topology, weight=num_taxa/5.0) )

#specifying monitored variables
monitors.append( mnScreen(alpha, TL, printgen=1000) )
monitors.append( mnFile(psi, filename="output/concatenated.trees", printgen=100) )
monitors.append( mnModel(filename="output/concatenated.log", printgen=100) )

#running the analysis
mymcmc = mcmc(mymodel, moves, monitors, nruns=2, combine="mixed")
mymcmc.run(generations=300000, tuningInterval=200)

#make the treee!!
treetrace = readTreeTrace("output/concatenated.trees", treetype="non-clock", burnin=0.25)
consensus_tree = consensusTree(treetrace,"output/concatenated_MAP.tre")


#####This was ran in terminal with the following script:
# cd to my wd
# rb revbayesanalysis.R


