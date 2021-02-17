#############################################################################################
# Required Libraries 
#############################################################################################
library(rotl)
library(ape)
library(geiger)
library(phybase)
library(phytools)

#############################################################################################
# Get a tree from Open Tree of Life. 
#############################################################################################
phy <- get_study_tree("ot_485", "tree1")
plot(phy, cex=0.3)

#############################################################################################
# Drop Taxa from original to simplify
#############################################################################################
phy <- drop.random(phy, Ntip(phy) - 10)
plot(phy)
axisPhylo()


#############################################################################################
# Simulate Gene Tree 
#############################################################################################
gene.tree <- phybase::sim.coaltree.phylo(phy, pop.size=1e-12)
plot(gene.tree)

#############################################################################################
# Compare Gene tree to initial tree 
#############################################################################################
plot(cophylo(phy, gene.tree, cbind(sort(phy$tip.label), sort(gene.tree$tip.label))))


#############################################################################################
# Tree issues
#############################################################################################
species.tree <- rcoal(7)
species.tree$edge.length <- species.tree$edge.length / (10*max(branching.times(species.tree)))
gene.tree <- phybase::sim.coaltree.phylo(species.tree)
plot(cophylo(species.tree, gene.tree, cbind(sort(species.tree$tip.label), sort(gene.tree$tip.label))))


#############################################################################################
# Increasing the length of the species tree
#############################################################################################
tip.rows <- which(species.tree$edge[,2]<=Ntip(species.tree))
species.tree2 <- species.tree
species.tree2$edge.length[tip.rows] <- 100 + species.tree2$edge.length[tip.rows]
gene.tree2 <- phybase::sim.coaltree.phylo(species.tree2)
plot(cophylo(species.tree2, gene.tree2, cbind(sort(species.tree2$tip.label), sort(gene.tree2$tip.label))))

#############################################################################################
# Plot the Cladogram 
#############################################################################################

species.tree2.clado <- compute.brlen(species.tree2)
gene.tree2.clado <- compute.brlen(gene.tree2)
plot(cophylo(species.tree2.clado, gene.tree2.clado, cbind(sort(species.tree2.clado$tip.label),
                                                          sort(gene.tree2.clado$tip.label))))
