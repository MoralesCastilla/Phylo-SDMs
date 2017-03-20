############################################
#
# Script with an example to test statistical performance of phylo-SDMs vs. non-phylo-SDMs
# using functions to:
#   a) Simulate community data matrices 
#   b) Compute Phylogenetic predictors
#   c) Compute Trait predictors (similar logic as Phylogenetic predictors)
#
#
#
# Ignacio Morales-Castilla, Jonathan Davies, William D. Pearse & Pedro Peres-Neto
############################################




####################################
# Load packages and functions
####################################
lapply(c("caper","phytools","geiger","mvMORPH","dismo","randomForest","SDMTools"), require, character.only=T)
source('phyloSDM_headers.R')

####################################
# Simulate and plot data
####################################
# Initial values
nsps <- 100
nsites <- 250
tree <- sim.bdtree(n=nsps) 
mu <- 2     #weight to Environment, E
gamma <- 2  #weight to Trait, T (phylogenetic or random)
omega <- 0  #weight to E*T


# Simulate Environmental gradients and phylogenetic and non-phylogenetic traits
E <- matrix(rnorm(nsites),nrow=nsites,ncol=1) # environmental matrix
Tr <-  scale(rTraitCont(rescale(tree, "delta", 0.01))) # very conserved
Tr2 <-  scale(rTraitCont(rescale(tree, "delta", 10))) # very labile
Tr3 <-  scale(rTraitCont(rescale(tree, "delta", 1))) # Brownian
Tr.rnd <-  as.numeric(scale(rnorm(nsps))) # random non-phylogenetic


# Simulate community data matrices according to each trait 
L <- generate_community(E, Tr, mu, gamma, omega)
L2 <- generate_community(E, Tr2, mu, gamma, omega)
L3 <- generate_community(E, Tr3, mu, gamma, omega)
L.rnd <- generate_community(E, Tr.rnd, mu, gamma, omega)

#Plot it out
phy.dendro  <-  as.dendrogram(as.hclust.phylo(tree))

heatmap(L,scale="none", Rowv=E,Colv=phy.dendro,col=c("white","black"),labCol=FALSE,labRow=FALSE)
heatmap(L2,scale="none", Rowv=E,Colv=phy.dendro,col=c("white","black"),labCol=FALSE,labRow=FALSE)
heatmap(L3,scale="none", Rowv=E,Colv=phy.dendro,col=c("white","black"),labCol=FALSE,labRow=FALSE)
heatmap(L.rnd,scale="none", Rowv=E,Colv=phy.dendro,col=c("white","black"),labCol=FALSE,labRow=FALSE)

################################################
# Calculate phylogenetic and trait predictors
################################################
#Assume we're working with species 'j'
j <- 1
name.target.species <- tree$tip.label[j]

# computing phylogenetic predictor for species j within L
Phylo.pred <- generate_predictor(L,tree,name.target.species,type="transform")
PhyVector <- scale(Phylo.pred$predictor) # Phylogenetic predictors at varying alpha levels P

# computing non-phylogenetic random predictor for species j within L
names(Tr.rnd) <- tree$tip.label
name.target.species <- names(Tr.rnd)[j]
trait.dists <- dist(Tr.rnd)
Trait.pred <- generate_predictor(L,trait.dists,name.target.species,type="transform")
TraitVector <- scale(Trait.pred$predictor) # Random predictors at varying alpha levels R

# model and compare explained deviance
colnames(L) <- tree$tip.label
L <- L[which(rowSums(L)>0),]
  
  
D2.phylo.rnd <- array(0,dim=c(ncol(TraitVector),2))
                      
for(i in 1:ncol(TraitVector)){
    phylo.SDM1 <- glm(L[,name.target.species]~PhyVector[,i], family="binomial",na.action=na.omit) # phylo SDM
    rnd.SDM1 <- glm(L[,name.target.species]~TraitVector[,i], family="binomial",na.action=na.omit) # non-phylo SDM
    
    D2.phylo.rnd[i,1] <- 1-phylo.SDM1$deviance/phylo.SDM1$null.deviance
    D2.phylo.rnd[i,2] <- 1-rnd.SDM1$deviance/rnd.SDM1$null.deviance
}                         
                           
# plot distribution of explained deviances by P and R across different alpha levels
boxplot(D2.phylo.rnd)

