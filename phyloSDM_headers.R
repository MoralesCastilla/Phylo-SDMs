############################################
# Functions to:
# * Simulate community data matrices 
# * Compute Phylogenetic predictors
# * Compute Trait predictors (similar logic as Phylogenetic predictors)
# 
# Ignacio Morales-Castilla, T. Jonathan Davies, William D. Pearse & Pedro Peres-Neto
#
# Please cite:
# "Combining phylogeny and co-occurrence to improve single species distribution models"
#
#
############################################

#' Simulate site-by-species community data matrices
#'
#' The probability of presence of species jth in site ith is a
#' function of an environmental gradient, a trait, and the
#' interactions between the two. Weighting parameters (mu, gamma, and
#' omega) define the importance of each of these three aspects
#' respectively.
#'
#' @param env numerical vector (same length as number of sites) of
#' environmental gradient
#' @param trait numerical vector of species' traits (same length as
#' number of species)
#' @param mu weight for environmental gradient (env)
#' @param gamma weight for species' traits (trait)
#' @param omega weight for interaction of environmental gradient and species' traits
#' @return site-by-species matrix
#' @examples
#' env.trait.additive <- generate_community(rnorm(100), rnorm(20), c(1, 1, 0)
generate_community <- function(env, trait, mu, gamma, omega){
    #Setup
    n.sites <- length(env)
    n.spp <- length(trait)

    #Do it
    L <- matrix(data=0,nrow=n.sites,ncol=n.spp)
    for(j in seq_along(trait))
        L[,j] <- rbinom(n.sites,1,(1/(1+exp(-mu*env+gamma*trait[j]+omega*env*trait[j]))))
    return(L)
}

#' Compute predictors (of phylogeny, traits, etc.)
#'
#' The predictors are determined by a community data matrix (L), the
#' phylogenetic tree of all species or functional trait distances
#' (sp.info, and the choice of phylogenetic distance measure (see
#' 'type' and 'alpha'). The user must also specify the target species
#' ('target.sp').
#' @param L site-by-species matrix
#' @param sp.info Either an ape::phylo phylogenetic tree, or a matrix
#' of distances among species (based on functional traits, phylogeny,
#' etc.)
#' @param name of the target species for which to calculate the
#' phylogenetic predictors
#' @param target.sp target species for which to calculate the distance
#' predictor
#' @param type either 'original' (use equation 1 in Morales-Castilla
#' et al.) or 'transform' (use equation 2 in Morales-Castilla et
#' al.). See 'alpha'.
#' @param alpha vector of length 2; determines how 'alpha' is
#' calculated. First element is the first value, second value is the
#' spacing between values. Maximum value is determined by phylogenetic
#' distance matrix and cannot be altered.
#' @return list with first element of predictor(s), and a second of
#' any site with only sites dropped from L
generate_predictor <- function(L, sp.info, target.sp, type=c("original","transform"), alpha=c(0.1,0.2)){
    #Setup
    if(!inherits(sp.info, c("dist","phylo")))
       stop("'sp.info' must be either a phylogeny or a species' distance matrix")
    if(is.phylo(sp.info)) dist <- cophenetic.phylo(sp.info) else dist <- as.matrix(sp.info)
    if(nrow(as.matrix(dist)) != ncol(L))
        stop("'sp.info' must be of same dimensions as L")
    n.sites <- nrow(L)
    n.species <- ncol(L)
    colnames(L) <- colnames(dist)
    type <- match.arg(type)
    distance.with.other.species <- dist[target.sp,] # vector of distances of the target species with others  
    distribution.other.species <- L[,-which(colnames(L)==target.sp)]
    distance.with.other.species <- distance.with.other.species[-which(names(distance.with.other.species)==target.sp)]
    sites.with.zero.sums <- which(apply(distribution.other.species,1,sum)==0)
    if (length(sites.with.zero.sums)>0)
        distribution.other.species <- distribution.other.species[-sites.with.zero.sums,]

    #Do the work
    if (type=="original"){
        Phylo.predictor <- distribution.other.species%*%(1/distance.with.other.species)/apply(distribution.other.species,1,sum)} else {
            alpha <- seq(alpha[1],max(distance.with.other.species),alpha[2])
            Phylo.predictor <- matrix(0,n.sites,length(alpha))
            for (k in seq_along(alpha))
                Phylo.predictor[,k] <- distribution.other.species%*%(exp(-distance.with.other.species/alpha[k]))/apply(distribution.other.species,1,sum)
        }
    
    result <- list(predictor=Phylo.predictor,sites.with.zero.sums=sites.with.zero.sums)
    return(result)
}



