library(ggplot)
library(ape)
library(cluster)

set.seed(13)

tree <- ape::read.tree("/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/Ce_phy_file_LD_0.9.phy.contree")

# Patristic distance between all strains
D <- ape::cophenetic.phylo(tree)

# PAM (k-medoids) chooses k representative leaves (actual strains)
k <- 64
pam_fit <- pam(as.dist(D), k = k, diss = TRUE)

core64 <- rownames(D)[pam_fit$medoids]