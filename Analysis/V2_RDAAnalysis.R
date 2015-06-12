###############
# V2 Phytocosms
# RDA Analysis
# Quentin Schorpp
# 02.06.2015
##############

source("Analysis/V2_RequiredPackages.R")
source("Data/V2_LoadData.R")

# Sample 11 has strange Cmic values
full.frame <- full.frame[-11,]

# Divide data frame in subsets ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
soil <- full.frame[,c(12,14:18)]
bio <- full.frame[,8:12]
env <- full.frame[,c(8:12,14:18)]
groups <- full.frame[,1:7]
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
summary(soil)
is.na(soil) <- 


# Z-transformation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
soil$NO3 <- soil$NO3*(-1)
for (i in 1:10) {
        na.action=na.omit(soil)
        soil[,i] <- soil[,i] + abs(min(soil[,i], na.rm=TRUE))   
}

soil.log <- log1p(soil)
soil.z <- data.frame(scale(soil.log))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Test for normal distribution ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par(mfrow=c(3,4))
for (i in 1:10) {
        qqnorm(soil.z[,i], xlab=colnames(soil.z[i]))
        qqline(soil.z[,i])
        
           
}

par(mfrow=c(3,4))
for (i in 1:10) {
        hist(soil.z[,i], xlab=colnames(soil.z[i]))
}
summary(soil.z)

soil.mnorm <- t(soil.z)
mshapiro.test(soil.mnorm) # NO Normal Distribution!!!!!
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

soil.z[is.na(soil.z)] <- 0 

# PCA
soil.pca <- rda(soil.z, scale=FALSE)
summary(soil.pca)
summary(soil.pca, scaling=1)


# Inertia: in vegan's language this is the term for "variation" in the data. This term comes from the word of CA (Sect. 5.4). 
# In PCA, the "inertia" is either the sum of the variances of the variables (PCA on a covariance matrix) or the sum of the diagonal values 
# of the correlation matrix (PCA on a correlation matrix), i.e. the sum of all correlations of the variables with themselves, 
# which corresponds to the number of variables.

# Constrained and unconstrained (constrained = gezwungen, abhängig): In PCA the  analysis is unconstrained and so are the results.

# Eigenvalues: symbolized as lambda, these are measures of the importance of an axis. By dividing them by the total inertia, they can be expressed as 
# Proportion explained / proportion of variation accounted for 

# Scaling: 1 = distance biplot: the eigenvectors are scaled to unit length. (1) Distances among objects in the biplot are approximations
# of their Euclidean distances in multidimensional space. (2) The angles among descriptor vectors are meaningless
#          2 = correlation biplot: ech eigenvector is scaled to the squareroot of the eigenvalue. (1) Distances among objects are not 
# approximations of their Euclidean distances in multidimensional space. (2) The angles between descriptors  in the biplot reflect their
# correlations

# Species scores: coordinates of the arrow heads of the variables. For historical rasons response variables are always calles "species" in vegan

# Site scores: coordinates of the sites in the ordination diagram. Objects are always called "Sites" in vegan output files.




# Loads of variables:
scores(soil.pca, choices = 1:4, display = "species", scaling = 0)

# Biplot

source("Analysis/evplot.R")
source("Analysis/cleanplot.pca.R")

eigenvalues <- soil.pca$CA$eig
evplot(eigenvalues) # One interprets only the axes whose eigenvalues are larger than the length of corresponding piece of the sticks


par(mfrow=c(1,2))
biplot(soil.pca, scaling=1,  main="PCA - Scaling 1")
pcacircle(soil.pca)
biplot(soil.pca, scaling=2,  main="PCA - Scaling 2")

cleanplot.pca(soil.pca, ahead=0)

par(mfrow=c(1,1))
pl <- biplot(soil.pca, type = "points")
identify(pl,"sp", plot=TRUE, atpen=TRUE, labels=names(soil.z), cex=1.0)

points(pl, "sites", pch=25, bg=interaction(groups$soil, groups$plant), cex=0.7)

rda(env, scale=TRUE)
