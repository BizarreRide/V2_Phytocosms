###############
# V2 Phytocosms
# RDA Analysis
# Quentin Schorpp
# 02.06.2015
##############

# Description:
# This is an attempt for mutlivariate analyses of Phytocosms from the phytocosms experiment of Valentina Sandor
# All measures are calculated as Differences between beginning and end of the experiment, i.e. changes
        


source("Analysis/V2_RequiredPackages.R") #packages agricolae and vegan are colliding
source("Data/V2_LoadData.R")

# Sample 11 has strange Cmic values
full.frame <- full.frame[-11,] # Phytocosm 11 was strange?

# Divide data frame in subsets ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
soil <- full.frame[,c(12,14:18)] # abiotic soil measures, including Cmic
colnames(soil)[2:3] <- c("carbon","Nitrogen")
bio <- full.frame[,8:12] # biotic soil characteristics
env <- full.frame[,c(8:12,14:18)] # the same as phyto
groups <- full.frame[,1:7]
phyto <- cbind(groups,bio, soil[,2:6])  # merge of bio and soil variables 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Soil ####        
## Z-transformation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Here is an approach, that tries to handle negative values in some variables
# Since No3 is strictly negative, i just turn it inot a positive variable
# Then all variables get their minimum added,...

#soil$NO3 <- soil$NO3*(-1)
#for (i in 1:6) {
#        na.action=na.omit(soil)
#        soil[,i] <- soil[,i] + abs(min(soil[,i], na.rm=TRUE))   
#}

#soil.log <- log1p(soil)
#soil.z <- data.frame(scale(soil.log))
#soil.z[is.na(soil.z)] <- 0

# However, without caring for multinormality:
# soil.z <- data.frame(scale(soil))
#  -> rda() can do this 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Test for normal distribution ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par(mfrow=c(3,4))
for (i in 1:6) {
        qqnorm(soil.z[,i], xlab=colnames(soil.z[i]))
        qqline(soil.z[,i])
        
           
}

par(mfrow=c(3,4))
for (i in 1:6) {
        hist(soil.z[,i], xlab=colnames(soil.z[i]))
}
summary(soil.z)

par(mfrow=c(1,1))

scatter.smooth(soil$NH4)
ggplot(full.frame, aes(y=NH4,x=c(1:79), col=plant, shape=treat)) + geom_point(size=8)

soil.mnorm <- t(soil.z)
mshapiro.test(soil.mnorm) # NO Multi-Normal Distribution!!!!!
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## PCA of soil variables ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# PCA
soil.pca <- rda(soil.z, scale=FALSE)
summary(soil.pca) # scaling 2 by default
summary(soil.pca, scaling=1)

# Loads of variables:
scores(soil.pca, choices = 1:4, display = "species", scaling = 0)

## Some Explanations of PCA Terms: ####
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

## Biplots ####

source("Analysis/evplot.R")
source("Analysis/cleanplot.pca.R")

eigenvalues <- soil.pca$CA$eig
evplot(eigenvalues) # One interprets only the axes whose eigenvalues are larger than the length of corresponding piece of the sticks

cleanplot.pca(soil.pca, ahead=0)

par(mfrow=c(1,2))
biplot(soil.pca, scaling=1,  main="PCA - Scaling 1")
pcacircle(soil.pca)
biplot(soil.pca, scaling=2,  main="PCA - Scaling 2")

### Interpretation ####
# The Proportion of variance accounted for by the first two axes is 0.7375 or 73.8% This high value makes us confident 
# that our interpretation of the first pair of axes extracts most relevant information from the data. 

# First, the scaling 1 biplot displays a feature that must be explained. The circle is called a circle of equilibrium contribution. Its radius is equal to sqrt(d/p), 
# where d is the number of axes represented in the biplot (usually d=2) and p is the number of dimensions of the PCA space (i.e. the number of variables of the data matrix,
#                                                                                                                          6 in this case)
# The radius of this circle represents the length of the vector representing a variable that contributes equally to all the dimanesions of the PCA space. Therefore, 
# for any given pair of axis the variables that have vectors longer than this radius make a higher contribution than average and can be interpreted with confidence.

# The top right quarter of the scaling 1 biplot shows Microcosms with the highest values for pH, C and N content of the soil
# Ther is a gradient from topleft to buttom right qarter, The gradient displays Microcosms with either high
# values for No3 (topleft) of high values of NH4 (buttomright). The buttom right Microcosms have the highest values for both cmic and NH4, 
# And also the lowes values of NO3. Keep in Mind, that No3 was multiplicated with -1 prior to analysis.

# The scaling 2 biplot shows, that variables NO3 and nH4 are highly negatively correlated, as well as Cmic and No3. N, C and pH have nearly othogonal arrows indicating 
# a correlation close to 0. 
#####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Another way of plotting ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

soil.pca <- prcomp(soil, scale=TRUE)

fit <- hclust(dist(soil.pca$x[,1:2]), method="complete") 
plot(fit)
groups <- cutree(fit, k=4) 
plot(groups)

plotPCA(soil.pca$x[,1:3],4)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# PCA of bio variables ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bio[is.na(bio)] <- 0
bio.pca <- rda(bio, scale=TRUE)
cleanplot.pca(bio.pca, ahead=0)


library(FactoMineR)
PCA(bio)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Bio and soil variables together
phyto[complete.cases(phyto$Ltbm),!names(phyto) %in% "col"] # subset only L. terrestris biomass without NAs
phyto[complete.cases(phyto$col),!names(phyto) %in% "Ltbm"] # subset only collembolan counts without NAs
phyto[is.na(phyto)] <- 0 # All NAs get 0, in particular this refers to abundance and weights of soil organisms [No weight = zero weight????]


# variation in the samples according to variables 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Lineplot for each individual samples and its according measurement
phyto <- phyto[!names(phyto) %in% c("Ltbm","col")] # exclude Lt biomass and collembolans
phyto.melt <- melt(phyto, id.vars=1:7)
ggplot(phyto.melt, aes(y=value, x=variable)) + geom_line(aes(group=ID, col=treat))

# Lineplot for each individual sample and its according centered and variance stabilized measurement
#phyto <- phyto[,-grep("ID", names(phyto))]
phyto.z <- data.frame(cbind(phyto[,1:7],scale(phyto[,8:15])))
phyto.melt <- melt(phyto.z, id.vars=1:7)
ggplot(phyto.melt, aes(y=value, x=variable)) + geom_line(aes(group=ID, col=treat))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# PCA of a merge of bio and soil variables ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phyto <- phyto[,8:15]

phyto[is.na(phyto)] <- 0
phyto.pca <- vegan::rda(phyto, scale=TRUE)  #, scale=TRUE)

cleanplot.pca(phyto.pca)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Combining Clustering and Ordination Results ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Clustering the objects using environmental data: Euclidean distance after standardization
# the variables, followed by Ward clustering


phyto.ward <- hclust(dist(scale(phyto)), "ward.D") # partition in groups with smallest variance

plot(phyto.ward)

# Cut the dendrogram to yield 4 groups
gr <- cutree(phyto.ward, k=5)
grl <- levels(factor(gr))

# Get the site scores, scaling 1
sit.sc1 <- scores(phyto.pca, display="wa", scaling=1)

# Plot the sites with cluster symbols and colours (scaling 1)
par(mfrow=c(1,1))
phyto.p1 <- plot(phyto.pca, display="wa", scaling=1, type="n", main="PCA correlation + clusters")
abline(v=0, lty="dotted")
abline(h=0, lty="dotted")

for(i in 1:length(grl)) {
        points(sit.sc1[gr==i,], pch=(14+i), cex=2, col=i+1)
}

text(sit.sc1, row.names(phyto), cex=0.7, pos=3)

# add the dendrogram
ordicluster(phyto.p1, phyto.ward, col="dark grey")
legend(locator(1), paste("Group", c(1:length(grl))), pch=14+c(1:length(grl)), col=14+c(1:length(grl)), pt.cex=2)

str(groups)
groups$cosm <- factor(groups$cosm)

phyto.rda <- rda(phyto ~ ., scale=TRUE,groups[,-c(1,5,6,7)])
summary(phyto.rda)
coef(phyto.rda)

# Unadjusted R^2 retrieved from the rda result
(R2 <- RsquareAdj(phyto.rda)$r.squared)
(R2 <- RsquareAdj(phyto.rda)$adj.r.squared)


# Scaling 1: distance triplot
plot(phyto.rda, scaling=1, main="Distance Triplot RDA")
phyto.sc <- scores(phyto.rda, choices=1:2, scaling=1, display="sp")
arrows(0,0, phyto.sc[,1], phyto.sc[,2], length=0, lty=1, col="red")

# Scaling 2: correlation triplot
plot(phyto.rda, main="Correlation Triplot RDA")
phyto2.sc <- scores(phyto.rda, choices=1:2, display="sp")
arrows(0,0, phyto2.sc[,1], phyto2.sc[,2], length=0, lty=1, col="red")


anova.cca(phyto.rda, step=1000)
anova.cca(phyto.rda, by="axis", step=1000)













