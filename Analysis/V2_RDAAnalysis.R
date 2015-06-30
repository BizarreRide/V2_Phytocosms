###############
# V2 phytocosms
# RDA Analysis
# Quentin Schorpp
# 02.06.2015
##############

# Description:
# This is an attempt for mutlivariate analyses of phytocosms from the phytocosms experiment of Valentina Sandor
# All measures are calculated as Differences between beginning and end of the experiment, i.e. changes
        

# Load data and packages and Divide data frame in subsets ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Analysis/V2_RequiredPackages.R") #packages agricolae and vegan are colliding
source("Data/V2_LoadData.R")
source("Analysis/evplot.R")
source("Analysis/cleanplot.pca.R")

# Sample 11 has strange Cmic values compared to other replicates, see stacked bars in LoadData.R
# full.frame <- full.frame[-11,] # phytocosm 11 was strange?

soil <- full.frame[,c(12,14:18)] # abiotic soil measures, including Cmic
colnames(soil)[2:3] <- c("carbon","Nitrogen")
bio <- full.frame[,8:12] # biotic soil characteristics
env <- full.frame[,c(8:12,14:18)] # the same as env
groups <- full.frame[,1:7]
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Soil ####        
## Z-transformation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Here is an approach, that tries to handle negative values in some variables
# Since No3 is strictly negative, i just turn it into a positive variable
# Henc all changes are negative, but expressed as the absolute of the values
# this needs some time to think about, when looking at the PCA biplots

#soil$NO3 <- soil$NO3*(-1)

# In the next step i add the according minimum value to all variables, this leads to shift in the mesurement range
# of the minimum value, zeroes in the original measurements get that value.

#for (i in 1:6) {
#        na.action=na.omit(soil)
#        soil[,i] <- soil[,i] + abs(min(soil[,i], na.rm=TRUE))   
#}

# To achieve normal distributed variables, I calculate the log(x+1) of all values, 
# since there are still zeroes in the data

# soil.log <- log1p(soil)

# TO account for heterogeneity in the variances of the varaibles, i scale them 
# scaling before doing PCA is redundant, because rda(Scale=TRUE) can do this

#soil.z <- data.frame(scale(soil.log))

# All NAs get 0, in particular this refers to abundance and weights of soil organisms 
# [No weight = zero weight????]

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

## PCA  ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# PCA
soil.pca <- rda(soil.z, scale=FALSE)
summary(soil.pca) # scaling 2 by default
summary(soil.pca, scaling=1)

# Loads of variables:
scores(soil.pca, choices = 1:4, display = "species", scaling = 0)

### Some Explanations of PCA Terms: ####
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

### Biplots ####
eigenvalues <- soil.pca$CA$eig
evplot(eigenvalues) # One interprets only the axes whose eigenvalues are larger than the length of corresponding piece of the sticks

cleanplot.pca(soil.pca, ahead=0)

par(mfrow=c(1,2))
biplot(soil.pca, scaling=1,  main="PCA - Scaling 1")
pcacircle(soil.pca)
biplot(soil.pca, scaling=2,  main="PCA - Scaling 2")

#### Interpretation ####
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### Another way of plotting ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

soil.pca <- prcomp(soil, scale=TRUE)

fit <- hclust(dist(soil.pca$x[,1:2]), method="complete") 
plot(fit)
groups <- cutree(fit, k=4) 
plot(groups)

plotPCA(soil.pca$x[,1:3],4)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# bio variables ####
## PCA ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bio[is.na(bio)] <- 0
bio.pca <- rda(bio, scale=TRUE)
cleanplot.pca(bio.pca, ahead=0)


library(FactoMineR)
PCA(bio)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Environmental (env) without Lt biomass and collembolans ####
# This datasets omits CN ratio
env2 <- cbind(groups,env)  # merge of bio and soil variables 

# Subsetting
env.Lt <- env[complete.cases(env$Ltbm),!names(env) %in% "col"] # subset only L. terrestris biomass without NAs
env.col <- env[complete.cases(env$col),!names(env) %in% "Ltbm"] # subset only collembolan counts without NAs
env.total <- env[is.na(env)] <- 0 # All NAs get 0, in particular this refers to abundance and weights of soil organisms [No weight = zero weight????]



## variation in the samples according to variables ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Lineplot for each individual samples and its according measurement
env2 <- env2[!names(env2) %in% c("Ltbm","col")] # exclude Lt biomass and collembolans
env2.melt <- melt(env2, id.vars=1:7)
ggplot(env2.melt, aes(y=value, x=variable)) + geom_line(aes(group=ID, col=treat))

# Lineplot for each individual sample and its according centered and variance stabilized measurement
#env2 <- env2[,-grep("ID", names(env2))]
env2.z <- data.frame(cbind(env2[,1:7],scale(env2[,8:15])))
env2.melt <- melt(env2.z, id.vars=1:7)
ggplot(env2.melt, aes(y=value, x=variable)) + geom_line(aes(group=ID, col=treat))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## PCA ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

env <- env[,!names(env) %in% c("Ltbm","col")]
# env[is.na(env)] <- 0 #alternatively if Ltbm and col should stay in the dataframe
env.pca <- vegan::rda(env, scale=TRUE)  
summary(env.pca)
summary(env.pca, scaling=1)

eigenvalues <- env.pca$CA$eig
evplot(eigenvalues) # One interprets only the axes whose eigenvalues are larger than the length of corresponding piece of the sticks
# three Axes are important

cleanplot.pca(env.pca)

env.pca2 <- prcomp(env, scale=TRUE)
par(mfrow=c(1,1))
fit <- hclust(dist(env.pca2$x[,1:2]), method="complete") 
plot(fit)
groups.env <- cutree(fit, k=4) 

biplot(env.pca2)

plotPCA(env.pca2$x[,1:3],4)


# In der Piepeline
env.Lt.pca <- vegan::rda(env.Lt, scale=TRUE)  
env.col.pca <- vegan::rda(env.col, scale=TRUE)  
cleanplot.pca(env.Lt.pca)
cleanplot.pca(env.col.pca)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Combining Clustering and Ordination Results ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Clustering the objects using environmental data: Euclidean distance after standardization
# the variables, followed by Ward clustering


env.ward <- hclust(dist(scale(env)), "ward.D") # partition in groups with smallest variance

plot(env.ward)

# Cut the dendrogram to yield 4 groups
gr <- cutree(env.ward, k=5)
grl <- levels(factor(gr))

# Get the site scores, scaling 1
sit.sc1 <- scores(env.pca, display="wa", scaling=1)

# Plot the sites with cluster symbols and colours (scaling 1)
par(mfrow=c(1,1))
env.p1 <- plot(env.pca, display="wa", scaling=1, type="n", main="PCA correlation + clusters")
abline(v=0, lty="dotted")
abline(h=0, lty="dotted")

for(i in 1:length(grl)) {
        points(sit.sc1[gr==i,], pch=(14+i), cex=2, col=i+1)
}

text(sit.sc1, row.names(env), cex=0.7, pos=3)

# add the dendrogram
ordicluster(env.p1, env.ward, col="dark grey")
legend(locator(1), paste("Group", c(1:length(grl))), pch=14+c(1:length(grl)), col=14+c(1:length(grl)), pt.cex=2)

str(groups)
groups$cosm <- factor(groups$cosm)

env.rda <- rda(env ~ ., scale=TRUE,groups[,-c(1,5,6,7)])
summary(env.rda)
coef(env.rda)

# Unadjusted R^2 retrieved from the rda result
(R2 <- RsquareAdj(env.rda)$r.squared)
(R2 <- RsquareAdj(env.rda)$adj.r.squared)


# Scaling 1: distance triplot
plot(env.rda, scaling=1, main="Distance Triplot RDA")
env.sc <- scores(env.rda, choices=1:2, scaling=1, display="sp")
arrows(0,0, env.sc[,1], env.sc[,2], length=0, lty=1, col="red")

# Scaling 2: correlation triplot
plot(env.rda, main="Correlation Triplot RDA")
env2.sc <- scores(env.rda, choices=1:2, display="sp")
arrows(0,0, env2.sc[,1], env2.sc[,2], length=0, lty=1, col="red")


anova.cca(env.rda, step=1000)
anova.cca(env.rda, by="axis", step=1000)













