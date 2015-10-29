#####################
# PCA plots
# functions and internet search
# Quentin Schorpp
# 19.06.2015
####################


# from http://stackoverflow.com/questions/24282143/pca-multiplot-in-r
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

res.pca <- prcomp(data, scale=TRUE)

fit <- hclust(dist(res.pca$x[,1:2]), method="complete") 
plot(fit)
groups <- cutree(fit, k=4) 
plot(groups)


library(rgl)
plotPCA <- function(x, nGroup) {
        n <- ncol(x) 
        if(!(n %in% c(2,3))) { # check if 2d or 3d
                stop("x must have either 2 or 3 columns")
        }
        
        fit <- hclust(dist(x), method="complete") # cluster
        groups <- cutree(fit, k=nGroup)
        
        if(n == 3) { # 3d plot
                plot3d(x, col=groups, type="s", size=1, axes=F)
                axes3d(edges=c("x--", "y--", "z"), lwd=3, axes.len=2, labels=FALSE)
                grid3d("x")
                grid3d("y")
                grid3d("z")
        } else { # 2d plot
                maxes <- apply(abs(x), 2, max)
                rangeX <- c(-maxes[1], maxes[1])
                rangeY <- c(-maxes[2], maxes[2])
                plot(x, col=groups, pch=19, xlab=colnames(x)[1], ylab=colnames(x)[2], xlim=rangeX, ylim=rangeY)
                lines(c(0,0), rangeX*2)
                lines(rangeY*2, c(0,0))
        }
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




# from: http://www.r-bloggers.com/using-r-two-plots-of-principal-component-analysis/
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(reshape2)
library(ggplot2)

sample.groups <- c(rep(1, 10), rep(2, 10), rep(3, 10),
                   rep(4, 10), rep(5, 10))
variable.groups <- c(rep(1, 10), rep(2, 10), rep(3, 10),
                     rep(4, 10), rep(5, 10), rep(6, 10),
                     rep(7, 10))

data <- matrix(nrow=length(sample.groups), ncol=70)
base.data <- matrix(nrow=length(sample.groups), ncol=7)

for (j in 1:ncol(base.data)) {
        mu <- rnorm(1, 0, 4)
        sigma <- runif(1, 5, 10)
        base.data[,j] <- sample.groups*mu +
                rnorm(length(sample.groups), 0, sigma)
}

for (j in 1:ncol(data)) {
        mu <- runif(1, 0, 4)
        data[,j] <- base.data[,variable.groups[j]] +
                rnorm(length(sample.groups), mu, 10)
}


heatmap <- qplot(x=Var1, y=Var2, data=melt(cor(data)), geom="tile",
                 fill=value)

heatmap


pca <- prcomp(data, scale=T)
melted <- cbind(variable.groups, melt(pca$rotation[,1:9]))

barplot <- ggplot(data=melted) +
        geom_bar(aes(x=Var1, y=value, fill=variable.groups), stat="identity") +
        facet_wrap(~Var2)
barplot

scores <- data.frame(sample.groups, pca$x[,1:3])
pc1.2 <- qplot(x=PC1, y=PC2, data=scores, colour=factor(sample.groups)) +
        theme(legend.position="none")
pc1.3 <- qplot(x=PC1, y=PC3, data=scores, colour=factor(sample.groups)) +
        theme(legend.position="none")
pc2.3 <- qplot(x=PC2, y=PC3, data=scores, colour=factor(sample.groups)) +
        theme(legend.position="none")

pc1.2

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


