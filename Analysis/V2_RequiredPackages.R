###############
# V2 Phytocosms
# Required packages
# Quentin Schorpp
# 02.06.2015
##############


library(ggplot2)
library(reshape2)
library(plyr)
library(MethComp)
#library(agricolae)
library(foreign)
library(multcomp)
library(gridExtra)
library(extrafont)
library(effects)
library(MASS)
library(vegan)
library(mvnormtest)
library(splitstackshape)
library(FactoMineR)
library(rgl)


#### Function for Standard error
se <- function(x) sqrt(var(x)/length(x))

##### ggplot2 - Theme
mytheme = 
        theme_bw() + 
        theme(strip.background = element_rect(color = "light grey", fill="black", size=0.1),
              strip.text.x = element_text(size=8,  colour="white", face="italic"),
              strip.text.y = element_text(size=8,  colour="white", face="italic"),
              axis.text.x = element_text(size=7),
              axis.title.x = element_text(size=8,face="bold", family="Times New Roman"),
              axis.text.y = element_text(size=7),
              axis.title.y = element_text(size=8,face="bold", family="Times New Roman"),
              axis.line = element_line(size=0.25),
              axis.ticks = element_line(size=0.25),
              plot.title = element_text(size=11,face="bold", family="Times New Roman"),
              panel.margin = unit(0, "lines"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour="black", size=0.2, fill=NA))

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
