library(rgl)
library(scatterplot3d)
library(vegan)
library(vegan3d)

OTUs_Bulk<-read.table("OTUs_Bulk.txt", header=TRUE)
Phylum_Bulk<-read.table("Phylum_Bulk.txt", header=TRUE)
rda.vasc.0 <- rda(OTUs_Bulk ~ 1, data = Phylum_Bulk) # model containing only species matrix and intercept
rda.vasc.all <- rda(OTUs_Bulk ~ ., data = Phylum_Bulk) # model including all variables from matrix chem (the dot after tilda (~) means ALL!)
Bulk_Phylum<-ordistep (rda.vasc.0, scope = formula (rda.vasc.all), direction = 'both')
save(Bulk_Phylum, file="Bulk_Phylum.txt")
Bulk_Phylum.summary<-summary(Bulk_Phylum)
save(Bulk_Phylum.summary, file="Bulk_Phylum.summary.txt")
Bulk_Phylum.anova<-anova(Bulk_Phylum)
save(Bulk_Phylum.anova, file="Bulk_Phylum.anova.txt")
Bulk_Phylum.anova_term<-anova(Bulk_Phylum, by = "term", step=1000)    # Permutation test for rda under reduced model  Terms added sequentially (first to last)
save(Bulk_Phylum.anova_term, file="Bulk_Phylum.anova_term.txt")
Bulk_Phylum.anova_margin<-anova(Bulk_Phylum, by = "margin", step=1000)  # Permutation test for rda under reduced model  Marginal effects of terms
save(Bulk_Phylum.anova_margin, file="Bulk_Phylum.anova_margin.txt")
Bulk_Phylum.anova_axis<-anova(Bulk_Phylum, by = "axis", step=1000)
save(Bulk_Phylum.anova_axis, file="Bulk_Phylum.anova_axis.txt")

load("Bulk_Phylum.txt")
load("Bulk_Phylum.summary.txt")
load("Bulk_Phylum.anova.txt")
load("Bulk_Phylum.anova_term.txt")
load("Bulk_Phylum.anova_margin.txt")
load("Bulk_Phylum.anova_axis.txt")

bg3d(color="black")
ordirgl(Fractions_Class, display="sites", choices = 1:3, type = "p", ax.col = "white", arr.col = "yellow", text, col="black")
with(Class_Fractions, texts3d(scores(Fractions_Class, choices=1:3, display="sites"), texts=c("*", "^", "v", "#")[Fraction_], 
                              col=c("green", "deepskyblue", "red")[Soil_], font="bold", cex=c(2,2,1.25,1.25)[Fraction_]))
axes3d(edges = "bbox", labels = TRUE, tick = TRUE, box=FALSE, expand = 1.03, col="white")


bg3d(color="white")
ordirgl(Fr_Ge, display="sites", choices = 1:3, type = "p", ax.col = "grey", arr.col = "black", text, 
        col="transparent")
with(Class_Fractions, texts3d(scores(Fr_Ge, choices=1:3, display="sites"), texts=c("*", "^", "v", "#")[Fraction_], 
                              col=c("darkgreen", "blue", "red")[Soil_], font="bold", cex=c(2,2,1.25,1.25)[Fraction_]))
axes3d(edges = "bbox", labels = TRUE, tick = TRUE, box=FALSE, expand = 1.03, col="black")
texts3d(x=1.6, y=0, z=0, texts="15.9 %", col="grey")
texts3d(x=0, y=1.9, z=0, texts="6.8 %", col="grey")
texts3d(x=0, y=0, z=3.9, texts="4.1 %", col="grey")
with(Class_Fractions, legend3d("topleft", c("63-2000 µm", "20-    63 µm", "   2-    20 µm", "          < 2 µm", "UNF", "1.0 NPK", "1.5 AM"), 
                               pch=c("#","^","v","*","O","O","O"), pt.cex=c(1,1.5,1,1.5, 1, 1, 1), 
                               col=c("black", "black", "black","black", "red", "blue", "darkgreen")))


play3d(spin3d(axis=c(1,1,1), rpm=5), duration=Inf)