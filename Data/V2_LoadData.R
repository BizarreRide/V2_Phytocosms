###############
# V2 Phytocosms
# Load Data
# Quentin Schorpp
# 02.06.2015
##############


# load datasets ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index<- read.delim("Data/V2_index.txt")
Lt.bm.loam <- read.delim("Data/V2_Lt_bm_loam.txt")
Lt.bm <- read.delim("Data/V2_Lt_bm.txt")
col.pop <- read.delim("Data/V2_col.txt")
pH <- read.delim("Data/V2_pH.txt")
CN <- read.delim("Data/V2_CN.txt")
roots <- read.delim("Data/V2_roots.txt")
shoots <- read.delim("Data/V2_shoots.txt")
cmic <- read.delim("Data/V2_cmic.txt")
leaching <- read.delim("Data/V2_Leaching.txt")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Configure "Index"
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#index$plant <- revalue(index$plant, c(Sp ="italic(S.~perfoliatum)",Zm="italic(Z.~mays)"))
index$soil <- revalue(index$soil, c(loam ="Loam",sand="Sand"))
index$Lt <- as.factor(index$Lt)
index$Fc <- as.factor(index$Fc)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Process Microbial Biomass
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmic2 <- cmic[!(cmic$plant=="Wp"),] # exclude "Without Plant"
cmic2$ID1  <- as.factor(rep(c(1:15),16)) #  Id for ???
cmic2$treat <- factor(cmic2$treat, levels=c("Lt","Fc","Mix","Control")) # reorder factor levels

# Boxplot for Cmic
ggplot(cmic2, aes(x= ID1, y=Cmic),group=cosm, labeller=label_parsed) +
        geom_bar(fill="light blue", col="black", lwd=0.2, stat="identity") +
        ylab(expression(paste("µg-C /g soil"))) +
        xlab("Pseudo-Replicates") +
        ggtitle(expression(paste("Microbial Biomass - Outlier search"))) +
        facet_grid(soil~plant+treat, labeller=label_parsed) + 
        mytheme

# Boxplot for Airflow
ggplot(cmic2, aes(x= ID1, y=Airflow),group=cosm, labeller=label_parsed) +
        geom_bar(fill="light blue", col="black", lwd=0.2, stat="identity") +
        ylab(expression(paste("µg-C /g soil"))) +
        xlab("Pseudo-Replicates") +
        coord_cartesian(ylim=c(130,190)) +
        ggtitle(expression(paste("Microbial Biomass - Outlier search"))) +
        facet_grid(soil~plant+treat, labeller=label_parsed) + 
        mytheme
# How will Airflow be used in data processing? or is it even incorporated yet?

# Eliminate Outliers
a <- which(cmic2$Cmic<200)
a <- c(a,with(cmic2[(cmic2$treat=="Lt" & cmic2$plant=="Zm" & cmic2$soil=="Loam"),], which(cmic2$Cmic==min(Cmic))))
cmic2 <- cmic2[-a,]; rm(a)
cmic2 <- aggregate(Cmic~cosm+soil+plant+treat, cmic2, mean)
a <- with(cmic[(cmic$plant=="Wp"),], tapply(Cmic, list(soil,treat), mean))[1,6]
b <- with(cmic[(cmic$plant=="Wp"),], tapply(Cmic, list(soil,treat), mean))[2,6]
cmic2 <- cmic2[order(cmic2$soil,cmic2$cosm),]
cmic2$Cmic.start <- rep(c(a,b), each=40)
        
ggplot(cmic2, aes(x=treat, y=Cmic))
ggplot(cmic2, aes(x=treat, y=Cmic), labeller=label_parsed) +
        geom_bar(fill="light blue", col="black", lwd=0.2, stat="identity") +
        ylab(expression(paste("µg-C /g soil"))) +
        xlab("Pseudo-Replicates") +
        ggtitle(expression(paste("Microbial Biomass - Outlier search"))) +
        facet_grid(soil~plant, labeller=label_parsed) + 
        mytheme
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Create Full Frame with all variables ####
# Variables should resemble the final state or the value with most clear influence of factor index variables, i.e. "change"
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
col.pop[col.pop$col==0,]$col <- NA
Lt.bm2 <- merge(Lt.bm[,-c(2:5)], col.pop, by="ID", all.y =TRUE, incomparables=0)

full.frame <- cbind(index,
      ph.before=pH$pH1,
      ph.after=pH$pH5,
      CN[,-c(1:5)],
      roots[,-c(1:5)],
      shoots[,-c(1:5)],
      col=col.pop$col,
      Lt.bm2[,2:4])

full.frame <- full.frame[-c(41,82),]

full.frame <- cbind(full.frame, Cmic=cmic2$Cmic, leaching[,-c(1:7)])


full.frame$root <- full.frame$root.after - full.frame$root.before
full.frame$shoot <- full.frame$shoot.after - full.frame$shoot.before
full.frame$Ltbm <- full.frame$Lt.bm.diff
full.frame$col <- full.frame$col
full.frame$Cmic <- cmic2$Cmic - cmic2$Cmic.start # should "start values be subtracted?
full.frame$CN <- full.frame$CN.after - full.frame$CN.before
full.frame$N <- full.frame$N.after - full.frame$N.before
full.frame$C <- full.frame$C.after - full.frame$C.before
full.frame$NH4 <- full.frame$NH4.5 - full.frame$NH4.1
full.frame$NO3 <- full.frame$NO3.5 - full.frame$NO3.1
full.frame$pH <- full.frame$ph.after - full.frame$ph.before

full.frame <- cbind(full.frame[,1:7], full.frame[,c("root","shoot", "Ltbm", "col", "Cmic", "CN", "C", "N", "NH4", "NO3", "pH")])
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






