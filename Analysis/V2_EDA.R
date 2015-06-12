#################
# V2 Phytocosms
# Exploratory Data Analysis
# Quentin Schorpp
# 02.06.2015
#################

# Leaching water

str(leaching)

leaching.melt <- melt(leaching, id.vars=c(1:7))
leaching.melt <-concat.split.multiple(leaching.melt, "variable", ".")

leaching2 <- leaching.melt[!(leaching.melt$variable_2==1),]

ggplot(leaching2, aes(x=variable_2, y=value, col=treat), labeller=label_parsed) +
        geom_point(lwd=2) +
        ylab(expression(paste("?????"))) +
        xlab("sampling dates") +
        ggtitle(expression(paste("Percolation Water"))) +
        facet_grid(variable_1+plant~soil, labeller=label_parsed) + 
        mytheme + theme(legend.position="bottom")