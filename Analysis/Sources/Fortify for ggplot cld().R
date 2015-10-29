
# fortify is part of the ggplot2 package and 
# turns some output/objects into data tables

if (require("multcomp")) {
  amod <- aov(Cmic ~ treat, data = cmic2)
  wht <- glht(amod, linfct = mcp(treat = "Tukey"))
  
  fortify(wht)
  ggplot(wht, aes(lhs, estimate)) + geom_point()
  
  CI <- confint(wht)
  fortify(CI)
  ggplot(CI, aes(lhs, estimate, ymin = lwr, ymax = upr)) +
    geom_pointrange()
  
  fortify(summary(wht))
  ggplot(mapping = aes(lhs, estimate)) +
    geom_linerange(aes(ymin = lwr, ymax = upr), data = CI) +
    geom_point(aes(size = p), data = summary(wht)) +
    scale_size(trans = "reverse")
  
  cld <- cld(wht)
  fortify(cld)
}