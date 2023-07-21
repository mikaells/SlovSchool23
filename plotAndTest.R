#Read in libraries
library(vegan)
library(beeswarm)

####
#***Install if you dont have packages!
####

#read in data

#first the actual bacteria
tax=read.csv(file = "w3.csv")
#then their treatment group
groups=read.csv(file = "w3_groups.csv")

#####
#
#1) Alpha diversity
#
#####

#
SHANNON=diversity(x = tax,index = "shannon")
SIMPSON=diversity(x = tax,index = "shannon")
RICH=estimateR(tax)[1,]

#plot
beeswarm(CHAO~groups$x)
#test
t.test(SHANNON~groups$x)

####
#***repeat plot and test with the two other diversity metrics!
####

#####
#
#2) Beta diversity
#
#####

#visualization

#calculating nMDS
nmds=metaMDS(tax)

#plotting scores
#color by group
plot(scores(nmds, display = "sites"), col=factor(groups$x), pch=16)

#test for multivariate difference
adonis2(tax~groups2)

#####
#
#3) Investigating individual genera
#
#####

ENV=envfit(ord=nmds, env=tax)
plot(scores(nmds, display = "sites"), col=factor(groups$x),, pch=16)

names(which(ENV$vectors$r>0.5))

#plot
beeswarm(tax$Bacteria_Bacteroidota_Bacteroidia_Bacteroidales_Prevotellaceae_Prevotella~groups$x)
#test
t.test(tax$Bacteria_Bacteroidota_Bacteroidia_Bacteroidales_Prevotellaceae_Prevotella~groups$x)

####
#***repeat plot and test with other interesting bacteria!
####

