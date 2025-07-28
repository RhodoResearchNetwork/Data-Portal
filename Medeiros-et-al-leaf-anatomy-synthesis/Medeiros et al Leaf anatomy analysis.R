#Analysis and figures for Medeiros et al Leaf anatomy AOB Plants

#Table 1.Phylogenetic Canonical Correspondence Analysis.####
library(phytools)
speciesdata <-read.csv("Leaf anatomy PGLS and CCA data.csv", header = TRUE)

#get grand mean data
findspmeans<-grepl("GrandMean", speciesdata$DataType)
spmeans<-subset(speciesdata, findspmeans)
row.names(spmeans)<-spmeans$Species

#load & check the phylogeny
tree <- read.nexus("ML_constraint.nex")
plot(tree, edge.width = 2)
is.rooted(tree)
is.ultrametric(tree)#returns NO

#make the tree ultrametric, then write to file to use in analysis
tree = chronos(tree, lambda = 0, model = "correlated")
write.tree(tree,"utree.tre")

#load ultrametric tree
tree<-read.tree("utree.tre")
plotTree(tree)
is.rooted(tree)
is.ultrametric(tree)

#analysis of leaf economics, anatomy, climate and species traits
#remove columns with gas exchange data 
spmeans<-spmeans[-c(10:13)]

#makes sure only the species you want are on the tree
spmeans.tree<-drop.tip(tree, setdiff(tree$tip.label, spmeans$Species)) 
plotTree(spmeans.tree)

#create data matrices for CCA analysis
anatomyonly<-spmeans[c(5:9)]
anatomyonly<-data.matrix(anatomyonly, rownames.force=NA)

economicsonly<-spmeans[c(10:12)]
economicsonly<-data.matrix(economicsonly, rownames.force=NA)

leavesonly<-spmeans[c(5:12)]
leavesonly<-data.matrix(leavesonly, rownames.force=NA)

sptraits<-spmeans[c(13:20)]
sptraits<-data.matrix(sptraits, rownames.force=NA)

clim<-spmeans[c(21,22,24:29)]
clim<-data.matrix(clim, rownames.force=NA)

#pCCA analysis leaf economics versus climate
econclim<-phyl.cca(spmeans.tree, economicsonly, clim)

#pCCA analysis leaf anatomy versus leaf economics
anatecon<-phyl.cca(spmeans.tree, anatomyonly, economicsonly)

#pCCA analysis leaf anatomy versus climate
#remove Twet so that there are only 12 variables total across both matrices
clim2<-clim[,-c(5)]
anatclim<-phyl.cca(spmeans.tree, anatomyonly, clim2)

#pCCA analysis leaf anatomy versus species traits
#remove airspace so that there are only 12 variables total across both matrices 
anatomyonly2<-anatomyonly[,(c(1,2,4,5))]
anatomyonly2<-scale(anatomyonly2)
sptraits2<-scale(sptraits)
anatsp<-phyl.cca(spmeans.tree, sptraits2, anatomyonly2)

#pCCA analysis leaf economics versus species traits
economicsonly2<-scale(economicsonly)
econsp<-phyl.cca(spmeans.tree, economicsonly2, sptraits2, lambda=0,fixed=TRUE)

##pCCA with Gas Exchange data
#find data for Holden arboretum (HA) only
findHAmeans<-grepl("HA", speciesdata$DataType)
HAmeans<-subset(speciesdata, findHAmeans)
row.names(HAmeans)<-HAmeans$Species

#create matrices for CCA
anatomyHA<-HAmeans[c(5:9)]
anatomyHA<-data.matrix(anatomyHA, rownames.force=NA)

econHA<-HAmeans[c(14:16)]
econHA<-data.matrix(econHA, rownames.force=NA)

leavesonlyHA<-HAmeans[c(5:9,14:16)]
leavesonlyHA<-data.matrix(leavesonlyHA, rownames.force=NA)

GE<-HAmeans[c(10:13)]
GE<-data.matrix(GE, rownames.force=NA)

#pCCA gas exchange versus climate
clim3<-scale(clim)
GEclim<-phyl.cca(spmeans.tree, GE, clim3)

#pCCA gas exchange versus species traits
GEsp<-phyl.cca(spmeans.tree, GE, sptraits)

#Table 2. Multiple comparisons of species means####
spdata <-read.csv("plant means for multicomp species.csv", header = TRUE)

library(multcomp)

spdata$Species<-as.factor(spdata$Species)
amod<-aov(TOTEP~Species, data=spdata)
summary(amod)
summary(glht(amod, linfct=mcp(Species="Tukey")))

amod2<-aov(AIR~Species, data=spdata)
summary(amod2)
summary(glht(amod2, linfct=mcp(Species="Tukey")))

amod3<-aov(MRCD~Species, data=spdata)
summary(amod3)
summary(glht(amod3, linfct=mcp(Species="Tukey")))

amod4<-aov(MRXA~Species, data=spdata)
summary(amod4)
summary(glht(amod4, linfct=mcp(Species="Tukey")))

amod5<-aov(A~Species, data=spdata)
summary(amod5)#no sig diffs across species

amod6<-aov(E~Species, data=spdata)
summary(amod6) #no sig diffs across species

amod7<-aov(gs~Species, data=spdata)
summary(amod7) #no sig diffs across species

amod8<-aov(WUE~Species, data=spdata)
summary(amod8) #no sig diffs across species

#Supplementary Tables. PGLS Model selection####
library(caper)
library(ape)
library(geiger)
library(nlme)
library(picante)
#create a comparative data object
anatomycompare <- comparative.data(spmeans.tree, spmeans, Species, vcv=TRUE)

#Table S1. Leaf anatomical traits versus leaf economics traits.####
modal <- pgls(spmeans$TOTEP~spmeans$SLA+spmeans$CN+spmeans$Habit, anatomycompare, lambda = 'ML')
summary(modal)
shapiro.test(modal$residuals)#ok
plot(modal)#ok

modal1 <- pgls(spmeans$PAL~spmeans$SLA+spmeans$CN+spmeans$Habit, anatomycompare, lambda = 'ML')
summary(modal1)
shapiro.test(modal1$residuals)#ok
plot(modal1)#ok

modal2 <- pgls(spmeans$AIR~spmeans$SLA+spmeans$CN+spmeans$Habit, anatomycompare, lambda = 'ML')
summary(modal2)
shapiro.test(modal2$residuals)#ok
plot(modal2)#ok

modal3 <- pgls(spmeans$MRCD~spmeans$SLA+spmeans$CN+spmeans$Habit, anatomycompare, lambda = 'ML')
summary(modal3)
shapiro.test(modal3$residuals)#ok
plot(modal3)#ok

modal4 <- pgls(spmeans$MRXA~spmeans$SLA+spmeans$CN+spmeans$Habit, anatomycompare, lambda = 'ML')
summary(modal4)
shapiro.test(modal4$residuals)#ok
plot(modal4)#ok

#Table S2. Climate versus leaf economics.####
modSLAa <- pgls(spmeans$SLA~spmeans$MAT, anatomycompare, lambda = 'ML')
summary(modSLAa)
shapiro.test(modSLAa$residuals)#ok
plot(modSLAa)#ok

modlife <- pgls(spmeans$LeafLifespan~spmeans$MAT, anatomycompare, lambda = 'ML')
summary(modlife)
shapiro.test(modlife$residuals)#ok
plot(modlife)#ok

modlife1 <- pgls(spmeans$CN~spmeans$Tseason, anatomycompare, lambda = 'ML')
summary(modlife1)
shapiro.test(modlife1$residuals)#ok
plot(modlife1)#ok

#Table S3. Species traits versus leaf anatomy.####
mod1ax <- pgls(spmeans$SLA~spmeans$SRL+spmeans$StemVesD+spmeans$CVSLA+spmeans$LeafNo+spmeans$BranchLeafArea+spmeans$Huber, anatomycompare, lambda = 'ML')
summary(mod1ax)
shapiro.test(mod1ax$residuals)#ok
plot(mod1ax)#ok

mod1bx <- pgls(spmeans$LeafLifespan~spmeans$SRL+spmeans$FOD+spmeans$StemVesD+spmeans$CVSLA+spmeans$LeafNo+spmeans$BranchLeafArea, anatomycompare, lambda = 'ML')
summary(mod1bx)
shapiro.test(mod1bx$residuals)#ok
plot(mod1bx)#ok

mod1cx <- pgls(spmeans$CN~spmeans$SRL+spmeans$FOD+spmeans$StemVesD+spmeans$BranchLeafArea+spmeans$Huber, anatomycompare, lambda = 'ML')
summary(mod1cx)
shapiro.test(mod1cx$residuals)#ok
plot(mod1cx)#ok

mod1a <- pgls(spmeans$TOTEP~spmeans$SRL+spmeans$FOD+spmeans$StemVesD+spmeans$CVSLA+spmeans$LeafSize+spmeans$LeafNo+spmeans$BranchLeafArea+spmeans$Huber, anatomycompare, lambda = 'ML')
summary(mod1a)
shapiro.test(mod1a$residuals)#ok
plot(mod1a)#ok

mod1b <- pgls(spmeans$PAL~spmeans$CVSLA, anatomycompare, lambda = 'ML')
summary(mod1b)
shapiro.test(mod1b$residuals)#ok
plot(mod1b)#ok

mod1c <- pgls(spmeans$AIR~spmeans$BranchLeafArea, anatomycompare, lambda = 'ML')
summary(mod1c)
shapiro.test(mod1c$residuals)#ok
plot(mod1c)#ok

mod1d <- pgls(spmeans$MRCD~spmeans$StemVesD, anatomycompare, lambda = 'ML')
summary(mod1d)
shapiro.test(mod1d$residuals)#ok
plot(mod1d)#ok

mod1e <- pgls(spmeans$MRXA~spmeans$LeafSize, anatomycompare, lambda = 'ML')
summary(mod1e)
shapiro.test(mod1e$residuals)#ok
plot(mod1e)#ok


#Table S4. Gas exchange versus climate.####
#find data for HA only
findHAmeans<-grepl("HA", speciesdata$DataType)
HAmeans<-subset(speciesdata, findHAmeans)
row.names(HAmeans)<-HAmeans$Species

#makes sure only the species you want are on the tree
HAmeans.tree<-drop.tip(tree, setdiff(tree$tip.label, HAmeans$Species)) 
plotTree(HAmeans.tree)
name.check(HAmeans.tree,HAmeans)

#create a comparative object
HAcompare <- comparative.data(HAmeans.tree, HAmeans, Species, vcv=TRUE)

#model selection for climate variables that explain gas exchange 
mod3c <- pgls(HAmeans$A~spmeans$Tseason, HAcompare, lambda = 'ML')
summary(mod3c)
shapiro.test(mod3c$residuals)#ok
plot(mod3c)#ok

mod4c<- pgls(HAmeans$E~spmeans$MAT+spmeans$Tseason+spmeans$Pseason+spmeans$TAR+spmeans$Twet, HAcompare, lambda = 'ML')
summary(mod4c)
shapiro.test(mod4c$residuals)#ok
plot(mod4c)#ok

mod5c <- pgls(HAmeans$gs~spmeans$Tseason, HAcompare, lambda = 'ML')
summary(mod5c)
shapiro.test(mod5c$residuals)#ok
plot(mod5c)#ok


mod5cd <- pgls(HAmeans$WUE~spmeans$Pwarm, HAcompare, lambda = 'ML')
summary(mod5cd)
shapiro.test(mod5cd$residuals)#ok
plot(mod5cd)#ok

