###Title: Analysis for unguarded X hypothesis
#Author: Suz Everingham
#Date: 8/8/19

##Contents:
#lines : reading data, creating log response ratio
#lines : running analysis of the effect of data type (mean/max/median) on lnRR and plotting results
#lines : cleaning zoe's data to just have means (where species have multiple data points)
#lines : running analysis of the effect of captive vs wild on lnRR and plotting results
#lines : testing for phylogenetic signal
#lines : running meta-analysis both with/w-out phylo (same answer)
#lines : plotting big graph of all species, lnRR on the phylo tree :)

library(metafor)
library(dplyr)
library(ggplot2)
library(ape)
library(phytools)
library(rotl)
library(tidyverse)
library(grid)
library(phylobase)
library(picante)
library(multcomp)
library(ggtree)
library(ggimage)

#read data
zoedata <- read.csv("lifespan.csv")

#create lnOR (log odds ratio)

zoedata <- zoedata %>%
  mutate(lnRR= log(homogametic.lifespan/heterogametic.lifespan)) %>% ##creating a lnRR variable
  transform(Species=as.character(Species))

# run an analysis to determine whether there is any difference
#between mean of max and max of max lifespan data in lnRR
meanormaxormedian <- lm(lnRR ~ data.type.age, data=zoedata)
summary(meanormaxormedian)
cont <- glht(meanormaxormedian, linfct = mcp(data.type.age = "Tukey"))
summary(cont)

#plot
zoedatatypeplot <- ggplot(aes(x=data.type.age, y=lnRR), data = zoedata) +
  geom_violin(fill="lightskyblue1") +
  stat_summary(fun.y=mean, geom="point", size=1, color="salmon") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="royalblue4", size=1) +
  theme_classic(base_size = 15) +
  xlab("Data Type") +
  ylab("ln(homogametic lifespan/heterogametic lifespan)")

plot(zoedatatypeplot)

##now try with just species where there is both mean and max data

zoemult <- zoedata %>%
  group_by(Species) %>% ### this whole function is subsetting only data that 
  filter(n()>1)         ##### has replicated spp

meanmaxmed <- lm(lnRR ~ data.type.age, data = zoemult)
summary(meanmaxmed)

##no difference in the species where there is multiple reps for lifespan data
##will just go for mean (because there is more data there already)
## cut out species with just mean vals and insert back into zoe data
zoemult <- filter(zoemult, data.type.age == 'mean') #selecting species that are only mean data in multiples


zoedata <- zoedata %>%
  group_by(Species) %>%
  filter(n()==1) ##removing any species that has any multiple data in zoes data

zoedata <- bind_rows(zoedata, zoemult)
zoedata <- zoedata[!is.na(zoedata$heterogametic.lifespan),]
#putting the selected means from the multiple data list with zoes list
##this dataset now gives me one point for each species, and where species used to have multiple points,
#they now only have one- and it is always mean.


###running an analysis for captive vs. wild 
captivevswild <- lm(lnRR ~ Population.source, data = zoedata)
summary(captivevswild)
captiveswildcont <- glht(captivevswild, linfct = mcp(Population.source = "Tukey"))
summary(captiveswildcont)


zoecaptiveswildplot <- ggplot(aes(x=Population.source, y=lnRR), data=zoedata) +
  geom_violin(fill="lightskyblue1") +
  stat_summary(fun.y=mean, geom="point", size=1, color="salmon") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="royalblue4", size=1) +
  theme_classic(base_size = 15) +
  xlab("") +
  ylab("ln(homogametic lifespan/heterogametic lifespan)")
  
plot(zoecaptiveswildplot)

###running a phylogenetic analysis
tree <- read.newick("phyloT_generated_tree.txt")
tree <- compute.brlen(tree, power = 1)
cor <- vcv(tree, cor = T)
plot(tree, cex = 0.6)
zoedata$Species <- as.character(zoedata$Species)
zoedata$lnRR <- as.numeric(zoedata$lnRR)


zoedata$Species_ <- gsub(" ", "_", zoedata$Species)
zoedata <- rename(zoedata, "tree$tip.label"=Species_)
zoedata <- right_join(zoedata, as.data.frame(tree$tip.label), by = "tree$tip.label")
zoedata <- rename(zoedata, Species_="tree$tip.label")
species.vector <- as.vector(zoedata$Species_)
zoedata[, "phylo"] <- species.vector
lnvector <- as.vector(zoedata$lnRR)
names(lnvector) <- zoedata$Species_

tree <- collapse.singles(tree)

lambda <- phylosig(tree, lnvector, method="lambda", test=TRUE)
print(lambda)

lambda2 <- phylosig(tree, lnvector, method="K", test=TRUE)
print(lambda2)

phylo4d <- phylo4d(tree, lnvector)


#run meta-analysis on all data

Xsexanalysis <- rma.mv(yi=lnRR, 
                       V = 1,
                       random = ~ 1|phylo,
                       R = list(phylo=cor),
                       method = "REML",
                       data = zoedata)
summary(Xsexanalysis)

# and without phylogeny

Xsexanalysis2 <- rma.uni(yi=lnRR,
                         vi=1,
                         method = "REML",
                         data=zoedata)
summary(Xsexanalysis2)

##or perform a simple t-test to see if mean is signficantly different from 0
t.test(zoedata$lnRR)
#both v.significant :)

####get magnitude of difference###
#overall mean

zoedata <- mutate(zoedata, dimorphism = (((homogametic.lifespan - heterogametic.lifespan)/homogametic.lifespan)*100)) #gives percent of lifespan difference between homogametic and heterogametic
mean(zoedata$dimorphism)
### average difference = 10.3% therefore on average, homogametic organism lives 11.5% longer than hetero

#also for female vs male hetero
aggregate(zoedata[, 5:21], list(zoedata$Sex.Determination), mean)

##female heterogametic = 4.1 % and male heterogametic = 16.7%

#also for class 
aggregate(zoedata[, 1:21], list(zoedata$Class), mean)
#Actinopterygii = 19.1%
#Amphibia = 2.1%
#Arachnida = -10.1%
#Aves = 5.1%
#Chondrichthyes = -298.7% ?? very weird but probably due to two species only
#Insecta = 16.1%
#Mammalia = 12.3%
#Reptilia = 14.7%

#get some data on numbers in total
classdata <- as.data.frame(summary(zoedata$Class))
orderdata <- as.data.frame(summary(zoedata$Order))
familydata <- as.data.frame(summary(zoedata$Family))

#is there a stastically significant difference between male or female heterogametic lifespan

maleorfemale <- lm(lnRR ~ Sex.Determination, data=zoedata)
summary(maleorfemale)
plot(maleorfemale)
malefemalecont <- glht(maleorfemale, linfct = mcp(Sex.Determination = "Tukey"))
summary(malefemalecont)

zoemalefemaleplot <- ggplot(aes(x=Sex.Determination, y=lnRR), data=zoedata) +
  geom_violin(fill="lightskyblue1") +
  stat_summary(fun.y=mean, geom="point", size=1, color="salmon") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="royalblue4", size=1) +
  theme_classic(base_size = 15) +
  xlab("") +
  ylab("ln(homogametic lifespan/heterogametic lifespan)")

plot(zoemalefemaleplot)


###plotting figures###


#tree for appendix
tiff("tree.tiff", width=1000, height=8000, res=200)
ggtree(tree)
dev.off()


##tree for manuscript
plottree <- ggtree(tree)
plot(plottree)
d1 <- data.frame(id=tree$tip.label, lnRR=zoedata$lnRR, sexdet=zoedata$Sex.Determination)
p2 <- facet_plot(plottree, panel="dot", data=d1, geom=geom_point, aes(x=lnRR, group=sexdet, color=sexdet, shape=sexdet), size=2.5) +
  scale_color_manual(values=c("firebrick1", "lightskyblue")) +
  scale_shape_manual(values=c(18,16)) +
  theme(axis.title.x=element_text(size=10),
        axis.text.x=element_text(size=10)) +
  geom_vline(xintercept=0) +
  theme(legend.position = c(0.85,0.95),
        legend.title=element_blank(),
        legend.text=element_blank())

tiff("plot1.tiff", height=3000, width=2000, res=300)
plot(p2)
dev.off()

ggplot(zoedata, aes(x=lnRR, y=Species, color=Sex.Determination)) + 
  geom_point(shape=18, size=8) +
  scale_color_manual(values=c("plum1", "forestgreen")) +
  theme_classic() + 
  theme(axis.text.x=element_text(size=30),
        axis.text.y=element_text(size=20),
       axis.title=element_text(size=50)) +
  xlab("ln(homogametic lifespan/heterogametic lifespan)") +
  geom_vline(xintercept = 0) +
  theme(legend.title = element_text(size=50),
        legend.text = element_text(size=35),
        legend.justification=c(1,0), 
        legend.position=c(0.98, 0.92),  
        legend.background = element_blank(),
        legend.key = element_blank()) ##plot!

#end

