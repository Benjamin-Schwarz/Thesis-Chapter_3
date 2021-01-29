

# To run computations and analyses in this script, load functions from Schwarz_et_al_Functions.R
source("Schwarz_et_al_Functions.R")





#####################################################################
## ------------------ Plot-level analyses ------------------------ ##
#####################################################################

# read plot-level data (recorded by second observer)
dat <- read.csv("plotlevel_data.csv", stringsAsFactors = FALSE)


## Is there a difference between treatment and control plots? ----------------------------------------------------------------------
# --> compare only afternoon data

# aggregate per plot: ALL CICHORIDEAE
per.plot <- subset(dat, time.of.day%in% c("pm") & plant %in% C.pool) %>% group_by(cdate,site,run,day,time.of.day,treatment,plot) %>% summarize(focal.abun=mean(focal.abun, na.rm=T), 
                                                                                                                                               visits=sum(visits, na.rm=T),
                                                                                                                                               obs.visits=sum(obs.visits, na.rm=T),
                                                                                                                                               rates=mean(rates, na.rm=T),
                                                                                                                                               obs.rates=mean(obs.rates, na.rm=T))
per.plot[is.nan(per.plot)] <- NA

mab <- glmmTMB(round(focal.abun) ~ treatment*day + (1|site/run), per.plot, ziformula=~0, family=nbinom2)  
summary(mab)
Anova(mab)

mvi <- glmmTMB((visits) ~ treatment*day + (1|site/run), per.plot, ziformula=~0, family=nbinom2) 
summary(mvi)
Anova(mvi)

# aggregate per plot: other plants
per.plot <- subset(dat, time.of.day%in% c("pm") & !(plant %in% C.pool)) %>% group_by(cdate,site,run,day,time.of.day,treatment,plot) %>% summarize(focal.abun=mean(focal.abun, na.rm=T), 
                                                                                                                                                  visits=sum(visits, na.rm=T),
                                                                                                                                                  obs.visits=sum(obs.visits, na.rm=T),
                                                                                                                                                  rates=mean(rates, na.rm=T),
                                                                                                                                                  obs.rates=mean(obs.rates, na.rm=T))
per.plot[is.nan(per.plot)] <- NA

mvi <- glmmTMB((visits) ~ treatment*day + (1|site/run), per.plot, ziformula=~0, family=nbinom2) 
summary(mvi)
Anova(mvi)


## Is there a difference between morning and afternoon? ----------------------------------------------------------------------
# --> compare only control data

# aggregate per plot: ALL CICHORIDEAE
per.plot <- subset(dat, treatment%in% c("no") & plant %in% C.pool) %>% group_by(cdate,site,run,day,time.of.day,treatment,plot) %>% summarize(focal.abun=mean(focal.abun, na.rm=T), 
                                                                                                                                             visits=sum(visits, na.rm=T),
                                                                                                                                             obs.visits=sum(obs.visits, na.rm=T),
                                                                                                                                             rates=mean(rates, na.rm=T),
                                                                                                                                             obs.rates=mean(obs.rates, na.rm=T))
per.plot[is.nan(per.plot)] <- NA

mab <- glmmTMB(round(focal.abun) ~ time.of.day*day + (1|site/run), per.plot, ziformula=~0, family=nbinom2)
summary(mab)
Anova(mab)

mvi <- glmmTMB((visits) ~ time.of.day*day + (1|site/run), per.plot, ziformula=~0, family=nbinom2)   
summary(mvi)
Anova(mvi)

# aggregate per plot: other plants
per.plot <- subset(dat, treatment%in% c("no") & !(plant %in% C.pool)) %>% group_by(cdate,site,run,day,time.of.day,treatment,plot) %>% summarize(focal.abun=mean(focal.abun, na.rm=T), 
                                                                                                                                                visits=sum(visits, na.rm=T),
                                                                                                                                                obs.visits=sum(obs.visits, na.rm=T),
                                                                                                                                                rates=mean(rates, na.rm=T),
                                                                                                                                                obs.rates=mean(obs.rates, na.rm=T))
per.plot[is.nan(per.plot)] <- NA

mvi <- glmmTMB((visits) ~ time.of.day*day + (1|site/run), per.plot, ziformula=~0, family=nbinom2) 
summary(mvi)
Anova(mvi)


## Is there a difference between treatments on reference days? ----------------------------------------------------------------------
# --> compare only reference day data

# aggregate per plot: ALL CICHORIDEAE
per.plot <- subset(dat, day%in% c(1) & plant %in% C.pool) %>% group_by(cdate,site,run,day,time.of.day,treatment,plot) %>% summarize(focal.abun=mean(focal.abun, na.rm=T), 
                                                                                                                                    visits=sum(visits, na.rm=T),
                                                                                                                                    obs.visits=sum(obs.visits, na.rm=T),
                                                                                                                                    rates=mean(rates, na.rm=T),
                                                                                                                                    obs.rates=mean(obs.rates, na.rm=T))
per.plot[is.nan(per.plot)] <- NA

mab <- glmmTMB(round(focal.abun) ~ time.of.day*treatment + (1|site/run), per.plot, ziformula=~0, family=nbinom2)   
summary(mab)
Anova(mab)

mvi <- glmmTMB((visits) ~ time.of.day*treatment + (1|site/run), per.plot, ziformula=~0, family=nbinom2)   
summary(mvi)
Anova(mvi)



# aggregate per plot: other plants
per.plot <- subset(dat, day==1 & !(plant %in% C.pool)) %>% group_by(cdate,site,run,day,time.of.day,treatment,plot) %>% summarize(focal.abun=mean(focal.abun, na.rm=T), 
                                                                                                                                 visits=sum(visits, na.rm=T),
                                                                                                                                 obs.visits=sum(obs.visits, na.rm=T),
                                                                                                                                 rates=mean(rates, na.rm=T),
                                                                                                                                 obs.rates=mean(obs.rates, na.rm=T))
per.plot[is.nan(per.plot)] <- NA

mvi <- glmmTMB((visits) ~ time.of.day*treatment + (1|site/run), per.plot, ziformula=~0, family=nbinom2) 
summary(mvi)
Anova(mvi)



#####################################################################
## ------------------ Network analyses --------------------------- ##
#####################################################################

## load interaction data (recorded by first observer) ---------------------------
obs_clean <- read.csv("interaction_data.csv", stringsAsFactors = FALSE)

## compile dataset (net.dat) for network-level analyses ---------------------------
# create data set with only control plots in the morning, but duplicated to function as pseudo-treatment plots 
# on day 1 and 2 treatment plots are NOT used for morning data (they are replaced by control plots)
pm.dat <- subset(obs_clean, time.of.day=="afternoon")
am.control <- subset(obs_clean, time.of.day=="morning" & treatment=="no")
am.treatment <- subset(obs_clean, time.of.day=="morning" & treatment=="no")
am.treatment$treatment <- "yes" # overwrite treatment column to have data for pseudo-treatment plots
net.dat <- rbind(am.treatment, am.control, pm.dat)
net.dat <- net.dat[!is.na(net.dat$freq), ] # duplicate rows with freq > 1 (does not work if freq has NAs)
net.dat <- data.frame(net.dat[rep(seq_len(dim(net.dat)[1]), net.dat$freq), , drop = FALSE], row.names=NULL)
net.dat$freq <- as.numeric(1)
net.dat$time.of.day <- as.factor(net.dat$time.of.day)


## generate result files and store them in wd ---------------------------

## compare results against null models ---------------------------
# function all.nullmodels() produces network indices and turnover measures for all null model randomizations
# function conf.interval() uses those results and gives lower and upper confidence intervals
  
## indices
# freq null model
freq.nulls <- all.nullmodels(data=net.dat, response="indices", null.model="freq", null.N=1000)
write.table(freq.nulls, "indices.freq.nullmodels1000.csv", sep=",",row.names=F)
freq.conf <- conf.interval(all.nulls=freq.nulls)
write.table(freq.conf, "indices.freq.conf1000.csv", sep=",",row.names=F)
# freq.time null model
freq.time.nulls <- all.nullmodels(data=net.dat, response="indices", null.model="freq.time", null.N=1000)
write.table(freq.time.nulls, "indices.freq.time.nullmodels1000.csv", sep=",",row.names=F)
freq.time.conf <- conf.interval(all.nulls=freq.time.nulls)
write.table(freq.time.conf, "indices.freq.time.conf1000.csv", sep=",",row.names=F)

## turnover
turnover.nulls <- all.nullmodels(data=net.dat, response="turnover", null.model=TRUE, null.N=1000)
write.table(turnover.nulls, "turnover.nullmodels1000.csv", sep=",",row.names=F)
turnover.conf <- conf.interval(all.nulls=turnover.nulls)
write.table(turnover.conf, "turnover.conf1000.csv", sep=",",row.names=F)
# without round 4
turnover.nulls <- all.nullmodels(data=subset(net.dat, round!=4), response="turnover", null.model=TRUE, null.N=1000)
write.table(turnover.nulls, "no4turnover.nullmodels1000.csv", sep=",",row.names=F)
turnover.conf <- conf.interval(all.nulls=turnover.nulls)
write.table(turnover.conf, "no4turnover.conf1000.csv", sep=",",row.names=F)


## run permutations tests to assess treatment effects ---------------------------
# shuffle 10 plots among 2 treatments: 252 combinations

## indices
indices3 <- indices.per.day(data=net.dat, option="Option3", shared=FALSE, perms=TRUE, z=FALSE, null.n=3, null.model=FALSE)
indices3$links.to.O <- indices3$link.richness-indices3$links.to.C
write.table(indices3, "indices.orig.perms3.csv", sep=",",row.names=F)
# global effects for all experimental runs combined
effects <- indices3 %>% group_by(site,run,day,plot.comb,true.t) %>% summarize_all(difference) 
effects <- subset(effects, select=-network)
global.effect <- global.conf(data=effects, which.day=2, n.perm=5000)
write.table(global.effect, "global.orig3.day2.effect.conf.csv", sep=",",row.names=F)
global.effect <- global.conf(data=effects, which.day=1, n.perm=5000)
write.table(global.effect, "global.orig3.day1.effect.conf.csv", sep=",",row.names=F)

## turnover
turnover3 <- turnover.per.day(data=net.dat, option="Option3", shared=FALSE, perms=TRUE, null.model=FALSE)
write.table(turnover3, "turnover.orig.perms3.csv", sep=",",row.names=F)
# global effects for all experimental runs combined
effects <- turnover3 %>% group_by(site,run,day,plot.comb,true.t) %>% summarize_all(difference) 
effects <- subset(effects, select=-network)
global.effect <- global.conf(data=effects, which.day=2, n.perm=5000)
write.table(global.effect, "turnover3.day2.effect.conf.csv", sep=",",row.names=F)
global.effect <- global.conf(data=effects, which.day=1, n.perm=5000)
write.table(global.effect, "turnover3.day1.effect.conf.csv", sep=",",row.names=F)
## turnover without round 4
no4.turnover3 <- turnover.per.day(data=subset(net.dat, round!=4), option="Option3", shared=FALSE, perms=TRUE, null.model=FALSE)
write.table(no4.turnover3, "no4turnover.orig.perms3.csv", sep=",",row.names=F)
# global effects for all experimental runs combined
effects <- no4.turnover3 %>% group_by(site,run,day,plot.comb,true.t) %>% summarize_all(difference) 
effects <- subset(effects, select=-network)
global.effect <- global.conf(data=effects, which.day=2, n.perm=5000)
write.table(global.effect, "no4turnover3.day2.effect.conf.csv", sep=",",row.names=F)
global.effect <- global.conf(data=effects, which.day=1, n.perm=5000)
write.table(global.effect, "no4turnover3.day1.effect.conf.csv", sep=",",row.names=F)

## plot results ---------------------------

# as an example, produce plant generality plot (can be adjusted for other variables)
all.mods <- read.csv("indices.orig.perms3.csv")
freq.conf <- read.csv("indices.freq.conf1000.csv")
freq.time.conf <- read.csv("indices.freq.time.conf1000.csv")
plant.gen <- plot.perm.test(data=all.mods, null.model=TRUE, conf.data=freq.conf, conf.data2=freq.time.conf, 
                    variable="vulnerability.LL", ylabel="Plant generality", ylow=2.6, yhigh=13.5, ref="n.s.", treat="*")
plant.gen

# as an example, produce link turnover plot (can be adjusted for other variables)
all.mods <- read.csv("turnover.orig.perms3.csv")
turnover.conf <- read.csv("turnover.conf1000.csv")
link.turnover <- plot.perm.test(data=all.mods, null.model=TRUE, conf.data=turnover.conf, conf.data2=NULL, 
                    variable="WN", ylabel="Link turnover", ylow=0.5, yhigh=1.1, ref="n.s.", treat="*")
link.turnover




