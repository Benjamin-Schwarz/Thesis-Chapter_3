

#################################################################
###    Functions needed to produce results and run analyses   ###
#################################################################



## packages etc ----
library(bipartite)
library(reshape2)
library(scales)
library(dplyr)
library(tidyverse)
library(tidyr)
library(tibble)
library(ggplot2) 
library(gridExtra)
library(grid)
library(scales)
library(magrittr)
library(ggpubr)
library(car)
library(glmmTMB)


## species of the Cichorieae
C.pool <- c("Cre_cap","Hyp_rad","Leo_aut","Leo_his")


## helper function to skip errors during modularity computation ---- 
setClass("ifError", representation(likelihood = "logical", modules = "logical"))
anError <- new("ifError", likelihood=NA, modules=NA) # need to create an object with slot "likelihood"
# function that returns "anError@likelihood" in case of error
eraser1 <- function (array) {return(tryCatch(computeModules(array), error=function(e) anError))}


## turn NaN into NA for whole data frame
is.nan.data.frame <- function(x) {do.call(cbind, lapply(x, is.nan))}


## calculate differences
difference <- function(x){as.numeric(x[1])-as.numeric(x[2])}


## combine two lists of matrices to one list of matrices
mat.sum <- function (...) Reduce("+", list(...))


## read data containing info about shared plants between treatment and control for each day
EXP.shared.plants <- read.csv("EXP.shared.plants.csv", na.strings=c("","NA"))



## functions to generate results -----------------------------------------------------------

## gives network level indices
indices.per.day <- function(data, option, shared, perms, z, null.n=50, null.model){
  
  results <- data.frame(site=NA, run=NA, day=NA, network=NA)
  results <- results[-1,]
  
  for(i in unique(data$site)){
    site.dat <- subset(data, site==i)
    for(j in unique(site.dat$run)){
      run.dat <- subset(site.dat, run==j)
      for(k in unique(run.dat$day)){
        dat <- subset(run.dat, day==k)
        
        if(perms==TRUE){plot.comb <- combn(unique(data$plot), 5, FUN = NULL, simplify = F)} else {plot.comb <- list(unique(subset(run.dat, time.of.day=="afternoon" & treatment=="yes")$plot))}
        
        for(l in 1:length(plot.comb)){
          
          dat$random.t <- ifelse(dat$plot %in% plot.comb[[l]], "yes", "no")
          # re-assign treatment to morning data manually as there are only 5 control plots
          # thus, for each network, morning data represent observations in the 5 control plots and only afternoon data are different combinations of 5 out of 10 plots
          dat$random.t[dat$time.of.day=="morning"] <- as.character(dat$treatment[dat$time.of.day=="morning"])
          if(option %in% c("Option1")) {dat$webID <- paste(dat$time.of.day, dat$random.t)}
          if(option %in% c("Option2", "Option3")) {dat$webID <- dat$random.t}
          
          # # optionally, exclude non-shared plant species between control and treatment
          # shared.plants <- subset(EXP.shared.plants, site==i & run==j & day==k)[,-c(1,2,3)]
          # shared.plants <- names(shared.plants[,colSums(is.na(shared.plants)) == 0])
          # if(shared==TRUE) {dat <- subset(dat, lower %in% shared.plants)}
          
          # select nullmodel
          if(null.model==FALSE) {webarray <- frame2webs(dat, varnames = c("lower", "higher", "webID", "freq"), type.out="array", empty=TRUE)}
          if(null.model=="freq") {
            array <- frame2webs(dat, varnames = c("lower", "higher", "webID", "freq"), type.out="array", empty=TRUE)
            cnull <- nullmodel(array[,,1], N=1, method='r2d')
            tnull <- nullmodel(array[,,2], N=1, method='r2d')
            webarray <- simplify2array(list(matrix(unlist(cnull), nrow=dim(array)[1], ncol=dim(array)[2]),
                                            matrix(unlist(tnull), nrow=dim(array)[1], ncol=dim(array)[2])))
            dimnames(webarray) <- dimnames(array)
          }
          if(null.model=="freq.time") {
            am.dat <- dat
            am.dat[am.dat$time.of.day=="afternoon",]$freq <- 0
            amarray <- frame2webs(am.dat, varnames = c("lower", "higher", "webID", "freq"), type.out="array", empty=F)
            amcnull <- nullmodel(amarray[,,1], N=1, method='r2d')
            amtnull <- nullmodel(amarray[,,2], N=1, method='r2d')
            
            pm.dat <- dat
            pm.dat[pm.dat$time.of.day=="morning",]$freq <- 0
            pmarray <- frame2webs(pm.dat, varnames = c("lower", "higher", "webID", "freq"), type.out="array", empty=F)
            pmcnull <- nullmodel(pmarray[,,1], N=1, method='r2d')
            pmtnull <- nullmodel(pmarray[,,2], N=1, method='r2d')
            
            cnull <- mapply(mat.sum, amcnull, pmcnull, SIMPLIFY = FALSE) 
            tnull <- mapply(mat.sum, amtnull, pmtnull, SIMPLIFY = FALSE)
            webarray <- simplify2array(list(matrix(unlist(cnull), nrow=dim(amarray)[1], ncol=dim(amarray)[2]),
                                            matrix(unlist(tnull), nrow=dim(amarray)[1], ncol=dim(amarray)[2])))
            dimnames(webarray) <- dimnames(amarray)
          }
          
          # calculate modularity and number of modules
          mod.list <- apply(webarray, 3, eraser1)
          mod.values <- unlist(lapply(mod.list, function(x) `@`(x ,likelihood)[[1]]))
          module.list <- lapply(mod.list, function(x) `@`(x ,modules))
          mod.number <- unlist(lapply(module.list, function(x) {if(!is.na(x)) {length(x[,1])-1} else {NA}}))
          
          # create empty vectors for the number of links that connect modules, the proportion of interactions made up by those links, and the CV of mean module time
          cross.number <- NULL
          cross.prop <- NULL
          mod.cv <- NULL
          # # fill empty vectors
          for (w in 1:length(dimnames(webarray)[[3]])){
            # module membership
            modules <- if(!is.na(module.list[[w]])) {module.list[[w]][-1, -c(1:2)]} else {NA}
            mod.df <- as.data.frame(t(modules))
            names(mod.df) <- letters[1:length(names(mod.df))]
            species.list <- data.frame(species=c(dimnames(bipartite::empty(webarray[,,w]))[[1]], dimnames(bipartite::empty(webarray[,,w]))[[2]]), module=NA)
            species.list$module <- if(is.na(mod.df[1,1])) {NA} else {apply(mod.df, 1, function(x) paste(names(mod.df)[x != 0]))}
            module.df <- left_join(subset(dat, webID==dimnames(webarray)[[3]][w]), species.list, by=c("lower"="species"))
            module.df <- left_join(module.df, species.list, by=c("higher"="species"))
            module.df <- module.df[!(is.na(module.df$module.y)),]
            module.df$cross.link <- ifelse(module.df$module.x==module.df$module.y, "no", "yes")
            cross <- subset(module.df, cross.link=="yes") %>% group_by(link) %>% summarize(freq=sum(freq))
            cross.links.number <- length(cross$link)
            cross.number <- c(cross.number, cross.links.number)
            cross.links.prop <- sum(cross$freq)/sum(module.df$freq)
            cross.prop <- c(cross.prop, cross.links.prop)
            time <- subset(module.df, cross.link=="no") %>% group_by(module.x) %>% summarize(time=mean(norm.time, na.rm=T))
            mod.cv <- c(mod.cv, sd(time$time)/mean(time$time))
          }
          
          
          # calculate other indices
          metric.values <- apply(webarray, 3, networklevel, index=c("connectance", "H2", "generality"))   #network- and group-level metrics
          mymetrics <- as.data.frame(t(metric.values))
          
          # calculate generality of Cichorideae visitors (remove pollinators with 0 visits to Cichorideae)
          if(i=="Freiburg"){foc.C <- "Cre_cap"} else{foc.C <- "Leo_his"}
          generality.C.visitors <- apply(webarray, 3, function(x) {if(sum(x[, apply(x[(dimnames(x)[[1]] %in% C.pool),,drop=F], 2, sum) != 0])>0) {grouplevel(x[, apply(x[(dimnames(x)[[1]] %in% C.pool),,drop=F], 2, sum) != 0], index="generality", level="higher")} else{NA}})
          C.generality <- apply(webarray, 3, function(x) {if(sum(x[(dimnames(x)[[1]] %in% C.pool),,drop=F])>0) {grouplevel(x[(dimnames(x)[[1]] %in% C.pool),,drop=F], index="generality", level="lower")} else{NA}})
          nonC.generality <-  apply(webarray, 3, function(x) {if(sum(x[!(dimnames(x)[[1]] %in% C.pool),,drop=F])>0) {grouplevel(x[!(dimnames(x)[[1]] %in% C.pool),,drop=F], index="generality", level="lower")} else{NA}})
          mymetrics <- cbind(mymetrics, generality.C.visitors, C.generality, nonC.generality)
          
          # optionally: calculate z-scores
          z.mod <- rep(NA, length(dimnames(webarray)[[3]]))
          z.H2 <- rep(NA, length(dimnames(webarray)[[3]]))
          
          time.z.mod <- rep(NA, length(dimnames(webarray)[[3]]))
          time.z.H2 <- rep(NA, length(dimnames(webarray)[[3]]))
          
          mod.lower <- rep(NA, length(dimnames(webarray)[[3]]))
          mod.upper <- rep(NA, length(dimnames(webarray)[[3]]))
          H2.lower <- rep(NA, length(dimnames(webarray)[[3]]))
          H2.upper <- rep(NA, length(dimnames(webarray)[[3]]))
          
          time.mod.lower <- rep(NA, length(dimnames(webarray)[[3]]))
          time.mod.upper <- rep(NA, length(dimnames(webarray)[[3]]))
          time.H2.lower <- rep(NA, length(dimnames(webarray)[[3]]))
          time.H2.upper <- rep(NA, length(dimnames(webarray)[[3]]))
          
          if(z==TRUE & perms==FALSE) {
            # make subwebs for morning and afternoon to randomize interactions within subwebs
            am.dat <- dat
            am.dat[am.dat$time.of.day=="afternoon",]$freq <- 0
            amarray <- frame2webs(am.dat, varnames = c("lower", "higher", "webID", "freq"), type.out="array", empty=F)
            pm.dat <- dat
            pm.dat[pm.dat$time.of.day=="morning",]$freq <- 0
            pmarray <- frame2webs(pm.dat, varnames = c("lower", "higher", "webID", "freq"), type.out="array", empty=F)
            
            for(n in 1:length(dimnames(webarray)[[3]])){
              nulls <- nullmodel(webarray[,,n], N=null.n, method='r2d')
              
              mods <- sapply(nulls, eraser1)
              mods.mean <- mean(unlist(lapply(mods, function(x) `@`(x ,likelihood)[[1]])), na.rm=T)
              mods.sd <- sd(unlist(lapply(mods, function(x) `@`(x ,likelihood)[[1]])), na.rm=T)
              z.mod[n] <- (mod.values[n]-mods.mean)/mods.sd
              mod.lower[n] <- quantile(unlist(lapply(mods, function(x) `@`(x ,likelihood)[[1]])), probs = 0.025, na.rm=TRUE)
              mod.upper[n] <- quantile(unlist(lapply(mods, function(x) `@`(x ,likelihood)[[1]])), probs = 0.975, na.rm=TRUE)
              
              hs <- sapply(nulls, networklevel, index="H2")
              z.H2[n] <- (mymetrics$H2[n]-mean(hs, na.rm=T))/sd(hs, na.rm=T)
              H2.lower[n] <- quantile(hs, probs = 0.025, na.rm=TRUE)
              H2.upper[n] <- quantile(hs, probs = 0.975, na.rm=TRUE)
              
              # include time of the day in the nullmodel
              am.nulls <- nullmodel(amarray[,,n], N=null.n, method='r2d')
              pm.nulls <- nullmodel(pmarray[,,n], N=null.n, method='r2d')
             
              time.nulls <- mapply(mat.sum, am.nulls, pm.nulls, SIMPLIFY = FALSE) 
              
              time.mods <- sapply(time.nulls, eraser1)
              time.mods.mean <- mean(unlist(lapply(time.mods, function(x) `@`(x ,likelihood)[[1]])), na.rm=T)
              time.mods.sd <- sd(unlist(lapply(time.mods, function(x) `@`(x ,likelihood)[[1]])), na.rm=T)
              time.z.mod[n] <- (mod.values[n]-time.mods.mean)/time.mods.sd
              time.mod.lower[n] <- quantile(unlist(lapply(time.mods, function(x) `@`(x ,likelihood)[[1]])), probs = 0.025, na.rm=TRUE)
              time.mod.upper[n] <- quantile(unlist(lapply(time.mods, function(x) `@`(x ,likelihood)[[1]])), probs = 0.975, na.rm=TRUE)
              
              time.hs <- sapply(time.nulls, networklevel, index="H2")
              time.z.H2[n] <- (mymetrics$H2[n]-mean(time.hs, na.rm=T))/sd(time.hs, na.rm=T)
              time.H2.lower[n] <- quantile(time.hs, probs = 0.025, na.rm=TRUE)
              time.H2.upper[n] <- quantile(time.hs, probs = 0.975, na.rm=TRUE)
              
            }
          } 
          
          # species, links
          plant.spp <- apply(webarray, 3, networklevel, index="number of species", level="lower")
          pollinator.spp <- apply(webarray, 3, networklevel, index="number of species", level="higher")
          freq <- apply(webarray, 3, sum, na.rm=T)
          links <- apply(webarray, 3, function(c)sum(c!=0))
          links.to.C <- apply(webarray[intersect(dimnames(webarray)[[1]], c("Cre_cap","Leo_his","Leo_aut","Hyp_rad")),,,drop=F], 3, function(c)sum(c!=0))
          prop.C.spec <- apply(webarray[intersect(dimnames(webarray)[[1]], c("Cre_cap","Leo_his","Leo_aut","Hyp_rad")),intersect(dimnames(webarray)[[2]], c("Panurgus_calcaratus","Eupeodes_corollae","Andrena_flavipes","Lasioglossum_malachurum","Lasioglossum_villosulum","Lasioglossum_leucozonium")),,drop=F], 3, sum, na.rm=T)/
            apply(webarray[intersect(dimnames(webarray)[[1]], c("Cre_cap","Leo_his","Leo_aut","Hyp_rad")),,,drop=F], 3, sum, na.rm=T)
          
          # count how many plots of random combination match true treatment plots to identify true results
          true.comb <- length(intersect(plot.comb[[l]], unique(subset(run.dat, time.of.day=="afternoon" & treatment=="yes")$plot)))
          
          # make results
          new.rows <- nrow(results)+(1:length(dimnames(webarray)[[3]]))
          results[new.rows, "plot.comb"] <- l
          results[new.rows, "true.t"] <- true.comb
          results[new.rows, "site"] <- i
          results[new.rows, "run"] <- j
          results[new.rows, "day"] <- k
          results[new.rows, "network"] <- dimnames(webarray)[[3]]
          results[new.rows, "modularity"] <- mod.values
          results[new.rows, "mod.lower"] <- mod.lower
          results[new.rows, "mod.upper"] <- mod.upper
          results[new.rows, "time.mod.lower"] <- time.mod.lower
          results[new.rows, "time.mod.upper"] <- time.mod.upper
          results[new.rows, "z.mod"] <- z.mod
          results[new.rows, "time.z.mod"] <- time.z.mod
          for(m in 1:length(mymetrics)){results[new.rows, names(mymetrics[m])] <- mymetrics[,m]}
          results[new.rows, "H2.lower"] <- H2.lower
          results[new.rows, "H2.upper"] <- H2.upper
          results[new.rows, "time.H2.lower"] <- time.H2.lower
          results[new.rows, "time.H2.upper"] <- time.H2.upper
          results[new.rows, "z.H2"] <- z.H2
          results[new.rows, "time.z.H2"] <- time.z.H2
          results[new.rows, "module.number"] <- mod.number
          results[new.rows, "cross.number"] <- cross.number
          results[new.rows, "cross.prop"] <- cross.prop
          results[new.rows, "mod.time.cv"] <- mod.cv
          results[new.rows, "plant.richness"] <- plant.spp
          results[new.rows, "pollinator.richness"] <- pollinator.spp
          results[new.rows, "link.richness"] <- links
          results[new.rows, "links.to.C"] <- links.to.C
          results[new.rows, "prop.C.spec"] <- prop.C.spec
          results[new.rows, "freq"] <- freq
        }
      }
    }
  }
  results[is.nan(results)] <- NA
  return(results)
}



## function to randomize data for null models
make_null <- function(data, null.model, model){
  if(null.model==FALSE){dat <- data}
  if(null.model==TRUE){
    if(model=="time"){   #assign time completely randomly 
     dat$time.of.day <- sample(c("morning","afternoon"), size=length(dat$time.of.day), replace=TRUE) 
    }
    if(model=="time.freq"){   #consider freq per time
      tdat <- subset(data, random.t=="yes")
      tdat$time.of.day <- sample(subset(data, random.t=="yes")$time.of.day)
      cdat <- subset(data, random.t=="no")
      cdat$time.of.day <- sample(subset(data, random.t=="no")$time.of.day)
      dat <- rbind(tdat, cdat)
    }
    if(model=="time.freq.plant"){   #consider freq per time and plant group
      Ctdat <- subset(data, random.t=="yes" & lower%in%C.pool)
      Ctdat$time.of.day <- sample(Ctdat$time.of.day)
      Ccdat <- subset(data, random.t=="no" & lower%in%C.pool)
      Ccdat$time.of.day <- sample(Ccdat$time.of.day)
      Otdat <- subset(data, random.t=="yes" & !(lower%in%C.pool))
      Otdat$time.of.day <- sample(Otdat$time.of.day)
      Ocdat <- subset(data, random.t=="no" & !(lower%in%C.pool))
      Ocdat$time.of.day <- sample(Ocdat$time.of.day)
      dat <- rbind(Ctdat, Ccdat, Otdat, Ocdat)
    }
    if(model=="plant"){   #assign plant group randomly but consider freq per plant group
      tdat <- subset(data, random.t=="yes")
      tdat$lower <- sample(tdat$lower)
      cdat <- subset(data, random.t=="no")
      cdat$lower <- sample(cdat$lower)
      dat <- rbind(tdat, cdat)
    }
    if(model=="plant.time"){   #assign plant group randomly but consider freq per plant group and plant composition per time of day
      amtdat <- subset(data, random.t=="yes" & time.of.day=="morning")
      amtdat$lower <- sample(amtdat$lower)
      amcdat <- subset(data, random.t=="no" & time.of.day=="morning")
      amcdat$lower <- sample(amcdat$lower)
      pmtdat <- subset(data, random.t=="yes" & time.of.day=="afternoon")
      pmtdat$lower <- sample(pmtdat$lower)
      pmcdat <- subset(data, random.t=="no" & time.of.day=="afternoon")
      pmcdat$lower <- sample(pmcdat$lower)
      dat <- rbind(amtdat, amcdat, pmtdat, pmcdat)
    }
  }
  return(dat)
}



## gives dissimilarities between morning and afternoon webs
turnover.per.day <- function(data, option, shared, perms, null.model){
  
  results <- data.frame(site=NA, run=NA, day=NA, network=NA)
  results <- results[-1,]
  
  for(i in unique(data$site)){
    site.dat <- subset(data, site==i)
    for(j in unique(site.dat$run)){
      run.dat <- subset(site.dat, run==j)
      for(k in unique(run.dat$day)){
        day.dat <- subset(run.dat, day==k)
        
        if(perms==TRUE){plot.comb <- combn(unique(data$plot), 5, FUN = NULL, simplify = F)} else {plot.comb <- list(unique(subset(run.dat, time.of.day=="afternoon" & treatment=="yes")$plot))}
        
        for(l in 1:length(plot.comb)){
          
          day.dat$random.t <- ifelse(day.dat$plot %in% plot.comb[[l]], "yes", "no")
          # re-assign treatment to morning data manually as there are only 5 control plots
          # thus, for each network, morning data represent observations in the 5 control plots and only afternoon data are different combinations of 5 out of 10 plots
          day.dat$random.t[day.dat$time.of.day=="morning"] <- as.character(day.dat$treatment[day.dat$time.of.day=="morning"])
          if(option %in% c("Option1")) {day.dat$webID <- paste(day.dat$time.of.day, day.dat$random.t)}
          if(option %in% c("Option2", "Option3")) {day.dat$webID <- day.dat$random.t}
          
          # # optionally, exclude non-shared plant species between control and treatment
          # shared.plants <- subset(EXP.shared.plants, site==i & run==j & day==k)[,-c(1,2,3)]
          # shared.plants <- names(shared.plants[,colSums(is.na(shared.plants)) == 0])
          # if(shared==TRUE) {dat <- subset(dat, lower %in% shared.plants)}
          
          # count how many plots of random combination match true treatment plots to identify true results
          true.comb <- length(intersect(plot.comb[[l]], unique(subset(run.dat, time.of.day=="afternoon" & treatment=="yes")$plot)))
          
          ## quantify within-day dynamics (only possible for Options 2 and 3)
          if(option %in% c("Option2", "Option3")) {

            ## quantitative dissimilarity
            # select appropriate nullmodel
            dat <- make_null(data=day.dat, null.model=null.model, model="time.freq")
            
            # select treatment plots
            tarray <- frame2webs(subset(dat, random.t=="yes"), varnames = c("lower", "higher", "time.of.day", "freq"), type.out="array")
            # calculate species turnover, rewiring, link turnover
            tbetavalues <- betalinkr.prop(tarray, index="jaccard", partitioning="poisot", binary=F, proportions=T)
            tbetavalues <- as.data.frame(t(unlist(tbetavalues)))
            # plant and pollinator turnover
            tplants <- rbind(apply(tarray[,,"afternoon", drop=F], 1, sum), apply(tarray[,,"morning", drop=F], 1, sum))
            tbetavalues$plant.turnover <- c(vegdist(decostand(tplants, method="total"), method="jaccard", binary=FALSE))
            tpolls <- rbind(apply(tarray[,,"afternoon", drop=F], 2, sum), apply(tarray[,,"morning", drop=F], 2, sum))
            tbetavalues$pollinator.turnover <- c(vegdist(decostand(tpolls, method="total"), method="jaccard", binary=FALSE))
            # select control plots
            carray <- frame2webs(subset(dat, random.t=="no"), varnames = c("lower", "higher", "time.of.day", "freq"), type.out="array")
            # calculate species turnover, rewiring, link turnover
            cbetavalues <- betalinkr.prop(carray, index="jaccard", partitioning="poisot", binary=F, proportions=T)
            cbetavalues <- as.data.frame(t(unlist(cbetavalues)))
            # plant and pollinator turnover
            cplants <- rbind(apply(carray[,,"afternoon", drop=F], 1, sum), apply(carray[,,"morning", drop=F], 1, sum))
            cbetavalues$plant.turnover <- c(vegdist(decostand(cplants, method="total"), method="jaccard", binary=FALSE))
            cpolls <- rbind(apply(carray[,,"afternoon", drop=F], 2, sum), apply(carray[,,"morning", drop=F], 2, sum))
            cbetavalues$pollinator.turnover <- c(vegdist(decostand(cpolls, method="total"), method="jaccard", binary=FALSE))
            
            betavalues <- rbind(cbetavalues, tbetavalues) # first control as it has to be same order as for mymetrics
            
            # binary dissimilarity
            # calculate species turnover, rewiring, link turnover
            tbetavalues.b <- betalinkr.prop(tarray, index="jaccard", partitioning="poisot", binary=T, proportions=T)
            tbetavalues.b <- as.data.frame(t(unlist(tbetavalues.b)))
            # plant and pollinator turnover
            tplants.b <- rbind(apply(tarray[,,"afternoon", drop=F], 1, sum), apply(tarray[,,"morning", drop=F], 1, sum))
            tbetavalues.b$plant.turnover.b <- c(vegdist(decostand(tplants, method="total"), method="jaccard", binary=TRUE))
            tpolls.b <- rbind(apply(tarray[,,"afternoon", drop=F], 2, sum), apply(tarray[,,"morning", drop=F], 2, sum))
            tbetavalues.b$pollinator.turnover.b <- c(vegdist(decostand(tpolls.b, method="total"), method="jaccard", binary=TRUE))
            # calculate species turnover, rewiring, link turnover
            cbetavalues.b <- betalinkr.prop(carray, index="jaccard", partitioning="poisot", binary=F, proportions=T)
            cbetavalues.b <- as.data.frame(t(unlist(cbetavalues.b)))
            # plant and pollinator turnover
            cplants.b <- rbind(apply(carray[,,"afternoon", drop=F], 1, sum), apply(carray[,,"morning", drop=F], 1, sum))
            cbetavalues.b$plant.turnover.b <- c(vegdist(decostand(cplants.b, method="total"), method="jaccard", binary=TRUE))
            cpolls.b <- rbind(apply(carray[,,"afternoon", drop=F], 2, sum), apply(carray[,,"morning", drop=F], 2, sum))
            cbetavalues.b$pollinator.turnover.b <- c(vegdist(decostand(cpolls.b, method="total"), method="jaccard", binary=TRUE))
            
            betavalues.b <- rbind(cbetavalues.b, tbetavalues.b) # first control as it has to be same order as for mymetrics
            names(betavalues.b) <- c("S.b", "OS.b", "WN.b", "ST.b", "plant.turnover.b", "pollinator.turnover.b")
            
            ## quantitative dissimilarity of C and O subwebs
            subdat <- make_null(data=day.dat, null.model, model="time.freq.plant")   #select appropriate nullmodel
            # Cichorideae subweb
            # select treatment plots
            Ctarray <- frame2webs(subset(subdat, random.t=="yes" & lower %in% c("Cre_cap","Leo_his", "Hyp_rad", "Leo_aut")), varnames = c("lower", "higher", "time.of.day", "freq"), type.out="array", emptylist=T)
            # calculate species turnover, rewiring, link turnover
            Ctbetavalues <- betalinkr.prop(Ctarray, index="jaccard", partitioning="poisot", binary=F, proportions=T)
            Ctbetavalues <- as.data.frame(t(unlist(Ctbetavalues)))
            # plant and pollinator turnover
            Ctplants <- rbind(apply(Ctarray[,,"afternoon",drop = F], 1, sum), apply(Ctarray[,,"morning",drop = F], 1, sum))
            Ctbetavalues$Cplant.turnover <- c(vegdist(decostand(Ctplants, method="total"), method="jaccard", binary=FALSE))
            Ctpolls <- rbind(apply(Ctarray[,,"afternoon",drop = F], 2, sum), apply(Ctarray[,,"morning",drop = F], 2, sum))
            Ctbetavalues$Cpollinator.turnover <- c(vegdist(decostand(Ctpolls, method="total"), method="jaccard", binary=FALSE))
            # select control plots
            Ccarray <- frame2webs(subset(subdat, random.t=="no" & lower %in% c("Cre_cap","Leo_his", "Hyp_rad", "Leo_aut")), varnames = c("lower", "higher", "time.of.day", "freq"), type.out="array")
            # calculate species turnover, rewiring, link turnover
            Ccbetavalues <- betalinkr.prop(Ccarray, index="jaccard", partitioning="poisot", binary=F, proportions=T)
            Ccbetavalues <- as.data.frame(t(unlist(Ccbetavalues)))
            # plant and pollinator turnover
            Ccplants <- rbind(apply(Ccarray[,,"afternoon",drop = F], 1, sum), apply(Ccarray[,,"morning",drop = F], 1, sum))
            Ccbetavalues$Cplant.turnover <- c(vegdist(decostand(Ccplants, method="total"), method="jaccard", binary=FALSE))
            Ccpolls <- rbind(apply(Ccarray[,,"afternoon",drop = F], 2, sum), apply(Ccarray[,,"morning",drop = F], 2, sum))
            Ccbetavalues$Cpollinator.turnover <- c(vegdist(decostand(Ccpolls, method="total"), method="jaccard", binary=FALSE))
            Cbetavalues <- rbind(Ccbetavalues, Ctbetavalues) # first control as it has to be same order as for mymetrics
            names(Cbetavalues) <- c("C.S", "C.OS", "C.WN", "C.ST", "C.plant.turnover", "C.pollinator.turnover")
            # others subweb
            # select treatment plots
            Otarray <- frame2webs(subset(subdat, random.t=="yes" & !(lower %in% c("Cre_cap","Leo_his", "Hyp_rad", "Leo_aut"))), varnames = c("lower", "higher", "time.of.day", "freq"), type.out="array")
            # calculate species turnover, rewiring, link turnover
            Otbetavalues <- betalinkr.prop(webarray=Otarray, index="jaccard", partitioning="poisot", binary=F, proportions=T)
            Otbetavalues <- as.data.frame(t(unlist(Otbetavalues)))
            # plant and pollinator turnover
            Otplants <- rbind(apply(Otarray[,,"afternoon",drop = F], 1, sum), apply(Otarray[,,"morning",drop = F], 1, sum))
            Otbetavalues$Oplant.turnover <- c(vegdist(decostand(Otplants, method="total"), method="jaccard", binary=FALSE))
            Otpolls <- rbind(apply(Otarray[,,"afternoon",drop = F], 2, sum), apply(Otarray[,,"morning",drop = F], 2, sum))
            Otbetavalues$Opollinator.turnover <- c(vegdist(decostand(Otpolls, method="total"), method="jaccard", binary=FALSE))
            # select control plots
            Ocarray <- frame2webs(subset(subdat, random.t=="no" & !(lower %in% c("Cre_cap","Leo_his", "Hyp_rad", "Leo_aut"))), varnames = c("lower", "higher", "time.of.day", "freq"), type.out="array")
            # calculate species turnover, rewiring, link turnover
            Ocbetavalues <- betalinkr.prop(webarray=Ocarray, index="jaccard", partitioning="poisot", binary=F, proportions=T)
            Ocbetavalues <- as.data.frame(t(unlist(Ocbetavalues)))
            # plant and pollinator turnover
            Ocplants <- rbind(apply(Ocarray[,,"afternoon",drop = F], 1, sum), apply(Ocarray[,,"morning",drop = F], 1, sum))
            Ocbetavalues$Oplant.turnover <- c(vegdist(decostand(Ocplants, method="total"), method="jaccard", binary=FALSE))
            Ocpolls <- rbind(apply(Ocarray[,,"afternoon",drop = F], 2, sum), apply(Ocarray[,,"morning",drop = F], 2, sum))
            Ocbetavalues$Opollinator.turnover <- c(vegdist(decostand(Ocpolls, method="total"), method="jaccard", binary=FALSE))
            Obetavalues <- rbind(Ocbetavalues, Otbetavalues) # first control as it has to be same order as for mymetrics
            names(Obetavalues) <- c("O.S", "O.OS", "O.WN", "O.ST", "O.plant.turnover", "O.pollinator.turnover")
            
            ## rewiring between Cichorieae and other plant species
            plant.dat <- day.dat
            plant.dat$lower <- ifelse(plant.dat$lower %in% C.pool, "Cichorieae", "other")
            null.plant <- make_null(data=plant.dat, null.model, model="plant.time")   #select appropriate nullmodel
            # select treatment plots
            tparray <- frame2webs(subset(null.plant, random.t=="yes"), varnames = c("lower", "higher", "time.of.day", "freq"), type.out="array")
            # calculate species turnover, rewiring, link turnover
            tpbetavalues <- betalinkr.prop(tparray, index="jaccard", partitioning="poisot", binary=F, proportions=T)
            tpbetavalues <- as.data.frame(t(unlist(tpbetavalues)))
            # select control plots
            cparray <- frame2webs(subset(null.plant, random.t=="no"), varnames = c("lower", "higher", "time.of.day", "freq"), type.out="array")
            # calculate species turnover, rewiring, link turnover
            cpbetavalues <- betalinkr.prop(cparray, index="jaccard", partitioning="poisot", binary=F, proportions=T)
            cpbetavalues <- as.data.frame(t(unlist(cpbetavalues)))
            
            ## assign pm visits to different origin (from C, from O, new) --> use dat from time.freq null model
            #treatment plots
            am.prop <- tparray[,,"morning"] / rbind(colSums(tparray[,,"morning"]), colSums(tparray[,,"morning"]))
            am.prop[is.nan(am.prop)] <- 0
            new <- ifelse(colSums(am.prop)==0, 1, 0)
            am.prop <- rbind(am.prop, new)   #matrix with probability of visiting a plant group in the morning for each pollinator
            fromC <- tparray[,,"afternoon"] * rbind(am.prop[1,], am.prop[1,])
            fromO <- tparray[,,"afternoon"] * rbind(am.prop[2,], am.prop[2,])
            fromNew <- tparray[,,"afternoon"] * rbind(am.prop[3,], am.prop[3,])
            tpoll.move <- data.frame(C.C=rowSums(fromC)[1], C.O=rowSums(fromC)[2], O.C=rowSums(fromO)[1], O.O=rowSums(fromO)[2], new.C=rowSums(fromNew)[1], new.O=rowSums(fromNew)[2])
            #tpoll.move <- tpoll.move / sum(tpoll.move)   #relative proportion of all visits
            #control plots
            cam.prop <- cparray[,,"morning"] / rbind(colSums(cparray[,,"morning"]), colSums(cparray[,,"morning"]))
            cam.prop[is.nan(cam.prop)] <- 0
            cnew <- ifelse(colSums(cam.prop)==0, 1, 0)
            cam.prop <- rbind(cam.prop, cnew)   #matrix with probability of visiting a plant group in the morning for each pollinator
            cfromC <- cparray[,,"afternoon"] * rbind(cam.prop[1,], cam.prop[1,])
            cfromO <- cparray[,,"afternoon"] * rbind(cam.prop[2,], cam.prop[2,])
            cfromNew <- cparray[,,"afternoon"] * rbind(cam.prop[3,], cam.prop[3,])
            cpoll.move <- data.frame(C.C=rowSums(cfromC)[1], C.O=rowSums(cfromC)[2], O.C=rowSums(cfromO)[1], O.O=rowSums(cfromO)[2], new.C=rowSums(cfromNew)[1], new.O=rowSums(cfromNew)[2])
            #cpoll.move <- cpoll.move / sum(cpoll.move)   #relative proportion of all visits
            poll.move <- rbind(cpoll.move, tpoll.move)   #first control as it has to be same order as for mymetrics
            poll.move$Cswitcher <- poll.move$C.O / (poll.move$C.C+poll.move$C.O)
            poll.move$Oswitcher <- poll.move$O.C / (poll.move$O.O+poll.move$O.C)
            
            ## pollinator dissimilarity (C-->C, C-->O, O-->C, O-->O, amCO, pmCO)
            plant.dat <- day.dat
            plant.dat$lower <- ifelse(plant.dat$lower %in% C.pool, "Cichorieae", "other")
           
            # C-->C
            CC <- rbind(subset(plant.dat, time.of.day=="morning" & lower=="Cichorieae"), subset(plant.dat, time.of.day=="afternoon" & lower=="Cichorieae"))
            CC <- make_null(data=CC, null.model, model="time.freq")   #select appropriate nullmodel
            CCpolls <- frame2webs(CC, varnames = c("time.of.day", "higher", "random.t", "freq"), type.out="array")
            # treatment plots (if permutation leads to zero treatment plots, NA is produced)
            if(!("yes" %in% dimnames(CCpolls)[[3]])) {tpbetavalues$amC.pmC <- NA} else{tpbetavalues$amC.pmC <- c(vegdist(decostand(CCpolls[,,"yes"], method="total"), method="jaccard", binary=FALSE))}
            # control plots (if permutation leads to zero control plots, NA is produced)
            if(!("no" %in% dimnames(CCpolls)[[3]])) {cpbetavalues$amC.pmC <- NA} else{cpbetavalues$amC.pmC <- c(vegdist(decostand(CCpolls[,,"no"], method="total"), method="jaccard", binary=FALSE))}

            # C-->O
            CO <- rbind(subset(plant.dat, time.of.day=="morning" & lower=="Cichorieae"), subset(plant.dat, time.of.day=="afternoon" & lower=="other"))
            CO <- make_null(data=CO, null.model, model="time.freq")   #select appropriate nullmodel
            COpolls <- frame2webs(CO, varnames = c("time.of.day", "higher", "random.t", "freq"), type.out="array")
            # treatment plots (if permutation leads to zero treatment plots, NA is produced)
            if(!("yes" %in% dimnames(COpolls)[[3]])) {tpbetavalues$amC.pmO <- NA} else{tpbetavalues$amC.pmO <- c(vegdist(decostand(COpolls[,,"yes"], method="total"), method="jaccard", binary=FALSE))}
            # control plots (if permutation leads to zero control plots, NA is produced)
            if(!("no" %in% dimnames(COpolls)[[3]])) {cpbetavalues$amC.pmO <- NA} else{cpbetavalues$amC.pmO <- c(vegdist(decostand(COpolls[,,"no"], method="total"), method="jaccard", binary=FALSE))}
            
            # O-->C
            OC <- rbind(subset(plant.dat, time.of.day=="morning" & lower=="other"), subset(plant.dat, time.of.day=="afternoon" & lower=="Cichorieae"))
            OC <- make_null(data=OC, null.model, model="time.freq")   #select appropriate nullmodel
            OCpolls <- frame2webs(OC, varnames = c("time.of.day", "higher", "random.t", "freq"), type.out="array")
            # treatment plots (if permutation leads to zero treatment plots, NA is produced)
            if(!("yes" %in% dimnames(OCpolls)[[3]])) {tpbetavalues$amO.pmC <- NA} else{tpbetavalues$amO.pmC <- c(vegdist(decostand(OCpolls[,,"yes"], method="total"), method="jaccard", binary=FALSE))}
            # control plots (if permutation leads to zero control plots, NA is produced)
            if(!("no" %in% dimnames(OCpolls)[[3]])) {cpbetavalues$amO.pmC <- NA} else{cpbetavalues$amO.pmC <- c(vegdist(decostand(OCpolls[,,"no"], method="total"), method="jaccard", binary=FALSE))}
            
            # O-->O
            OO <- rbind(subset(plant.dat, time.of.day=="morning" & lower=="other"), subset(plant.dat, time.of.day=="afternoon" & lower=="other"))
            OO <- make_null(data=OO, null.model, model="time.freq")   #select appropriate nullmodel
            OOpolls <- frame2webs(OO, varnames = c("time.of.day", "higher", "random.t", "freq"), type.out="array")
            # treatment plots (if permutation leads to zero treatment plots, NA is produced)
            if(!("yes" %in% dimnames(OOpolls)[[3]])) {tpbetavalues$amO.pmO <- NA} else{tpbetavalues$amO.pmO <- c(vegdist(decostand(OOpolls[,,"yes"], method="total"), method="jaccard", binary=FALSE))}
            # control plots (if permutation leads to zero control plots, NA is produced)
            if(!("no" %in% dimnames(OOpolls)[[3]])) {cpbetavalues$amO.pmO <- NA} else{cpbetavalues$amO.pmO <- c(vegdist(decostand(OOpolls[,,"no"], method="total"), method="jaccard", binary=FALSE))}
            
            
            # amCO
            amCO <- rbind(subset(plant.dat, time.of.day=="morning" & lower=="Cichorieae"), subset(plant.dat, time.of.day=="morning" & lower=="other"))
            amCO <- make_null(data=amCO, null.model, model="plant")   #select appropriate nullmodel
            amCOpolls <- (frame2webs(amCO, varnames = c("lower", "higher", "random.t", "freq"), type.out="array"))
            # treatment plots (if permutation leads to zero treatment plots, NA is produced)
            if(!("yes" %in% dimnames(amCOpolls)[[3]])) {tpbetavalues$amC.amO <- NA} else{tpbetavalues$amC.amO <- c(vegdist(decostand(amCOpolls[,,"yes"], method="total"), method="jaccard", binary=FALSE))}
            # control plots (if permutation leads to zero control plots, NA is produced)
            if(!("no" %in% dimnames(amCOpolls)[[3]])) {cpbetavalues$amC.amO <- NA} else{cpbetavalues$amC.amO <- c(vegdist(decostand(amCOpolls[,,"no"], method="total"), method="jaccard", binary=FALSE))}
            
            
            # pmCO
            pmCO <- rbind(subset(plant.dat, time.of.day=="afternoon" & lower=="Cichorieae"), subset(plant.dat, time.of.day=="afternoon" & lower=="other"))
            pmCO <- make_null(data=pmCO, null.model, model="plant")   #select appropriate nullmodel
            pmCOpolls <- (frame2webs(pmCO, varnames = c("lower", "higher", "random.t", "freq"), type.out="array"))
            # treatment plots (if permutation leads to zero treatment plots, NA is produced)
            if(!("yes" %in% dimnames(pmCOpolls)[[3]])) {tpbetavalues$pmC.pmO <- NA} else{tpbetavalues$pmC.pmO <- c(vegdist(decostand(pmCOpolls[,,"yes"], method="total"), method="jaccard", binary=FALSE))}
            # control plots (if permutation leads to zero control plots, NA is produced)
            if(!("no" %in% dimnames(pmCOpolls)[[3]])) {cpbetavalues$pmC.pmO <- NA} else{cpbetavalues$pmC.pmO <- c(vegdist(decostand(pmCOpolls[,,"no"], method="total"), method="jaccard", binary=FALSE))}
            
            betavalues.p <- rbind(cpbetavalues, tpbetavalues) # first control as it has to be same order as for mymetrics
            names(betavalues.p) <- c("pS","pOS","pWN","pST","amC.pmC","amC.pmO","amO.pmC","amO.pmO","amC.amO","pmC.pmO")
            
            # determine pollinator switching
            betavalues.p$CtoO <- betavalues.p$amC.pmO/(betavalues.p$amC.amO)
            betavalues.p$OtoC <- betavalues.p$amO.pmC/(betavalues.p$amC.amO)
            
            mymetrics <- as.data.frame(cbind(betavalues, betavalues.b, Cbetavalues, Obetavalues, betavalues.p, poll.move))
          }
          
          # make results
          new.rows <- nrow(results)+(1:2)
          results[new.rows, "plot.comb"] <- l
          results[new.rows, "true.t"] <- true.comb
          results[new.rows, "site"] <- i
          results[new.rows, "run"] <- j
          results[new.rows, "day"] <- k
          results[new.rows, "network"] <- c("no","yes")
          for(m in 1:length(mymetrics)){results[new.rows, names(mymetrics[m])] <- mymetrics[,m]}
        }
      }
    }
  }
  results[is.nan(results)] <- NA
  return(results)
}



## generate null model results
all.nullmodels <- function(data, response, null.N, null.model, model){
  
  if(response=="indices"){
    all.nulls <- indices.per.day(data, option="Option3", shared=FALSE, perms=FALSE, z=FALSE, null.n=1, null.model)
    all.nulls$n <- 1
    
    for(n in 1:(null.N-1)){
      nulls <- indices.per.day(data=net.dat, option="Option3", shared=FALSE, perms=FALSE, z=FALSE, null.n=1, null.model)
      nulls$n <- n+1
      all.nulls <- rbind(all.nulls, nulls)
    }
  }
  
  if(response=="turnover"){
    all.nulls <- turnover.per.day(data=net.dat, option="Option3", shared=FALSE, perms=FALSE, null.model=TRUE)
    all.nulls$n <- 1
    
    for(n in 1:(null.N-1)){
      nulls <- turnover.per.day(data, option="Option3", shared=FALSE, perms=FALSE, null.model=TRUE)
      nulls$n <- n+1
      all.nulls <- rbind(all.nulls, nulls)
    }
  }
  
  if(response=="niche"){
    all.nulls <- temp.niche(data=net.dat, option="Option3", shared=FALSE, perms=FALSE, null.n=1, null.model=TRUE)
    all.nulls$n <- 1
    
    for(n in 1:(null.N-1)){
      nulls <- temp.niche(data=net.dat, option="Option3", shared=FALSE, perms=FALSE, null.n=1, null.model=TRUE)
      nulls$n <- n+1
      all.nulls <- rbind(all.nulls, nulls)
    }
  }
  return(all.nulls)
}



## provide confidence intervals for null models
conf.interval <- function(all.nulls){
  
  # determine means of six experimental runs per day-treatment combination
  all.means <- all.nulls %>% group_by(n, day, network) %>% summarize_all(mean, na.rm=TRUE)
  
  # determine confidence intervals for means
  lower <- all.means[,!(colnames(all.means) %in% c("site", "run", "plot.comb", "true.t", "perm"))] %>% group_by(day, network) %>% 
    summarize_all(quantile, probs = 0.025, na.rm=TRUE)
  lower$n <- "lower"
  
  median <- all.means[,!(colnames(all.means) %in% c("site", "run", "plot.comb", "true.t", "perm"))] %>% group_by(day, network) %>% 
    summarize_all(quantile, probs = 0.5, na.rm=TRUE)
  median$n <- "median"
  
  upper <- all.means[,!(colnames(all.means) %in% c("site", "run", "plot.comb", "true.t", "perm"))] %>% group_by(day, network) %>% 
    summarize_all(quantile, probs = 0.975, na.rm=TRUE)
  upper$n <- "upper"
  
  conf <- rbind(lower,median,upper)
  conf <- conf[with(conf, order(day, network, n)), ]
  return(conf)
}




## gives results of permutation tests in relation to true values
global.conf <- function(data, which.day, n.perm){
  set.seed(55)
  
  all_perms <- as.data.frame(data[0,!(colnames(data) %in% c("site", "run", "day", "plot.comb", "true.t", "perm"))])
  
  for(i in 1:n.perm){
    comb <- sample(1:252, 6, replace=T)
    
    perm <- rbind(
      subset(data, site=="Ebnet" & run==1 & day==which.day & plot.comb==comb[1]),
      subset(data, site=="Ebnet" & run==2 & day==which.day & plot.comb==comb[2]),
      subset(data, site=="Freiburg" & run==1 & day==which.day & plot.comb==comb[3]),
      subset(data, site=="Freiburg" & run==2 & day==which.day & plot.comb==comb[4]),
      subset(data, site=="Kirchzarten" & run==1 & day==which.day & plot.comb==comb[5]),
      subset(data, site=="Kirchzarten" & run==2 & day==which.day & plot.comb==comb[6])
    )
    
    mean <- as.data.frame(t(colMeans(perm[sapply(perm, is.numeric)], na.rm=T)))
    
    # adding ID for each row
    mean$perm <- i
    
    all_perms <- rbind(all_perms, mean)
  }
  
  # determine confidence intervals for index values across all sites and runs
  
  lower <- all_perms[,!(colnames(all_perms) %in% c("site", "run", "day", "plot.comb", "true.t", "perm"))] %>% 
    summarize_all(quantile, probs = 0.025, na.rm=TRUE)
  
  median <- all_perms[,!(colnames(all_perms) %in% c("site", "run", "day", "plot.comb", "true.t", "perm"))] %>% 
    summarize_all(quantile, probs = 0.5, na.rm=TRUE)
  
  upper <- all_perms[,!(colnames(all_perms) %in% c("site", "run", "day", "plot.comb", "true.t", "perm"))] %>% 
    summarize_all(quantile, probs = 0.975, na.rm=TRUE)
  
  global.conf <- rbind(lower,median,upper)
  
  # merge with real index value for treatment
  true <- rbind(
    subset(data, site=="Ebnet" & run==1 & day==which.day & true.t==5),
    subset(data, site=="Ebnet" & run==2 & day==which.day & true.t==5),
    subset(data, site=="Freiburg" & run==1 & day==which.day & true.t==5),
    subset(data, site=="Freiburg" & run==2 & day==which.day & true.t==5),
    subset(data, site=="Kirchzarten" & run==1 & day==which.day & true.t==5),
    subset(data, site=="Kirchzarten" & run==2 & day==which.day & true.t==5)
  )
  
  true[is.nan(true)] <- NA
  true.mean <- as.data.frame(t(colMeans(true[sapply(true, is.numeric)], na.rm=T)))
  global.effect <- rbind(global.conf, true.mean[,colnames(true.mean) %in% colnames(global.conf)])
  global.effect <- data.frame(index=names(global.effect),lower=t(global.effect)[,1],median=t(global.effect)[,2],upper=t(global.effect)[,3],true=t(global.effect)[,4], row.names=NULL)
  
  return(global.effect)
}



## nice plot to compare network-level indices across treatments
plot.perm.test <- function(data, null.model, conf.data, conf.data2, variable, ylabel, ylow, yhigh, ref, treat){
  
  ggdat <- subset(subset(data, true.t==5) %>% group_by(day, network) %>% 
                    summarise_each(funs(mean(., na.rm=T), sd = sd(., na.rm=T), se = sd(., na.rm=T)/sqrt(sum(!is.na(.)))), variable),
                  network %in% c("yes", "no"))
  
  if(null.model==TRUE){
    null.dat <- dcast(conf.data, day + network ~ n, value.var=variable)
    ggdat <- left_join(ggdat, null.dat, by = c("day", "network"))
    if(!is.null(conf.data2)){
      null.dat2 <- dcast(conf.data2, day + network ~ n, value.var=variable)
      null.dat2 <- null.dat2 %>% rename(lower.time = lower, upper.time = upper)
      ggdat <- left_join(ggdat, null.dat2, by = c("day", "network"))
    }
  }
  
  
  ggdat$day[ggdat$day==1] <- "Reference day"
  ggdat$day[ggdat$day==2] <- "Treatment day"
  ggdat$network <- gsub("no", "Control", ggdat$network)
  ggdat$network <- gsub("yes", "Treatment", ggdat$network)
  
  
  
  a <- ggplot(ggdat, aes(x=day, y=mean, fill=network)) +
    geom_errorbar(position=position_dodge(.7), width=.25, aes(ymin=mean-se, ymax=mean+se)) +
    geom_point(position=position_dodge(.7), shape=21, size=5) +
    geom_errorbar(position=position_dodge(.7), color="#CCCCCC", width=0, aes(ymin=mean+se+0.03*(yhigh-ylow), ymax=ylow+0.9*(yhigh-ylow))) +
    geom_errorbar(position=position_dodge(.0), color="#CCCCCC", width=0.7, aes(ymin=ylow+0.9*(yhigh-ylow), ymax=ylow+0.9*(yhigh-ylow))) +
    annotate("text",x="Reference day",y=ylow+0.95*(yhigh-ylow),label=ref) +
    annotate("text",x="Treatment day",y=ylow+0.95*(yhigh-ylow),label=treat) +
    scale_fill_manual(values=c("#CCCCCC","#FFFFFF")) +
    scale_y_continuous(name=ylabel) +
    coord_cartesian(ylim=c(ylow, yhigh)) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
  
  if(null.model==TRUE){
    
    b <- ggplot(ggdat) +
      geom_rect(xmin = c(0.5,1,1.5,2), xmax = c(1,1.5,2,2.5), ymin = ggdat$lower, ymax = ggdat$upper, fill="black", alpha = 0.4) +
      geom_errorbar(aes(x=day, ymin=mean-se, ymax=mean+se, group=network), position=position_dodge(1), width=.25) +
      geom_point(aes(x=day, y=mean, fill=network), position=position_dodge(1), shape=21, size=5) +
      geom_errorbar(aes(x=day, ymin=mean+se+0.03*(yhigh-ylow), ymax=ylow+0.9*(yhigh-ylow), group=network), position=position_dodge(1), color="#CCCCCC", width=0) +
      geom_errorbar(aes(x=day, ymin=ylow+0.9*(yhigh-ylow), ymax=ylow+0.9*(yhigh-ylow), group=network), position=position_dodge(.0), color="#CCCCCC", width=1) +
      #geom_segment(x=2-0.1, xend=2+0.1, y=ggdat$mean[3]+0.1*(yhigh-ylow), yend=ggdat$mean[4]+0.1*(yhigh-ylow), 
      #             size = 1, color="white", arrow = arrow(length = unit(0.3, "cm"), type="closed")) +
      annotate("text",x="Reference day",y=ylow+0.95*(yhigh-ylow),label=ref) +
      annotate("text",x="Treatment day",y=ylow+0.95*(yhigh-ylow),label=treat) +
      scale_fill_manual(values=c("#CCCCCC","#FFFFFF")) +
      geom_vline(xintercept = 1.5) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_continuous(name=ylabel) +
      coord_cartesian(ylim=c(ylow, yhigh)) +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            panel.grid = element_blank(),
            legend.title = element_blank(),
            legend.position = "none")
    
    c <- ggplot(ggdat) +
      geom_rect(xmin = c(0.5,1,1.5,2), xmax = c(0.95,1.45,1.95,2.45), ymin = ggdat$lower, ymax = ggdat$upper, fill="black", alpha = 0.4) +
      geom_rect(xmin = c(0.55,1.05,1.55,2.05), xmax = c(1,1.5,2,2.5), ymin = ggdat$lower.time, ymax = ggdat$upper.time, fill="darkblue", alpha = 0.4) +
      geom_errorbar(aes(x=day, ymin=mean-se, ymax=mean+se, group=network), position=position_dodge(1), width=.25) +
      geom_point(aes(x=day, y=mean, fill=network), position=position_dodge(1), shape=21, size=5) +
      geom_errorbar(aes(x=day, ymin=mean+se+0.03*(yhigh-ylow), ymax=ylow+0.9*(yhigh-ylow), group=network), position=position_dodge(1), color="#CCCCCC", width=0) +
      geom_errorbar(aes(x=day, ymin=ylow+0.9*(yhigh-ylow), ymax=ylow+0.9*(yhigh-ylow), group=network), position=position_dodge(.0), color="#CCCCCC", width=1) +
      #geom_segment(x=2-0.1, xend=2+0.1, y=ggdat$mean[3]+0.1*(yhigh-ylow), yend=ggdat$mean[4]+0.1*(yhigh-ylow), 
      #             size = 1, color="white", arrow = arrow(length = unit(0.3, "cm"), type="closed")) +
      annotate("text",x="Reference day",y=ylow+0.95*(yhigh-ylow),label=ref) +
      annotate("text",x="Treatment day",y=ylow+0.95*(yhigh-ylow),label=treat) +
      scale_fill_manual(values=c("#CCCCCC","#FFFFFF")) +
      geom_vline(xintercept = 1.5) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_continuous(name=ylabel) +
      coord_cartesian(ylim=c(ylow, yhigh)) +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            panel.grid = element_blank(),
            legend.title = element_blank(),
            legend.position = "none")
  }
  if(null.model==TRUE & is.null(conf.data2)) {return(b)} else if(null.model==TRUE & !is.null(conf.data2)) {return(c)} else {return(a)}
}



## nice plot to compare network-level indices across treatments --> ONLY REFERENCE DAYS FOR MAIN TEXT
plot.perm.test <- function(data, null.model, conf.data, conf.data2, variable, ylabel, ylow, yhigh, ref, treat){
  
  ggdat <- subset(subset(data, true.t==5 & day==2) %>% group_by(day, network) %>% 
                    summarise_each(funs(mean(., na.rm=T), sd = sd(., na.rm=T), se = sd(., na.rm=T)/sqrt(sum(!is.na(.)))), variable),
                  network %in% c("yes", "no"))
  
  if(null.model==TRUE){
    null.dat <- dcast(conf.data, day + network ~ n, value.var=variable)
    ggdat <- left_join(ggdat, null.dat, by = c("day", "network"))
    if(!is.null(conf.data2)){
      null.dat2 <- dcast(conf.data2, day + network ~ n, value.var=variable)
      null.dat2 <- null.dat2 %>% rename(lower.time = lower, upper.time = upper)
      ggdat <- left_join(ggdat, null.dat2, by = c("day", "network"))
    }
  }
  
  
  ggdat$day[ggdat$day==1] <- "Reference day"
  ggdat$day[ggdat$day==2] <- "Treatment day"
  ggdat$network <- gsub("no", "Control", ggdat$network)
  ggdat$network <- gsub("yes", "Treatment", ggdat$network)
  
  
  
  a <- ggplot(ggdat, aes(x=network, y=mean, fill=network)) +
    geom_errorbar(position=position_dodge(.7), width=.25, aes(ymin=mean-se, ymax=mean+se)) +
    geom_point(position=position_dodge(.7), shape=21, size=5) +
    geom_errorbar(position=position_dodge(.7), color="#CCCCCC", width=0, aes(ymin=mean+se+0.03*(yhigh-ylow), ymax=ylow+0.9*(yhigh-ylow))) +
    geom_segment(color="#CCCCCC", aes(x=1, xend=2, y=ylow+0.9*(yhigh-ylow), yend=ylow+0.9*(yhigh-ylow))) +
    #annotate("text",x="Reference day",y=ylow+0.95*(yhigh-ylow),label=ref) +
    annotate("text",x=1.5,y=ylow+0.95*(yhigh-ylow),label=treat) +
    scale_fill_manual(values=c("#CCCCCC","#FFFFFF")) +
    scale_y_continuous(name=ylabel) +
    coord_cartesian(ylim=c(ylow, yhigh)) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
  
  if(null.model==TRUE){

    b <- ggplot(ggdat) +
      geom_rect(xmin = c(0.5,1.5), xmax = c(1.5,2.5), ymin = ggdat$lower, ymax = ggdat$upper, fill="black", alpha = 0.4) +
      geom_errorbar(aes(x=network, ymin=mean-se, ymax=mean+se, group=network), position=position_dodge(1), width=.25) +
      geom_point(aes(x=network, y=mean, fill=network), position=position_dodge(1), shape=21, size=5) +
      geom_errorbar(aes(x=network, ymin=mean+se+0.03*(yhigh-ylow), ymax=ylow+0.9*(yhigh-ylow), group=network), position=position_dodge(1), color="#CCCCCC", width=0) +
      geom_segment(color="#CCCCCC", aes(x=1, xend=2, y=ylow+0.9*(yhigh-ylow), yend=ylow+0.9*(yhigh-ylow))) +
      #annotate("text",x="Reference day",y=ylow+0.95*(yhigh-ylow),label=ref) +
      annotate("text",x=1.5,y=ylow+0.95*(yhigh-ylow),label=treat) +
      scale_fill_manual(values=c("#CCCCCC","#FFFFFF")) +
      #geom_vline(xintercept = 1.5) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_continuous(name=ylabel) +
      coord_cartesian(ylim=c(ylow, yhigh)) +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            panel.grid = element_blank(),
            legend.title = element_blank(),
            legend.position = "none")
    
    c <- ggplot(ggdat) +
      geom_rect(xmin = c(0.5,1.5), xmax = c(1.4,2.4), ymin = ggdat$lower, ymax = ggdat$upper, fill="black", alpha = 0.4) +
      geom_rect(xmin = c(0.6,1.6), xmax = c(1.5,2.5), ymin = ggdat$lower.time, ymax = ggdat$upper.time, fill="darkblue", alpha = 0.4) +
      geom_errorbar(aes(x=network, ymin=mean-se, ymax=mean+se, group=network), position=position_dodge(1), width=.25) +
      geom_point(aes(x=network, y=mean, fill=network), position=position_dodge(1), shape=21, size=5) +
      geom_errorbar(aes(x=network, ymin=mean+se+0.03*(yhigh-ylow), ymax=ylow+0.9*(yhigh-ylow), group=network), position=position_dodge(1), color="#CCCCCC", width=0) +
      geom_segment(color="#CCCCCC", aes(x=1, xend=2, y=ylow+0.9*(yhigh-ylow), yend=ylow+0.9*(yhigh-ylow))) +
      #annotate("text",x="Reference day",y=ylow+0.95*(yhigh-ylow),label=ref) +
      annotate("text",x=1.5,y=ylow+0.95*(yhigh-ylow),label=treat) +
      scale_fill_manual(values=c("#CCCCCC","#FFFFFF")) +
      #geom_vline(xintercept = 1.5) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_continuous(name=ylabel) +
      coord_cartesian(ylim=c(ylow, yhigh)) +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            panel.grid = element_blank(),
            legend.title = element_blank(),
            legend.position = "none")
  }
  if(null.model==TRUE & is.null(conf.data2)) {return(b)} else if(null.model==TRUE & !is.null(conf.data2)) {return(c)} else {return(a)}
}



## modification of bipartite::betalinkr() to calculate species and link dissimilarity (turnover) between two networks (here: morning and afternoon)
betalinkr.prop <- function(webarray, index = "bray", binary=TRUE, partitioning="commondenom", proportions=!binary, function.dist="vegdist", distofempty="zero", partition.st=FALSE, partition.rr=FALSE){
  
  if (class(webarray)=="list") {webarray <- webs2array(webarray)}
  if (dim(webarray)[[3]]!=2) warning("function is designed for a single pair of two webs; unclear output")
  
  webarray <- webarray[apply(webarray,1,sum)>0, apply(webarray,2,sum)>0, , drop=FALSE] # removing species not observed in either web: improves speed :-)
  
  # standardizing to proportions if wanted (now on webarray with marginal totals; decostand only used for method!=commondenom)
  if (proportions){
    if (binary){warning("standardizing to proportions for binary index; do you really want to do this?!?")}
    array.of.sums <- webarray
    
    #############################################################################
    array.of.sums[,,1] <- ifelse(sum(webarray[,,1])==0, 1, sum(webarray[,,1]))   # avoid division by 0
    array.of.sums[,,2] <- ifelse(sum(webarray[,,2])==0, 1, sum(webarray[,,2]))
    #############################################################################
    
    webarray <- webarray / array.of.sums
  }
  
  if (binary==TRUE){
    webarray[webarray>0] <- 1
  }  
  
  # for "shared species subweb", set non-shared species links to zero
  array.sharedsp <- webarray
  array.sharedsp[rowSums(apply(webarray, MARGIN=c(1,3), sum)>0)!=2, , ] <- 0
  array.sharedsp[, rowSums(apply(webarray, MARGIN=c(2,3), sum)>0)!=2, ] <- 0
  
  # all links
  linkmx <- array2linkmx(webarray)
  # only links of shared species
  linkmx.sharedsp <- array2linkmx(array.sharedsp)  
  
  # removing links never observed in complete webarray; same subset for both linkmx!
  #links.observed <- colSums(linkmx)>0
  #linkmx <- linkmx[, links.observed]
  #linkmx.sharedsp <- linkmx.sharedsp[, links.observed]
  #-> not using this as it actually made the function slower! (in one test case)
  
  # species community matrix (combining upper and lower of bipartite web)
  # now ensuring a single entry per species also in named foodwebs (adding ingoing and outgoing links, instead of counting them separately)
  specmx.lower <- apply(webarray, c(3,1), sum)
  specmx.higher <- apply(webarray, c(3,2), sum)
  
  specmx.higher.unique <- specmx.higher[, !(colnames(specmx.higher) %in% colnames(specmx.lower)), drop=F]  ### add drop=F so that one-column matrices are not turned into one-row matrices
  #specmx.higher.unique <- specmx.higher[, !(colnames(specmx.higher) %in% colnames(specmx.lower))]         ### original code by Jochen
  
  if (is.null(colnames(specmx.higher))) {specmx.higher.unique <- specmx.higher}  # if no species names are given, assuming bipartite webs
  specmx.all <- cbind(specmx.lower, specmx.higher.unique)  # e.g. sites X (plants, pollinators); 
  duplicolnames <- setdiff(colnames(specmx.higher), colnames(specmx.higher.unique))
  specmx.all[, duplicolnames] <- specmx.all[, duplicolnames] + specmx.higher[, duplicolnames]
  
  if (partitioning=="poisot"){
    if (partition.st | partition.rr){warning("further partitioning only available with method partitioning='commondenom'")}
    
    # alternative subsets of linkmx (partitioning ST and OS) --
    # my approach is to set the species/links to zero instead of excluding them (which will be done anyways when dissimilarity is calculated)
    # this makes it easier to match species / links (even without names)
    # shared links of shared species (only LINKS occurring in both sites)
    linkmx.sharedli <- linkmx   
    linkmx.sharedli[, colSums(linkmx.sharedli>0)==1] <- 0
    # varying links of shared species
    linkmx.rewiring <- linkmx.sharedsp - linkmx.sharedli
    linkmx.RewSha <- linkmx.sharedsp  # all links excluding those from unique species
    # links of non-shared species
    linkmx.uniquesp <- linkmx - linkmx.sharedsp
    linkmx.UniSha <- linkmx.uniquesp + linkmx.sharedli # all links excluding rewiring links
    
    # standardizing to proportions if wanted
    if (proportions){
      linkmx.RewSha <- decostand(linkmx.RewSha, method="total")
      linkmx.UniSha <- decostand(linkmx.UniSha, method="total")
    }
    
    # calculating dissimilarity / the betalink components --
    if (function.dist=="vegdist"){
      b_s <- vegdist(specmx.all, method=index, binary=binary) # "S"
      b_wn <- vegdist(linkmx, method=index, binary=binary) # "WN"
      b_zero <- b_wn  # preparation for "distofempty": first get structure of distance matrix
      b_zero[] <- 0   # preparation for "distofempty": define a zero distance matrix
      if (distofempty=="zero" & any(rowSums(linkmx.RewSha)==0)){  # set to conceptually correct value; avoids warning
        b_os.raw <- b_zero # no shared species means zero contribution of OS
      } else {
        b_os.raw <- vegdist(linkmx.RewSha, method=index, binary=binary) # "OS"
      }
    } 
    if (function.dist=="betadiver"){
      if (binary==FALSE) {
        warning("betadiver only uses binary data; for quantitative indices use vegdist")
      } else {
        b_s <- betadiver(specmx.all, method=index) # "S"
        b_wn <- betadiver(linkmx, method=index) # "WN"
        if (distofempty=="zero" & any(rowSums(linkmx.RewSha)==0)){  # set to the conceptually correct value; avoids warning
          b_os.raw <- b_zero # no shared species means zero contribution of OS
        } else {
          b_os.raw <- betadiver(linkmx.RewSha, method=index) # "OS"
        }
      }
    }
    # output (and final steps of calculation)
    b_os <- b_os.raw
    b_st <- b_wn - b_os.raw
    return(c(S=b_s, OS=b_os, WN=b_wn, ST=b_st))
  }
  
  if (partitioning=="commondenom"){
    # here, index must be explicitly specified as one of Sorensen or Jaccard (lower case)
    # quantitative equivalents available with binary=F; I actually use the quantitative formulas, which simplify to binary index if given binary data 
    
    # A, B, C follows Legendre 2014; index "tot"=total differences, "rew"=rewiring differences, "uni"=unique species differences
    # I can calculate all components with linkmx and linkmx.sharedsp, and their difference
    pmins <- pmin(linkmx[1,], linkmx[2,])
    A <- sum(pmins)
    B.tot <- sum(linkmx[1,] - pmins)
    B.rew <- sum(linkmx.sharedsp[1,] - pmin(linkmx.sharedsp[1,], linkmx.sharedsp[2,])) # note: frequency-changes of shared interactions also included in rewiring here!
    B.uni <- B.tot - B.rew  # here it works with subtraction
    C.tot <- sum(linkmx[2,] - pmins)
    C.rew <- sum(linkmx.sharedsp[2,] - pmin(linkmx.sharedsp[1,], linkmx.sharedsp[2,]))
    C.uni <- C.tot - C.rew
    if (index == "bray") {index <- "sorensen"}  # for convenience, bray is also allowed (but avoiding this in first place to be consistent with Legendre-terminology)
    if (index == "sorensen"){denominator <- 2*A + B.tot + C.tot}
    if (index == "jaccard"){denominator <- A + B.tot + C.tot}
    b_wn <- (B.tot + C.tot) / denominator
    b_os <- (B.rew + C.rew) / denominator
    b_st <- (B.uni + C.uni) / denominator 
    b_s <- vegdist(specmx.all, method=switch(index, "jaccard"="jaccard", "sorensen"="bray"), binary=binary) # "S"
    if (partition.st){
      # partition ST further into ST.l, ST.h and ST.lh
      # this seems easier to do when calculating with the arrays first, converting to linkmx only later
      # lh = all - sharedhighORlow
      # l = all - sharedlow - sharedboth - lh
      # h = all - sharedhigh - sharedboth - lh
      array.sharedlow <- webarray
      array.sharedlow[rowSums(apply(webarray, MARGIN=c(1,3), sum)>0)!=2, , ] <- 0
      array.sharedhigh <- webarray
      array.sharedhigh[, rowSums(apply(webarray, MARGIN=c(2,3), sum)>0)!=2, ] <- 0
      array.sharedhighORlow <- webarray
      array.sharedhighORlow[rowSums(apply(webarray, MARGIN=c(1,3), sum)>0)!=2, rowSums(apply(webarray, MARGIN=c(2,3), sum)>0)!=2, ] <- 0
      linkmx.lh <- array2linkmx(webarray - array.sharedhighORlow)  # interactions with both partners unique to one of the webs
      linkmx.l <- array2linkmx(webarray - array.sharedlow) - linkmx.lh # interactions with lower sp unique to one of the webs (but higher species shared)
      linkmx.h <- array2linkmx(webarray - array.sharedhigh) - linkmx.lh # interactions with higher sp unique to one of the webs (but lower species shared)
      # Note: the pmin-part can be omitted here, as it will always be zero!
      B.l <- sum(linkmx.l[1,])
      B.h <- sum(linkmx.h[1,])
      B.lh <- sum(linkmx.lh[1,])
      # B.l + B.h + B.lh == B.uni  # check, must be true
      C.l <- sum(linkmx.l[2,])
      C.h <- sum(linkmx.h[2,])
      C.lh <- sum(linkmx.lh[2,])
      # C.l + C.h + C.lh == C.uni  # check, must be true
      b_st.l <- (B.l + C.l) / denominator 
      b_st.h <- (B.h + C.h) / denominator 
      b_st.lh <- (B.lh + C.lh) / denominator 
    }
    if (partition.rr){
      # partition WN and OS further into link Replacement component and link Richness difference component
      # for binary=FALSE, link richness difference component is the dissimilarity component due to differences in network totals
      if (proportions==TRUE){warning("partitionining into replacement and richness (abundance) difference components may be meaningless with proportions")}
      b_wn.repl <- 2*min(B.tot, C.tot) / denominator
      b_os.repl <- 2*min(B.rew, C.rew) / denominator
      b_wn.rich <- abs(B.tot - C.tot) / denominator
      b_os.rich <- abs(B.rew - C.rew) / denominator
    }
    # output (concetenated for secondary partitionings)
    output <- c(S=b_s, OS=b_os, WN=b_wn, ST=b_st)
    if (partition.st==TRUE){
      output <- c(output, ST.l=b_st.l, ST.h=b_st.h, ST.lh=b_st.lh)
    }
    if (partition.rr==TRUE){
      output <- c(output, WN.repl=b_wn.repl, OS.repl=b_os.repl, WN.rich=b_wn.rich, OS.rich=b_os.rich)
    }
    return(output)
  }
}


