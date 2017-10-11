

#' @title Transmission Module for Heat Bath Exposure
#'
#' @description Stochastically simulates disease transmission from the heat bath given the current
#'              state of the active HIV negative population.
#'
#' @inheritParams aging.mard
#'
#' @details
#' This is tack-on to trans.mard which is unsually the final substantive function that occurs within the time loop at
#' each time step. This function takes all active HIV negative nodes and calculates a
#' transmission probability for each individual based on age specific exposure to the heatbath
#' After transmission events, individual-level attributes for the infected
#' persons are updated and summary statistics for incidence calculated.
#'
#' The per-act transmission probability depends on the following elements:
#' an age and risk catagory specific probability of aquiring HIV from the heat bath per time step
#' THe input values \code{dat$param$heat.prob} is an Riskg by age matrix.
#'
#' @return
#' For each new infection, the disease status, infection time, and related
#' HIV attributes are updated for the infected node. Summary statistics for
#' disease incidence overall, and by race and age groups are calculated and
#' stored on \code{dat$epi}.
#'
#' @keywords module
#' @export
#'

heatbath.adol <- function(dat, at){
  
  
  ## Infected from the heat bath
  #Parameters
  age<-floor(dat$attr$age)
  riskg<-dat$attr$riskg
  race<-dat$attr$race
  active<-dat$attr$active
  status<-dat$attr$status
  debut<-dat$attr$has.debuted
  heat.contact<-dat$param$heat.prob.contact
  heat.intensity<-dat$param$heat.multiplier
  heat.trans<-dat$param$heat.trans
  everAI<-dat$attr$everAI
  pce <- dat$param$prep.class.effect
  prepStat<-dat$attr$prepStat
  prepClass<-dat$attr$prepClass
  
  #Scale the intesity of the heatbath with multiplier
  heat<-heat.contact*heat.intensity
  
  #Zero out the attributes for current number of adult contacts.
  dat$attr$AI.adult.count.t<-rep(0,length(dat$attr$AI.adult.count.t))
  
  #Select active who have debuted
  ids<-which(active==1 & dat$attr$has.debuted==1)
  
  #Determine who has Contact and who becomes infected 
  age<-age[ids]
  status<-status[ids]
  riskg<-riskg[ids]
  prepStat<-prepStat[ids]
  prepClass<-prepClass[ids]
  AI<-everAI[ids]
  infect.new<-rep(0,length(age))
  
  infect.new.mat.in<-cbind(ids,status,riskg,age,prepStat,prepClass,AI,infect.new)
  infect.new.mat.out<-NULL
  age.vec<-sort(unique(age))
  riskg.vec<-sort(unique(riskg))
  
#Calculate the number of Contacts.
  for(i in 1:length(riskg.vec)){
      for(j in 1:length(age.vec)){
        
      infect.new.mat.temp<-subset(infect.new.mat.in,infect.new.mat.in[,"riskg"]==riskg.vec[i] & infect.new.mat.in[,"age"]==age.vec[j])
      infect.new.mat.temp[,"AI"]<-rpois(length(infect.new.mat.temp[,"AI"]),heat[i,j])
      infect.new.mat.out<-rbind(infect.new.mat.out,infect.new.mat.temp)

      }
    }
  

 
    #Select those that had AI.
    AI.mat.out<-subset(infect.new.mat.out,infect.new.mat.out[,"AI"]>=1)
    AI<-AI.mat.out[,"ids"]
    AI.mat.out<-as.data.frame(AI.mat.out)

    #Calculate the transmissions
    #Include PrEP
    #Select those that had AI and are negative
    infect.new.mat.out<-subset(infect.new.mat.out,infect.new.mat.out[,"AI"]>=1 & infect.new.mat.out[,"status"]==0)
    infected<-NULL
    infect.new.mat.out<-as.data.frame(infect.new.mat.out)

    if (length(infect.new.mat.out$infect.new)>0){
      
    infect.new.mat.out.n<-subset(infect.new.mat.out, infect.new.mat.out[,"prepStat"] == 1 & infect.new.mat.out[,"prepClass"] == "n")
    infect.new.mat.out.n$infect.new<-sapply(infect.new.mat.out.n$ids, function(x) 
      max(rbinom(infect.new.mat.out.n$AI[infect.new.mat.out.n$ids==x],1,heat.trans*(1 - pce[1])),0)) 
    
    infect.new.mat.out.l<-subset(infect.new.mat.out, infect.new.mat.out[,"prepStat"] == 1 & infect.new.mat.out[,"prepClass"] == "l")
    infect.new.mat.out.l$infect.new<-sapply(infect.new.mat.out.l$ids, function(x) 
      max(rbinom(infect.new.mat.out.l$AI[infect.new.mat.out.l$ids==x],1,heat.trans*(1 - pce[2])),0))    
      
      
    infect.new.mat.out.m<-subset(infect.new.mat.out, infect.new.mat.out[,"prepStat"] == 1 & infect.new.mat.out[,"prepClass"] == "m")
    infect.new.mat.out.m$infect.new<-sapply(infect.new.mat.out.m$ids, function(x) 
      max(rbinom(infect.new.mat.out.m$AI[infect.new.mat.out.m$ids==x],1,heat.trans*(1 - pce[3])),0))   
      
    infect.new.mat.out.h<-subset(infect.new.mat.out, infect.new.mat.out[,"prepStat"] == 1 & infect.new.mat.out[,"prepClass"] == "h")
    infect.new.mat.out.h$infect.new<-sapply(infect.new.mat.out.h$ids, function(x) 
      max(rbinom(infect.new.mat.out.h$AI[infect.new.mat.out.h$ids==x],1,heat.trans*(1 - pce[4])),0))   
    
    infect.new.mat.out.none<-subset(infect.new.mat.out,infect.new.mat.out[,"prepStat"] == 0) 
    infect.new.mat.out.none$infect.new<-sapply(infect.new.mat.out.none$ids, function(x) 
      max(rbinom(infect.new.mat.out.none$AI[infect.new.mat.out.none$ids==x],1,heat.trans),0))   
    
    
    infect.new.mat.out<-rbind(infect.new.mat.out.n,infect.new.mat.out.l,infect.new.mat.out.m,infect.new.mat.out.h,infect.new.mat.out.none)
    infected<-as.vector(infect.new.mat.out$ids[infect.new.mat.out$infect.new==1])
    }
      
  # Update attributes
    AI<-as.numeric(AI)
 if (length(AI) >= 1){
     dat$attr$AI.time[AI]<-ifelse(dat$attr$everAI[AI]==0,at,dat$attr$AI.time[AI])
     dat$attr$everAI[AI]<-1
     dat$attr$AI.adult.count[AI]<-dat$attr$AI.adult.count[AI] + as.numeric(AI.mat.out$AI[AI.mat.out$ids==AI])
     dat$attr$AI.adult.count.t[AI]<- as.numeric(AI.mat.out$AI[AI.mat.out$ids==AI])
 }
    infected<-as.numeric(infected)
 if (length(infected) >= 1){
    dat$attr$status[infected] <- 1
    dat$attr$inf.time[infected] <- at
    dat$attr$vl[infected] <- 0
    dat$attr$stage[infected] <- "AR"
    dat$attr$stage.time[infected] <- 0
    dat$attr$diag.status[infected] <- 0
    dat$attr$tx.status[infected] <- 0
      
    
    dat$attr$infector[infected] <- "heat"
    dat$attr$inf.role[infected] <- "heat"
    dat$attr$inf.type[infected] <- "heat"
    dat$attr$inf.diag[infected] <- "heat" 
    dat$attr$inf.tx[infected] <- "heat"
    dat$attr$inf.stage[infected] <- "heat"
    
    dat$attr$cum.time.on.tx[infected] <- 0
    dat$attr$cum.time.off.tx[infected] <- 0
    
    
    # Summary Output
    dat$epi$incid[at] <- dat$epi$incid[at] + length(infected)
    dat$epi$incid.B[at] <- dat$epi$incid.B[at] + sum(race[infected] == "B")
    dat$epi$incid.W[at] <- dat$epi$incid.W[at] + sum(race[infected] == "W")
    dat$epi$incid.heat[at] <- length(infected)
     }
    
    if (length(infected) < 1){ dat$epi$incid.heat[at] <-0}
  
  
    
  return(dat)
}





