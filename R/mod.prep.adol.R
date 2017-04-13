
#' @title PrEP Module
#'
#' @description Module function for implementation and uptake of pre-exposure
#'              prophylaxis (PrEP) to prevent HIV infection.
#'
#' @inheritParams aging.mard
#'
#' @export
#'
prep.adol <- function(dat, at) {

  if (at < dat$param$prep.start) {
    return(dat)
  }
  

  ## Variables
  active <- dat$attr$active
  status <- dat$attr$status
  diag.status <- dat$attr$diag.status
  debuted<-dat$attr$has.debuted
  debuted.time<-dat$attr$has.debuted.time
  active.time<-dat$attr$active.time
  AI.time<-dat$attr$AI.time
  lnt <- dat$attr$last.neg.test
  prep.delay<-dat$param$prep.delay
  everAI<-dat$attr$everAI
  of.age<-dat$attr$of.age
  uaicount<-dat$attr$uaicount
  
  prepElig <- dat$attr$prepElig
  prepStat <- dat$attr$prepStat
  prepEver <- dat$attr$prepEver
  prepClass <- dat$attr$prepClass
  prepStart.time<-dat$attr$prepStart.time

  prep.elig.model <- dat$param$prep.elig.model
  prep.coverage <- dat$param$prep.coverage
  prep.cov.method <- dat$param$prep.cov.method
  prep.cov.rate <- dat$param$prep.cov.rate
  prep.class.prob <- dat$param$prep.class.prob
  prep.risk.reassess <- dat$param$prep.risk.reassess
  prepSpell<- dat$param$prepSpell
  prepDrop<- dat$param$prepDrop
  prep.uaicount.thresh<-dat$param$prep.uaicount.thresh

  ## Eligibility ---------------------------------------------------------------

  # Base eligibility
  idsEligStart <- which(active == 1 & status == 0 & prepStat == 0)

  idsEligStop <- NULL
  if (prep.risk.reassess == TRUE) {
    idsEligStop <- which(active == 1 & prepStat == 1)
  }

  if (prep.elig.model == "none")  prepElig[idsEligStart] <-0
  
  # Core eligiblity scenarios
  if (prep.elig.model != "none") {
    if (substr(prep.elig.model, 1, 4) == "adol") {
      if (prep.elig.model == "adol.entry") {
        c1 <- active
        c2 <- active
        c3 <- active
      } else if (prep.elig.model == "adol.debuted") {
        c1 <- active
        c2 <- debuted
        c3 <- debuted
      } else if (prep.elig.model == "adol.AI") {
        c1 <- active
        c2 <- debuted
        c3 <- everAI
      } else if (prep.elig.model == "adol.entry.older") {
        c1 <- active
        c2 <- of.age
        c3 <- of.age
      } else if (prep.elig.model == "adol.debuted.older") {
        c1 <- active
        c2 <- debuted
        c3 <- of.age
      } else if (prep.elig.model == "adol.AI.older") {
        c1 <- active
        c2 <- everAI
        c3 <- of.age
      } else if (prep.elig.model == "adol.entry.time") {
        c1 <- active
        c2 <- ifelse (at - active.time > prep.delay,1,0) 
        c3 <- ifelse (at - active.time > prep.delay,1,0)
      } else if (prep.elig.model == "adol.debuted.time") {
        c1 <- debuted
        c2 <- ifelse (at - debuted.time > prep.delay,1,0) 
        c3 <- ifelse (at - debuted.time > prep.delay,1,0) 
      } else if (prep.elig.model == "adol.AI.time") {
        c1 <- everAI
        c2 <- ifelse (at - AI.time > prep.delay,1,0) 
        c3 <- ifelse (at - AI.time > prep.delay,1,0)
      } else if (prep.elig.model == "adol.entry.older.time") {
        c1 <- active
        c2 <- of.age
        c3 <- ifelse (at - active.time > prep.delay,1,0)
      } else if (prep.elig.model == "adol.debuted.older.time") {
        c1 <- debuted
        c2 <- of.age
        c3 <- ifelse (at - debuted.time > prep.delay,1,0)
      } else if (prep.elig.model == "adol.AI.older.time") {
        c1 <- everAI
        c2 <- of.age
        c3 <- ifelse (at - AI.time > prep.delay,1,0)
      } else if (prep.elig.model == "adol.riskhist") {
        c1 <- ifelse (uaicount> prep.uaicount.thresh,1,0)
        c2 <- ifelse (uaicount> prep.uaicount.thresh,1,0)
        c3 <- ifelse (uaicount> prep.uaicount.thresh,1,0)
      } else if (prep.elig.model == "adol.riskhist.older") {
        c1 <- ifelse (uaicount> prep.uaicount.thresh,1,0)
        c2 <- of.age
        c3 <- of.age
      } else if (prep.elig.model == "adol.riskhist.time") {
        c1 <- ifelse (uaicount> prep.uaicount.thresh,1,0)
        c2 <- ifelse (at - AI.time > prep.delay,1,0)
        c3 <- ifelse (at - AI.time > prep.delay,1,0)
      } else if (prep.elig.model == "adol.riskhist.older.time") {
        c1 <- ifelse (uaicount> prep.uaicount.thresh,1,0)
        c2 <- of.age
        c3 <- ifelse (at - AI.time > prep.delay,1,0)
      }
      
      
      idsEligStart <- which(c1 > 0 & c2 > 0 & c3 > 0)
      idsEligStop <- which(c1 == 0 | c2 == 0 | c3 == 0)
    } 
    
    prepElig[idsEligStart] <- 1
    prepElig[idsEligStop] <- 0  
  }

  

  ## Stoppage ------------------------------------------------------------------

  # Diagnosis
  idsStpDx <- which(active == 1 & prepStat == 1 & diag.status == 1)

  # Death
  idsStpDth <- which(active == 0 & prepStat == 1)
  
  #Drop out
  idsStpDrop<-integer(0)
  if (prepSpell == TRUE){
  idsStpDrop <-which(active==1 & prepStat == 1)
  
  if (length(idsStpDrop>=1)){
          idsStpDrop <-ssample(idsStpDrop, prepDrop*length(idsStpDrop), replace=FALSE)
  }
  }
  # Transition to ineligibility
  idsStpInelig <- idsEligStop

  # Reset PrEP status
  idsStp <- c(idsStpDx, idsStpDth, idsStpDrop, idsStpInelig)
  prepStat[idsStp] <- 0
  prepClass[idsStp] <-NA
  
  #Drops are added back to eligible list.
  prepElig[idsStpDrop] <- 1


  ## Initiation ----------------------------------------------------------------

  if (prep.cov.method == "curr") {
    prepCov <- sum(prepStat == 1, na.rm = TRUE)/sum(prepElig == 1, na.rm = TRUE)
  }
  if (prep.cov.method == "ever") {
    prepCov <- sum(prepEver == 1, na.rm = TRUE)/sum(prepElig == 1, na.rm = TRUE)
  }
  prepCov <- ifelse(is.nan(prepCov), 0, prepCov)

  idsEligSt <- which(prepElig == 1)
  nEligSt <- length(idsEligSt)

  nStart <- max(0, min(nEligSt, round((prep.coverage - prepCov) *
                                             sum(prepElig == 1, na.rm = TRUE))))
  idsStart <- NULL
  if (nStart > 0) {
    if (prep.cov.rate >= 1) {
      idsStart <- ssample(idsEligSt, nStart)
    } else {
      idsStart <- idsEligSt[rbinom(nStart, 1, prep.cov.rate) == 1]
    }
  }

  # Attributes
  if (length(idsStart) > 0) {
    prepStat[idsStart] <- 1
    prepEver[idsStart] <- 1
    prepStart.time[idsStart]<-at

    # PrEP class is fixed over PrEP cycles
    needPC <- which(is.na(prepClass[idsStart]))
    prepClass[idsStart[needPC]] <- sample(x = c("n", "l", "m", "h"), size = length(needPC),
                                          replace = TRUE, prob = prep.class.prob)
  }


  ## Output --------------------------------------------------------------------

  # Attributes
  dat$attr$prepElig <- prepElig
  dat$attr$prepStat <- prepStat
  dat$attr$prepEver <- prepEver
  dat$attr$prepClass <- prepClass
  dat$attr$prepStart.time<-prepStart.time
  
  # Summary Statistics
  dat$epi$prepCov[at] <- prepCov
  dat$epi$prepStart[at] <- length(idsStart)
 

  return(dat)
}
