
#' @title Sexual Debut Module
#'
#' @description Module for sexual debut over time for active nodes in the population.  Having debuted is a condition for AI but it does not
#' require that AI occur.  It is a debut in terms of identity and willingness not action.  
#'
#' @param dat Master data list object of class \code{dat} containing networks,
#'        individual-level attributes, and summary statistics.
#' @param at Current time step.
#'
#' @inheritParams aging.mard
#' 
#' @return
#' This function returns \code{dat} after updating the nodal attribute
#' \code{has.debuted}. The \code{has.debuted} vertex attribute is also
#' updated on the three networks.
#'
#' @keywords module
#' @export
#'

debut.adol <- function(dat, at) {
  

  debutElig<-which(dat$attr$active == 1 & dat$attr$has.debuted == 0)
  debut <-rbinom(length(debutElig),1,dat$init$debut.prob)
  debut.time<-debut
  debut.time<-ifelse(debut.time==1,at,0)
  

  dat$attr$has.debuted[debutElig]<-debut
  dat$attr$has.debuted.time[debutElig]<-debut.time

  has.debuted<-dat$attr$has.debuted
  

  
  
    for (i in seq_along(dat$nw)) {
      dat$nw[[i]] <- set.vertex.attribute(dat$nw[[i]], "has.debuted", has.debuted)
    
  }

  return(dat)
}



