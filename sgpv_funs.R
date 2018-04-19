### FUNCTIONS ##################################################
sgpvFun <- function(ci.lb, ci.ub, null.lb, null.ub, return_dg = TRUE) {
  ## updated 1/5/2018: can now take vectors of confidence interval bounds (but single null interval bounds)
  
  # ci.lb <- pull(out_reg, ci.lb)
  # ci.ub <- pull(out_reg, ci.ub)
  
  length.ci = ci.ub-ci.lb
  length.null = null.ub-null.lb
  
  pval_delta = sapply(1:length(ci.lb), function(i) {(max(null.lb,min(null.ub,ci.ub[i]))-min(null.ub,max(null.lb,ci.lb[i])))/length.ci[i]*max(length.ci[i]/(2*length.null),1)})
  
  if(length.null == 0) {  ## alternately, if null.lb == null.ub
    
    null <- null.lb
    # null <- null.ub    ## if length.null == 0, this means that null.lb == null.ub, thus null = null.lb = null.ub
    
    pval_delta <- sapply(1:length(ci.lb), function(i) {if(null >= ci.lb[i] & null <= ci.ub[i]) {pd <- 0.5
    } else if(null <= ci.lb[i] | null >= ci.ub[i]) {pd <- 0};
      return(pd)})
    
    
  }
  
  
  if(return_dg == TRUE) {
    
    deltagap <- sapply(1:length(ci.lb), function(i) {
      if(is.na(pval_delta[i])) {
        dg = NA
      } else if(pval_delta[i] > 0) {
        dg = NA
      } else if(pval_delta[i] == 0) {
        dg = (max(ci.lb[i],null.lb)-min(null.ub,ci.ub[i]))/((null.ub-null.lb)/2)
        if(length.null == 0) {
          null <- null.lb
          if(null <= ci.lb[i]) {
            dg = ci.lb[i] - null
          } else if(null >= ci.ub[i]) {
            dg = null - ci.ub[i]
          }
        }
      } 
      return(dg)
    })
    
    out <- list(sgpv = pval_delta, delta.gap = deltagap)
  } else{out <- c(sgpv = pval_delta)}
  
  return(out)
  
}
sgpvPowerFun <- function(beta, beta.se, type = c("d", "c"), null.lb, null.ub) {
  
  if(type == "d") {
    pnorm((null.lb)/beta.se - beta/beta.se - qnorm(1-0.05/2)) + pnorm(-(null.ub)/beta.se + beta/beta.se - qnorm(1-0.05/2))
  } else if(type == "c") {
    pnorm((null.ub)/beta.se - beta/beta.se - qnorm(1-0.05/2)) - pnorm((null.lb)/beta.se - beta/beta.se + qnorm(1-0.05/2))
  }
  
}
sgpvFdrFun <- function(beta, beta.se, pi0 = 0.5, null, null.lb, null.ub) {
  
  piA = 1 - pi0
  
  powerA = sgpvPowerFun(beta, beta.se, type = "d", null.lb = null.lb, null.ub = null.ub)
  power0 = sgpvPowerFun(null, beta.se, type = "d", null.lb = null.lb, null.ub = null.ub)
  
  fdr = (1 + powerA/power0 * piA/pi0)^(-1)
  
}
sgpvFcrFun <- function(beta, beta.se, pi0 = 0.5, null, null.lb, null.ub) {
  
  piA = 1 - pi0
  
  powerA = sgpvPowerFun(beta, beta.se, type = "c", null.lb = null.lb, null.ub = null.ub)
  power0 = sgpvPowerFun(null, beta.se, type = "c", null.lb = null.lb, null.ub = null.ub)
  
  fcr = (1 + power0/powerA * pi0/piA)^(-1)
  
}
################################################################
