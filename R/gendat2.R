#' simulation data generating function (Adding a measured confounder)
#'
#' generating the simulation data to apply in \code{brace}
#' @param nsims number of simulation datasets
#' @param nobs number of observations for one dataset
#' @param f1 parameter for generating unmeasured binary confounder
#' @param f2 parameter for generating measured binary confounder
#' @param g parameter for generating group assignment
#' @param b1 unmeasured confounder effect on group assignment
#' @param b2 measured confounder effect on group assignment
#' @param w1 shape parameter in generating survival time for event 1
#' from weibull distribution
#' @param w2 shape parameter in generating survival time for event 2
#' from weibull distribution
#' @param omegaplus multiplier on the baseline hazard for event 1
#' @param theta1 multiplier on the baseline hazard for event 1
#' @param theta2 multiplier on the baseline hazard for event 2
#' @param k3 multiplier on the baseline hazard for event 2
#' @return a matrix of nsims*nobs row, which consists of nsims datasets
#' @examples
#' nsims = 1; nobs = 1500
#' f1 = f2 = 0.5; g = 0.333; b1 = 8; b2 = 4; w1 = w2 = 0.667
#' theta1 = 0.5; theta2 = 1; omegaplus = 1; k3 = 0.333
#' sim1 = gendat2(nsims,nobs,f1,f2,g,b1,b2,w1,w2,omegaplus,theta1,theta2,k3)
#' @export
gendat2 <- function(nsims,nobs,f1,f2,g,b1,b2,w1,w2,omegaplus,theta1,theta2,k3){
  dat_output = matrix(NA, nrow = nobs, ncol = 17)
  for(i in 1:nsims){
    for(j in 1:nobs){
      # factor1 and factor2 are confounders
      factor1 = rbinom(1, 1, f1) # try f = 0.5 or 0.3
      factor2 = rbinom(1, 1, f2)

      # group assignment depends on value of factors and g and b parameters
      # b1, b2 are confounder effects on group assignment, b1 unobserved
      group = rbinom(1, 1, 1/(1+exp(-(log(g/(1-g))+log(b1)*factor1+
                                        log(b2)*factor2))))

      # theta1 and omegaplus are multipliers on the baseline hazard for event 1
      # theta2 and k3 are multipliers on the baseline hazard for event 2
      scale1 = (omegaplus*(1-group*(1-theta1)))^(-1/w1)
      scale2 = ((1-group*(1-theta2))*(1-factor1*(1-k3)))^(-1/w2)

      t1 = rweibull(1, w1, scale1)
      t2 = rweibull(1, w2, scale2)
      time = min(t1, t2)

      if(factor1 == 0 & group == 0 & (t1 < t2))
        report = 1
      else if(factor1 == 0 & group == 0 & (t1 >= t2))
        report = 2
      else
        report = 0

      event1 = ifelse(t1 < t2, 1, 0)
      event2 = ifelse(t1 >= t2, 1, 0)

      # impose random censoring
      censored = ifelse(runif(1) > .66, 1, 0)
      event1 = ifelse(censored, 0, event1)
      event2 = ifelse(censored, 0, event2)

      event = ifelse(event1 == 1 | event2 == 1, 1, 0)
      pfs_ci = 0
      if(event1 == 1)
        pfs_ci = 1
      else if(event2 == 1)
        pfs_ci = 2

      #scale and map event times;
      time = round(6000*time)/100 + 6;

      trueomega = omegaplus/(omegaplus+1-f1+f1*k3)
      truetheta = (theta1*omegaplus+theta2-theta2*f1+theta2*f1*k3)/(omegaplus+1-f1+f1*k3)

      r = i
      dat = c(r, factor1, factor2, group, scale1, scale2, t1, t2, time, report,
                  event1, event2, censored, event, pfs_ci, trueomega, truetheta)

      dat_output[j,] = dat
    }
  }

  dat_output = data.frame(dat_output)
  colnames(dat_output) = c("r", "factor1", "factor2", "group", "scale1", "scale2", "t1", "t2",
                           "time", "report", "event1", "event2", "censored",
                           "event", "pfs_ci", "trueomega", "truetheta")
  return(dat_output)
}
