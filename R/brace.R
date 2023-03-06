#' Bias Reduction through Analysis of Competing Events
#'
#' \code{brace} is used to estimate the treatment effect with adjusted
#' confounders on the composite hazard
#' for primary or competing events, and adjust for bias from residual
#' confounding in non-randomized data by BRACE method
#' @param ftime vector of failure/censoring times
#' @param fstatus vector with a unique code for each failure type and a
#' separate code for censored observations (default is primary event = 1,
#' competing event = 2, censored = 0)
#' @param covs matrix (nobs x ncovs) of fixed covariates. If no covariates,
#' set covs = NA (default is NA)
#' @param trt vector of treatment indicator (1 for treatment group)
#' @param failcode code of fstatus that denotes the failure type of interest
#' @param cencode code of fstatus that denotes censored observations
#' @param PS whether to use propensity score method for adjusting the
#' confounding effect (1 for propensity score method, default is 0)
#' @param B bootstrap sample size for calculating the Confidence interval,
#' default is 1000
#' @return a list of class \code{brace}, with components:
#' \item{$Summary}{summary table of BRACE method}
#' \item{$`BRACE HR Distribution`}{the estimated regression coefficients in
#' each bootstrap sample}
#' \item{$`Omega Estimate`}{estimate of relative hazards for
#' primary events vs. combined events}
#' \item{$Epsilon}{the estimated bias}
#' \item{$`Combined Endpoint Model`}{the regression model for combined events}
#' \item{$`Primary Endpoint Model`}{the regression model for primary events}
#' \item{$`Competing Endpoint Model`}{the regression model for competing events}
#' \item{$`Omega Curve`}{estimate of omega over time}
#' \item{$`Combined Endpoint Curve`}{survival curve for combined events}
#' \item{$`Primary Endpoint Curve`}{survival curve for primary events}
#' \item{$`Competing Endpoint Curve`}{survival curve for competing events}
#' @import survival,survminer,ipw
#' @references Williamson, Casey W., et al. "Bias Reduction through Analysis of
#' Competing Events (BRACE) Correction to Address Cancer Treatment Selection
#' Bias in Observational Data." Clinical Cancer Research 28.9 (2022): 1832-1840.
#' @examples
#' nsims = 1; nobs = 1500
#' f = 0.5; g = 0.333; b = 8; w1 = w2 = 0.667
#' theta1 = 0.5; theta2 = 1; omegaplus = 1; k3 = 0.333
#' sim1 = gendat(nsims,nobs,f,g,b,w1,w2,omegaplus,theta1,theta2,k3)
#' ftime = sim1$time
#' fstatus = sim1$pfs_ci
#' covs = NA
#' trt = sim1$group
#' braceoutput = brace(ftime, fstatus, covs, trt, PS=0, B=10)
#'
#' nsims = 1; nobs = 1500
#' f1 = f2 = 0.5; g = 0.333; b1 = 8; b2 = 4; w1 = w2 = 0.667
#' theta1 = 0.5; theta2 = 1; omegaplus = 1; k3 = 0.333
#' sim1 = gendat2(nsims,nobs,f1,f2,g,b1,b2,w1,w2,omegaplus,theta1,theta2,k3)
#' ftime = sim1$time
#' fstatus = sim1$pfs_ci
#' covs = sim1$factor2
#' trt = sim1$group
#' braceoutput = brace(ftime, fstatus, covs, trt, PS=1, B=10)
#' @export
brace <- function(ftime, fstatus, covs = NA, trt, failcode = 1, cencode = 0,
                  PS = 0, B = 1000){
  if(all(is.na(covs)) == 0){
    covs = data.frame(covs)
    if(!var(c(length(ftime), length(fstatus), nrow(covs), length(trt))) == 0){
      stop("Make sure input dimension matches")
    }
  }else if(PS == 1){
    stop("covariates are needed for propensity model")
  }else{
    covs = NULL
  }

  ## create dummy variables
  cmpcode = 3-failcode-cencode
  prim = ifelse(fstatus == failcode, 1, 0)
  cmp = ifelse(fstatus == cmpcode, 1, 0)

  ## combine the data
  dat_comb = data.frame(cbind(ftime, prim, cmp, trt, covs))
  colnames(dat_comb) = c("time", "primary", "competing", "treatment", colnames(covs))
  dat_comb$combine = ifelse(dat_comb$primary == 1 | dat_comb$competing == 1, 1, 0)
  dat_comb$treatment = as.factor(dat_comb$treatment)

  ## run the multiple regression model for combined endpoint
  if(PS == 0){
    if(is.null(covs) == 0){
      fml_combine = as.formula(paste0("Surv(time, combine) ~ treatment + ",
                                      paste(colnames(covs), collapse = "+")))
    }else{
      fml_combine = as.formula(paste0("Surv(time, combine) ~ treatment"))
    }

    fit_combine = coxph(fml_combine, dat_comb)
    theta = summary(fit_combine)$coef[1,2]
    thetase = theta*summary(fit_combine)$coef[1,3] # apply delta method
    #thetase = summary(fit_combine)$coef[1,3]
    thetaCI = exp(confint(fit_combine))[1,]

    comb_fit = surv_fit(Surv(time, combine) ~ treatment, dat_comb)
    comb_survcurv = ggsurvplot(comb_fit, risk.table=TRUE, conf.int = TRUE, pval = TRUE,
                               title=c('Univariable Treatment Effect on Combined Endpoint'))

    ## fit the baseline data
    dat_comb_base = subset(dat_comb, dat_comb$treatment == 0) #baseline data
    fit_base = coxph(Surv(time, combine) ~ 1, dat_comb_base)
    base_t = basehaz(fit_base)[,2]
    mean_t = mean(dat_comb_base$time)
    if(any(dat_comb_base$time < 0)){
      stop("Failure/censoring time should be non-negative")
    }else{
      ind = which.min(abs(base_t - mean_t)) #select the record with the time that is closest to the mean event time in the baseline data set
    }

    ## primary events
    if(is.null(covs) == 0){
      fml_primary = as.formula(paste0("Surv(time, primary) ~ treatment + ",
                                    paste(colnames(covs), collapse = "+")))
    }else{
      fml_primary = as.formula(paste0("Surv(time, primary) ~ treatment"))
    }
    fit_primary = coxph(fml_primary, dat_comb)
    theta1 = summary(fit_primary)$coef[1,2]
    theta1CI = exp(confint(fit_primary))[1,]

    prim_fit = surv_fit(Surv(time, primary) ~ treatment, dat_comb)
    prim_survcurv = ggsurvplot(prim_fit, risk.table=TRUE, conf.int = TRUE, pval = TRUE,
                               title=c('Univariable Treatment Effect on Primary Event'))

    ## competing events
    if(is.null(covs) == 0){
      fml_compete = as.formula(paste0("Surv(time, competing) ~ treatment + ",
                                      paste(colnames(covs), collapse = "+")))
    }else{
      fml_compete = as.formula(paste0("Surv(time, competing) ~ treatment"))
    }

    fit_compete = coxph(fml_compete, dat_comb)
    theta2 = summary(fit_compete)$coef[1,2]
    #theta2se = summary(fit_compete)$coef[1,3]
    theta2se = theta2*summary(fit_compete)$coef[1,3] # apply delta method
    theta2P = summary(fit_compete)$coef[1,5]
    theta2CI = exp(confint(fit_compete))[1,]

    comp_fit = surv_fit(Surv(time, competing) ~ treatment, dat_comb)
    comp_survcurv = ggsurvplot(comp_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE,
                               title=c('Univariable Treatment Effect on Competing Event'))

    ## produce error if theta2 not smaller than 1 significantly
    if(!(theta2 < 1 & theta2P < 0.05)){
      warning("No Significant Benefit on Competing Event from Treatment, BRACE Not Indicated")
    }

    ## estimate omega
    ev1 = (survfit(coxph(Surv(dat_comb_base$time, dat_comb_base$primary) ~ 1), ctype = 1))
    ev = (survfit(coxph(Surv(dat_comb_base$time, dat_comb_base$combine) ~ 1), ctype = 1))
    omega = ev1$cumhaz/ev$cumhaz
    omega_fig_dat = data.frame(cbind(omega), 1:length(ev1$cumhaz), ev1$time)
    colnames(omega_fig_dat) = c("omega", "time")
    omega_curve = ggplot(aes(x = time, y = omega), data = omega_fig_dat) +
      geom_line() +
      geom_point() +
      ggtitle("Omega over time")

    ev1.haz = ev1$cumhaz[ind]
    ev.haz = ev$cumhaz[ind]
    omega_base = ev1.haz/ev.haz
    # ev1.se = ev1$std.err[ind]
    # ev.se = ev$std.err[ind]

    ## Bootstrap to get the adjusted estimator
    thetaB <- rep(0, B)

    for(b in 1:B){
      set.seed(b)
      # thetaB[b] = theta + (1-rnorm(1,mean=ev1.haz,sd=ev1.se)/rnorm(1,mean=ev.haz,sd=ev.se))*
      #   (1-rnorm(1,mean=theta2,sd=theta2se))
      n = nrow(dat_comb)
      boot_ind = sample(1:n, size = n, replace = TRUE)
      dat_comb_boot = dat_comb[boot_ind,]
      dat_comb_boot_base = subset(dat_comb_boot, dat_comb_boot$treatment == 0)
      fit_boot_base = coxph(Surv(time, combine) ~ 1, dat_comb_boot_base)
      base_boot_t = basehaz(fit_boot_base)[,2]
      ind_boot = which.min(abs(base_boot_t - mean_t))

      ev1_boot = (survfit(coxph(Surv(dat_comb_boot_base$time, dat_comb_boot_base$primary) ~ 1), ctype = 1))
      ev_boot = (survfit(coxph(Surv(dat_comb_boot_base$time, dat_comb_boot_base$combine) ~ 1), ctype = 1))
      ev1_haz_boot = ev1_boot$cumhaz[ind_boot]
      ev_haz_boot = ev_boot$cumhaz[ind_boot]

      if(is.null(covs) == 0){
        fml_compete_boot = as.formula(paste0("Surv(time, competing) ~ treatment + ",
                                             paste(colnames(covs), collapse = "+")))
      }else{
        fml_compete_boot = as.formula(paste0("Surv(time, competing) ~ treatment"))
      }
      fit_compete_boot = coxph(fml_compete_boot, dat_comb_boot)
      theta2_boot = summary(fit_compete_boot)$coef[1,2]

      thetaB[b] = theta + (1-ev1_haz_boot/ev_haz_boot)*(1-theta2_boot)
    }

    theta_BRACE = mean(thetaB)
    adjlower95 = quantile(thetaB, 0.025)
    adjupper95 = quantile(thetaB, 0.975)
    epsilon = mean(theta_BRACE)-theta
  }

  if(PS == 1){
    ## run the propensity score method for combined endpoint
    fml_ps = as.formula(paste0("treatment ~ ", paste(colnames(covs), collapse = "+")))
    ps_model = glm(fml_ps, family = binomial, data = dat_comb)
    dat_comb$ps = ifelse(dat_comb$treatment == 1, ps_model$fitted.values, 1-ps_model$fitted.values)
    dat_comb$ps = pmin(pmax(dat_comb$ps, 0.01), 1-0.01)
    dat_comb$id = 1:nrow(dat_comb)

    fit_combine = coxph(Surv(time, combine) ~ treatment, weights = ps, dat_comb, cluster = id)
    theta = summary(fit_combine)$coef[1,2]
    thetaCI = exp(confint(fit_combine))[1,]

    comb_fit = surv_fit(Surv(time, combine) ~ treatment, data = dat_comb, weights = dat_comb$ps)
    comb_survcurv = ggsurvplot(comb_fit, risk.table=TRUE, conf.int = TRUE, pval = TRUE,
                                  title=c('Weighted Treatment Effect on Combined Endpoint'))

    ## fit the baseline data
    dat_comb_base = subset(dat_comb, dat_comb$treatment == 0) #baseline data
    fit_base = coxph(Surv(time, combine) ~ 1, dat_comb_base)
    base_t = basehaz(fit_base)[,2]
    mean_t = mean(dat_comb_base$time)
    if(any(dat_comb_base$time < 0)){
      stop("Failure/censoring time should be non-negative")
    }else{
      ind = which.min(abs(base_t - mean_t)) #select the record with the time that is closest to the mean event time in the baseline data set
    }

    ## primary events
    fit_primary = coxph(Surv(time, primary) ~ treatment, weights = ps, dat_comb, cluster = id)
    theta1 = summary(fit_primary)$coef[1,2]
    theta1CI = exp(confint(fit_primary))[1,]

    prim_fit = surv_fit(Surv(time, primary) ~ treatment, data = dat_comb, weights = dat_comb$ps)
    prim_survcurv = ggsurvplot(prim_fit, risk.table=TRUE, conf.int = TRUE, pval = TRUE,
                               title=c('Weighted Treatment Effect on Primary Event'))

    ## competing events
    fit_compete = coxph(Surv(time, competing) ~ treatment, weights = ps, dat_comb, cluster = id)
    theta2 = summary(fit_compete)$coef[1,2]
    #theta2se = theta2*summary(fit_compete)$coef[1,3] # apply delta method
    theta2P = summary(fit_compete)$coef[1,6]
    theta2CI = exp(confint(fit_compete))[1,]

    comp_fit = surv_fit(Surv(time, competing) ~ treatment, data = dat_comb, weights = dat_comb$ps)
    comp_survcurv = ggsurvplot(comp_fit, risk.table = TRUE, conf.int = TRUE, pval = TRUE,
                               title=c('Weighted Treatment Effect on Competing Event'))

    ## produce error if theta2 not smaller than 1 significantly
    if(!(theta2 < 1 & theta2P < 0.05)){
      warning("No Significant Benefit on Competing Event from Treatment, BRACE Not Indicated")
    }

    ## estimate omega
    ev1 = (survfit(coxph(Surv(dat_comb_base$time, dat_comb_base$primary) ~ 1), ctype = 1))
    ev = (survfit(coxph(Surv(dat_comb_base$time, dat_comb_base$combine) ~ 1), ctype = 1))
    omega = ev1$cumhaz/ev$cumhaz
    omega_fig_dat = data.frame(cbind(omega), 1:length(ev1$cumhaz), ev1$time)
    colnames(omega_fig_dat) = c("omega", "time")
    omega_curve = ggplot(aes(x = time, y = omega), data = omega_fig_dat) +
      geom_line() +
      geom_point() +
      ggtitle("Omega over time")

    ev1.haz = ev1$cumhaz[ind]
    ev.haz = ev$cumhaz[ind]
    omega_base = ev1.haz/ev.haz
    # ev1.se = ev1$std.err[ind]
    # ev.se = ev$std.err[ind]

    ## Bootstrap to get the adjusted estimator
    thetaB <- rep(0, B)

    for(b in 1:B){
      set.seed(b)
      n = nrow(dat_comb)
      boot_ind = sample(1:n, size = n, replace = TRUE)
      dat_comb_boot = dat_comb[boot_ind,]
      dat_comb_boot_base = subset(dat_comb_boot, dat_comb_boot$treatment == 0)
      fit_boot_base = coxph(Surv(time, combine) ~ 1, dat_comb_boot_base)
      base_boot_t = basehaz(fit_boot_base)[,2]
      ind_boot = which.min(abs(base_boot_t - mean_t))

      ps_model_boot = glm(fml_ps, family = binomial, data = dat_comb_boot)
      dat_comb_boot$ps = pmin(pmax(ps_model_boot$fitted.values, 0.01), 100)
      ev1_boot = (survfit(coxph(Surv(dat_comb_boot_base$time, dat_comb_boot_base$primary) ~ 1), ctype = 1))
      ev_boot = (survfit(coxph(Surv(dat_comb_boot_base$time, dat_comb_boot_base$combine) ~ 1), ctype = 1))
      ev1_haz_boot = ev1_boot$cumhaz[ind_boot]
      ev_haz_boot = ev_boot$cumhaz[ind_boot]

      fit_compete_boot = coxph(Surv(time, competing) ~ treatment, weights = ps, dat_comb_boot)
      theta2_boot = summary(fit_compete_boot)$coef[1,2]

      thetaB[b] = theta + (1-ev1_haz_boot/ev_haz_boot)*(1-theta2_boot)
    }

    theta_BRACE = mean(thetaB)
    adjlower95 = quantile(thetaB, 0.025)
    adjupper95 = quantile(thetaB, 0.975)
    epsilon = mean(theta_BRACE)-theta
  }



  ## output summary table
  Summary = data.frame(t(matrix(c(theta_BRACE, adjlower95, adjupper95,
                                theta, thetaCI,
                                theta1, theta1CI,
                                theta2, theta2CI), nrow = 3)))
  colnames(Summary) = c("HR Estimate", "Lower 95% CI", "Upper 95% CI")

  rownames(Summary) = c('BRACE Estimate', 'Adjusted Combined Event',
                        'Adjusted Primary Event', 'Adjusted Competing Event')

  output = list(Summary, thetaB, omega_base, epsilon, fit_combine, fit_primary, fit_compete,
                omega_curve, comb_survcurv, prim_survcurv, comp_survcurv)
  names(output) = c('Summary', 'BRACE HR Distribution', 'Omega Estimate', 'Epsilon',
                    'Combined Endpoint Model', 'Primary Endpoint Model','Competing Endpoint Model',
                    "Omega Curve" ,'Combined Endpoint Curve',
                    'Primary Endpoint Curve','Competing Endpoint Curve')
  return(output)
}
