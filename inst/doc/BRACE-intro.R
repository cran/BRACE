## ---- include = FALSE, echo = FALSE-------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(knitr)
library(BRACE)

## ----echo = FALSE-------------------------------------------------------------
nsims = 1; nobs = 1500
f1 = f2 = 0.5; g = 0.333; b1 = 8; b2 = 4; w1 = w2 = 0.667
theta1 = 0.5; theta2 = 1; omegaplus = 1; k3 = 0.333
sim1 = gendat2(nsims,nobs,f1,f2,g,b1,b2,w1,w2,omegaplus,theta1,theta2,k3)
ftime = sim1$time
fstatus = sim1$pfs_ci
covs = sim1$factor2
trt = sim1$group

## ----echo = TRUE--------------------------------------------------------------
braceoutput = brace(ftime, fstatus, covs, trt, PS=0, B=200)
braceoutput$Summary
braceoutput$`Omega Curve`

## ----echo = TRUE--------------------------------------------------------------
braceoutput = brace(ftime, fstatus, covs, trt, PS=1, B=200)
braceoutput$Summary
braceoutput$`Omega Curve`

