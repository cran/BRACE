nsims = 1; nobs = 1500
f = 0.5; g = 0.333; b = 8; w1 = w2 = 0.667
theta1 = 0.5; theta2 = 1; omegaplus = 1; k3 = 0.333
sim1 = gendat(nsims,nobs,f,g,b,w1,w2,omegaplus,theta1,theta2,k3)
ftime = sim1$time
fstatus = sim1$pfs_ci
covs = rnorm(nrow(sim1)) #simulate a column of uncorrelated covariate
trt = sim1$group
braceoutput = brace(ftime, fstatus, covs, trt, PS=0, B=200)
