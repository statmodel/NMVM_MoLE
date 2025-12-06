Stop.rule <- function(log.lik) {
  if (length(log.lik) >= 3) { 
    n = length(log.lik)
    l.new = log.lik[n]
    l.old  = log.lik[(n-1)]
    l.old2  = log.lik[(n-2)]
    ait = (l.new - l.old)/(l.old - l.old2)
    ln.Inf = l.old + (l.new - l.old)/(1 - ait)
    out = ln.Inf - l.new	
    if (!is.na(out) ) out = abs(out)	
    else out = 0
    
  } else out = (log.lik[2] - log.lik[1])/abs(log.lik[1])		
  return( out )
}


mix.sample = function(n, param, family =  c("BS", "Lindly", "gig", "Exp"))
{
  
  if(family == "BS") {
    alpha = param[1]
    U = rnorm(n)
    out = (alpha * U + sqrt(( alpha * U)^2 + 4))^2 / 4
  }
  
  if(family == "Lindly") {
    alpha = param[1]
    out = VGAM::rlind(n, alpha)
  }
  
  if(family == "gig") {
    kappa = param[1]
    chi = param[2]
    psi = param[3]
    out = GIGrvg::rgig(n, lambda = kappa, chi = chi, psi = psi)
  }
  
  if(family == "Exp") {
    out = rexp(n, 0.5)
  }
  return(out)
}

r.NMV.family = function(n, sigma2, lambda, theta,
                        family = c("normal", "SL", "GHST", 
                                   "VG", "NMVBS", "NMVL",
                                   "NIG", "GH"),
                        SN = FALSE){
  z = rnorm(n, mean = 0, sd = sqrt(sigma2))
  if(family == "SL") W = rexp(n, 0.5)
  if(family == "GHST") W = mix.sample(n, c(-theta/2, theta, 0), family =  "gig")
  if(family == "VG")  W = mix.sample(n, c(theta[1],0, theta[2]), family =  "gig")
  if(family == "NIG") W = mix.sample(n, c(-0.5, 1, theta^2), family =  "gig")
  if(family == "GH")  W = mix.sample(n, c(theta,theta,theta), family =  "gig")
  if(family == "NMVBS") W = mix.sample(n, theta, family =  "BS")
  if(family == "NMVL") W = mix.sample(n, theta, family =  "Lindly")
  if(family == "normal") W = 1; lambda = 0
  out = lambda * W + sqrt(W) * z
  if(SN){
    z = sn :: rsn(n, xi = 0, omega = 1, alpha = 1, tau = 0)
    out = lambda * W + sqrt(W) * z
  }
  return(out)
}


r.NMV.reg = function(x,beta, sigma2, lambda, theta,
                      family = c("normal", "SL", "GHST", 
                                 "VG", "NMVBS", "NMVL",
                                 "NIG", "GH"),
                      SN = FALSE){
  epsi = r.NMV.family(1, sigma2, lambda, theta, family = family, SN = SN)
  epsi = x%*%beta+epsi
  out = epsi
  return(out)
}

r.SMSN.reg= function(x,beta, sigma2, lambda, nu,
                       family = c("ESN", "ST" , "SSL" , "SCN")){
  p=length(beta)
  if(family == "ESN"){
      epsi = sn :: rsn(n = 1, xi = x%*%beta, omega = sigma2, alpha = lambda, tau = nu)
  }
  
  if(family == "ST"){
      epsi = sn :: rst(n = 1, xi = x%*%beta, omega = sigma2, alpha = lambda, nu = nu)
  }
  
  if(family == "SCN"){
      epsi = sn :: rsn(n = 1, xi = x%*%beta, omega = sigma2, alpha = lambda, tau=0)
      U=1
      if(runif(1,0,1)<nu[1]) U=nu[2]
      epsi = epsi*sqrt(1/U)
  }
  if(family == "SSL"){
      epsi = sn :: rsn(n = 1, xi = x%*%beta, omega = sigma2, alpha = lambda, tau=0)
      U=rbeta(1,nu,1)
      epsi = epsi*sqrt(1/U)
  }
  
  out =epsi
  return(out)
}


r.mix.reg.NMV = function(n,x,beta, sigma2, lambda, theta, alpha,r,
                          family = c("normal", "SL", "GHST", 
                                     "VG", "NMVBS", "NMVL",
                                     "NIG", "GH"),
                          SN = FALSE){
  Zi =  c()
  y = c()
  for(i in 1 : n){
    pI=c(exp(t(alpha)%*%r[i,])/(1+sum(c(exp(t(alpha)%*%r[i,])))))
    pI=c(pI,c(1/(1+sum(c(exp(t(alpha)%*%r[i,]))))))
    ni = c(rmultinom(1, 1, pI))
    j = Zi[i] = which.max(ni)
    ss=r.NMV.reg(x[i,],beta[,j], sigma2[j], lambda[j], theta[,j],family = family,SN = SN)
    y[i] =ss
      }
  out = list(y = y, class = Zi)
  return(out)
}

r.mix.reg.SMSN = function(n,x,beta, sigma2, lambda, nu, alpha,r,
                           family = c("ESN", "ST","SSL","SCN")){
  Zi =  c()
  y = c()
  for(i in 1 : n){
    pI=c(exp(t(alpha)%*%r[i,])/(1+sum(c(exp(t(alpha)%*%r[i,])))))
    pI=c(pI,c(1/(1+sum(c(exp(t(alpha)%*%r[i,]))))))
    ni = c(rmultinom(1, 1, pI))
    j = Zi[i] = which.max(ni)
    ss=r.SMSN.reg(x[i,],beta[,j], sigma2[j], lambda[j], nu[,j],family = family)
    y[i] = ss
  }
  out = list(y = y, class = Zi)
  return(out)
}

inv.mat = function(M){
  eg = eigen(M)
  val = diag(1/eg$val)
  vec = cbind(eg$vec)
  vec %*% val %*% solve(vec)
}



predict.mix.reg.NMV = function(n,x,beta, sigma2, lambda, theta, alpha,r,
                                family = c("normal", "SL", "GHST", 
                                           "VG", "NMVBS", "NMVL",
                                           "NIG", "GH"),
                                SN = FALSE){
  
  yhat=Zi=c()
  for(i in 1 : n){
    pI=c(exp(t(alpha)%*%r[i,])/(1+sum(c(exp(t(alpha)%*%r[i,])))))
    pI=c(pI,c(1/(1+sum(c(exp(t(alpha)%*%r[i,]))))))
    ni = c(rmultinom(1, 1, pI))
    j = Zi[i] = which.max(ni)
    epsi = r.NMV.family(1, sigma2[j], lambda[j], theta[,j], family = family, SN = SN)
    yhat[i] = x[i,]%*%beta[,j]+epsi
  }
  out = list(yhat = yhat, class = Zi)
  
  return(out)
}


predict.mix.reg.SN = function(n,x,beta, sigma2, lambda, theta, alpha,r,
                                family = c("ESN", "ST" , "SSL","SCN")){
  
  yhat=Zi=c()
  for(i in 1 : n){
    pI=c(exp(t(alpha)%*%r[i,])/(1+sum(c(exp(t(alpha)%*%r[i,])))))
    pI=c(pI,c(1/(1+sum(c(exp(t(alpha)%*%r[i,]))))))
    ni = c(rmultinom(1, 1, pI))
    j = Zi[i] = which.max(ni)
      if(family == "ESN") epsi = sn :: rsn(n = 1, xi = 0, omega = sigma2[j], alpha = lambda[j], tau = theta[j])
      if(family == "ST") epsi = sn :: rst(n = 1, xi = 0, omega = sigma2[j], alpha = lambda[j], nu = theta[j])
      if(family == "ESL") epsi = sn :: rsn(n = 1, xi = 0, omega = sigma2[j], alpha = lambda[j], tau=0)*sqrt(1/rbeta(1,theta[j],1))
      if(family == "ECN") {
        U=1
        if(runif(1,0,1)<theta[j][1]) U=theta[j][2]
        epsi = sn :: rsn(n = 1, xi = 0, omega = sigma2[j], alpha = lambda[j], tau=0)*sqrt(1/U)
      }
      yhat[i]  = x[i,]%*%beta[,j]+epsi
    }
  out = list(yhat = yhat, class = Zi)
  return(out)
}



predict2.mix.reg.NMV = function(n,x,beta, sigma2, lambda, theta, class,
                               family = c("normal", "SL", "GHST", 
                                          "VG", "NMVBS", "NMVL",
                                          "NIG", "GH"),
                               SN = FALSE){
  
  yhat=c()
  for(i in 1 : n){
    epsi = r.NMV.family(1, sigma2[class[i]], lambda[class[i]], theta[,class[i]], family = family, SN = SN)
    yhat[i] = x[i,]%*%beta[,class[i]]+epsi
  }
  out = yhat
  return(out)
}


predict2.mix.reg.SN = function(n,x,beta, sigma2, lambda, theta, class,
                              family = c("ESN", "ST" , "SSL","SCN")){
  
  yhat=c()
  for(i in 1 : n){
    if(family == "ESN") epsi = sn :: rsn(n = 1, xi = 0, omega = sigma2[class[i]], alpha = lambda[class[i]], tau = theta[class[i]])
    if(family == "ST") epsi = sn :: rst(n = 1, xi = 0, omega = sigma2[class[i]], alpha = lambda[class[i]], nu = theta[class[i]])
    if(family == "ESL") epsi = sn :: rsn(n = 1, xi = 0, omega = sigma2[class[i]], alpha = lambda[class[i]], tau=0)*sqrt(1/rbeta(1,theta[class[i]],1))
    if(family == "ECN") {
      U=1
      if(runif(1,0,1)<theta[class[i]][1]) U=theta[class[i]][2]
      epsi = sn :: rsn(n = 1, xi = 0, omega = sigma2[class[i]], alpha = lambda[class[i]], tau=0)*sqrt(1/U)
    }
    yhat[i]  = x[i,]%*%beta[,class[i]]+epsi
  }
  out = yhat
  return(out)
}
