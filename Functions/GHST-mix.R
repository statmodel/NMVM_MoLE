
mix.reg.ST.EM <- function(y, x, beta=NULL, sigma2 = NULL,
                                 lambda = NULL, nu = NULL,  r, alpha=NULL,
                                 g = NULL, Class = NULL, p=NULL, error = 0.00001,
                                 iter.max = 100, Stp.rule = c("Log.like", "Atiken"), 
                                 error.est = F, per = 1, print = T, fix.sigma = F)
{
  
  BesselK = function(x, kappa, log.val = T){
    F1 = log( besselK(x, kappa, expon.scaled = TRUE) ) - x
    F2 = besselK(x, kappa)
    ifelse(log.val == T, return(F1), return(F2))
  }
  
  f.ST = function(x, sigma2, lambda, nu, log.val = F){
    ob = ghyp :: ghyp(lambda = -nu/2, chi = nu, psi = 0,
                      sigma = sqrt(sigma2), gamma = lambda, data = NULL)
    Val= ghyp :: dghyp(x, object = ob, logvalue = log.val)
    Val
  }
  d.ST.mix = function(y, x, beta, sigma2, lambda, nu, r, alpha,class,log = F)
  {
    g = length(sigma2)
    PDF = 0
    for(j in 1:g){
      Cent=y[class==j]-x[class==j,]%*%beta[,j]
      rnew=r[class==j,]
      for(i in 1:length(Cent)){
        if(j!=g) pI=c(exp(rnew[i,]%*%alpha[,j])/(1+sum(c(exp(t(alpha)%*%rnew[i,])))))
        if(j==g) pI=c(1/(1+sum(c(exp(t(alpha)%*%rnew[i,])))))
        PDF = PDF + pI * f.ST(Cent[i], sigma2[j], lambda[j], nu[j], log.val = F)
      }
    }
    Out = sum(PDF)
    ifelse(log == T, return(log(Out)), return(Out))
  }
  
  if(isTRUE(print)){
    cat(paste(rep("-", 70), sep = "", collapse = ""), "\n")
    cat('Finite mixture of GHST-based reg model','\n')
  }
  if(is.null(sigma2) || is.null(alpha) || is.null(beta) || is.null(lambda)){
    class = Class
    if(is.null(g) && is.null(class) ) 
      stop("The model is not specified correctly.\n")
    if(!is.null(class)){
      tt = table(Class)
      if(g == 0 || max(as.numeric(labels(tt)$Class)) != g) 
        g = max(as.numeric(labels(tt)$Class))
    }
    
    
    if(is.null(class)) {
      minMaxScaler <- function(data) {  
        return ((data - min(data)) / (max(data) - min(data)))  
      }
      class= kmeans(minMaxScaler(y),g)[[1]]
    }
    library(moments)
    beta=matrix(NA,p,g)
    sigma2=lambda=rep(0,0)
    for(j in 1:g){
      beta[,j]=solve(t(x[class==j,])%*%x[class==j,])%*%t(x[class==j,])%*%y[class==j]
      sigma2[j]=var(y[class==j]-x[class==j,]%*%beta[,j])
    }
    pI = table(class)/n
    CC=1/pI[g]
    DD=log(pI*CC)[1]
    rstar=apply(r,2,mean)
    q=length(rstar)
    al<-function(par){
      ll=par[1]
      for(i in 2:(q-1)){
        ll=ll+par[i]*rstar[i]
      }
      ll-DD
    }
    opt<- nlminb(rep(1,q), al, lower = rep(0.1, q), 
                 upper = rep(5, q), control = list(iter.max = 2000))
    alpha=opt$par
    alpha=matrix(rep(alpha,g-1),ncol=(g-1))
    nu = rep(10, g)
    lambda=rep(2,g)
  }else{
    class = Class
    if(is.null(g) && is.null(class) ) 
      stop("The model is not specified correctly.\n")
    if(!is.null(class)){
      tt = table(Class)
      if(g == 0 || max(as.numeric(labels(tt)$Class)) != g) 
        g = max(as.numeric(labels(tt)$Class))
    }
    
    if(is.null(class)) {
      minMaxScaler <- function(data) {  
        return ((data - min(data)) / (max(data) - min(data)))  
      }
      class= kmeans(minMaxScaler(y),g)[[1]]
    }
    library(moments)
    beta=beta
    sigma2=sigma2
    lambda=lambda
    pI = table(class)/n
    CC=1/pI[g]
    alpha=alpha
  }
  
  g = length(pI)
  
  ##------------------------------------------------------##
  ##----------------- Expectations -----------------------##
  ##------------------------------------------------------##
  EUY.ST = function(y,x,beta, sigma2, lambda, nu)
  {
    R = function(cc, k, a) BesselK(cc, k+a, log.val = T) - BesselK(cc, k, log.val = T) 
    erii=x%*%beta
    chii = nu + (y-erii)^2 / sigma2
    psii = lambda^2 / sigma2
    kappa = -(nu+1)/2
    alphai = sqrt(chii * psii)
    w = exp(0.5 * log(chii / psii) + R(alphai, kappa, 1))
    t = exp(-0.5 * log(chii / psii) + R(alphai, kappa, -1))
    return(list(w = w, t = t))
  }
  
  start.time = Sys.time()
  lk = lk.old = sum(d.ST.mix(y, x, beta, sigma2, lambda, nu, r, alpha,class, log = T))
  criterio = 1;  count = 0
  if(isTRUE(print)){
    cat(paste(rep("-", 70), sep = "", collapse = ""), "\n")
    cat("iter =", count, "\t logli.old=",  lk.old, "\n")
    cat(paste(rep("-", 70), sep = "", collapse = ""), "\n")
  }
  library(optimx)
  repeat
  {
    count = count + 1
    tal = matrix(0, n, g)
    for(j in 1:g){
      for(i in which(class==j)){
        e=x[i,]%*%beta[,j]
        pI=c(exp(t(alpha)%*%r[i,])/(1+sum(c(exp(t(alpha)%*%r[i,])))))
        pI=c(pI,c(1/(1+sum(c(exp(t(alpha)%*%r[i,]))))))
        tal[i, j] = pI[j] * f.ST(y[i]-e, sigma2[j], lambda[j], nu[j], log.val = F)
      }
    }
    for(k in 1:n) if(all(tal[k,] == 0)) tal[k,] = .Machine$double.xmin
    tal = tal/rowSums(tal)
    class=apply(tal, 1, which.max)
    SS = 0
    for (j in 1:g){
      u=w=rep(0,0)
      for(i in which(class==j)){
        ### E-step: 
        Ewt = EUY.ST(y[i], x[i,],beta[,j],sigma2[j], lambda[j], nu[j])
        w = c(w,Ewt$w); u = c(u,Ewt$t)
      }
      
      ### M-step: 
      beta[,j]=solve(t(u*x[class==j,])%*%x[class==j,])%*%t(x[class==j,])%*%(u*y[class==j]-lambda[j])
      
      Cent=c(y[class==j]-beta[,j]%*%t(x[class==j,]))
      lambda[j]=sum(Cent)/ sum(w)
      
      bb = u *(Cent)^2 + (lambda[j])^2 * w - 2 * lambda[j] * (Cent)
      sigma2[j] =sum( bb ) / sum(tal[, j])
      if(fix.sigma) SS = sum(bb) + SS
      
      ss=0
      for(i in 1:n) ss=ss+(tal[i, j]-c(exp(t(alpha)%*%r[i,])/(1+sum(c(exp(t(alpha)%*%r[i,])))))[j])*r[i,]
      if(j!=g) alpha[,j]=((4/sum(diag(r%*%t(r))))*ss)+alpha[,j]
    }
    ft = function(nu) -sum(d.ST.mix(y, x, beta, sigma2, lambda, nu, r, alpha,class, log = T))
    opt = nlminb(nu, ft, lower = rep(0.01,g), upper = runif(g,1,100), 
                 control = list(iter.max = 2000))
    nu = opt$par
    
    lk.new = sum(d.ST.mix(y, x, beta, sigma2, lambda, nu, r, alpha,class, log = T))
    
    if(is.nan(lk.new)) {
      lk.new = lk.old
      break
    }
    
    lk = c(lk, lk.new)
    if(Stp.rule == "Log.like") criterio = (lk.new - lk.old)/abs(lk.old)
    else { criterio = Stop.rule(lk) }
    
    diff = lk.new -lk.old
    if(count %% per == 0 || is.na(diff))
    {
      if(isTRUE(print)){
        cat('iter =', count, '\t logli =', lk.new, '\t diff =', 
            diff, Stp.rule, "'s diff =", criterio, '\n')
        cat(paste(rep("-", 60), sep = "", collapse = ""), "\n")
      }
    }
    
    if(criterio < error | count == iter.max) break
    lk.old = lk.new
  }
  
  # End of the estimation process
  lk = lk.new
  end.time = Sys.time()
  time.taken = end.time - start.time
  
  m <- g*(nrow(beta)+3)+nrow(alpha)*(g-1) 
  if(fix.sigma) m = m - g + 1
  
  aic = -2 * lk + 2 * m
  bic = -2 * lk + log(n) * m
  edc = -2 * lk + 0.2 * sqrt(n) * m
  aic_c = -2 * lk + 2 * n * m / (n - m - 1)
  abic = -2 * lk + m * log((n + 2) / 24)
  
  Group = apply(tal, 1, which.max)
  
  if(!is.null(Class)){
    true.clus = Class
    km.clus = Group
    tab = table(true.clus, km.clus)
    MCR = 1 - sum(diag(tab))/sum(tab)
    RII = aricode :: clustComp(km.clus, true.clus)
    
    obj.out = list( time = time.taken, group = Group, m = m,
                    sigma2 = sigma2, lambda = lambda, nu = nu, beta = beta,alpha=alpha, 
                    loglike = lk, aic = aic, bic = bic, edc = edc, 
                    aic_c = aic_c, abic = abic, iter = count, MCR = MCR, 
                    similarity = t(as.matrix(RII)), 
                    cross_class = tab, Zij = tal,
                    convergence = criterio < error, crite = criterio)
    
  }
  else{obj.out = list( time = time.taken, group = Group, m = m,
                       sigma2 = sigma2, lambda = lambda, nu = nu, beta = beta, alpha=alpha, 
                       loglike = lk, aic = aic, bic = bic, edc = edc, 
                       aic_c = aic_c, abic = abic, iter = count, 
                       convergence = criterio < error, crite = criterio)} 
  
  obj.out
}
