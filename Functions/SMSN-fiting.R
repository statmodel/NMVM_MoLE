
dt.lsC <- function(x, loc, sigma2 = 1, shape = 1, nu = 4){
  k1 <- sqrt(nu / 2) * gamma((nu-1)/2) / gamma(nu/2)
  b <-  -sqrt(2/pi)*k1
  delta <- shape / (sqrt(1 + shape^2))
  Delta <- sqrt(sigma2)*delta
  med <- loc+b*Delta
  d <- (x - med)/sqrt(sigma2)
  dens <- 2 * dt(d, df = nu) * pt(sqrt((1+nu)/(d^2 + nu)) * d * shape, 1 + nu)/sqrt(sigma2)
  return(dens)
}

dSNCC <- function(y, mu, sigma2, shape, nu){
  k1 <- nu[1] / nu[2]^(1/2) + 1 - nu[1]
  b <-  -sqrt(2/pi) * k1
  delta <- shape / (sqrt(1 + shape^2))
  Delta <- sqrt(sigma2)*delta
  med<- mu + b * Delta
  dens <- 2 * (nu[1] * dnorm(y, med, sqrt(sigma2/nu[2])) *
                 pnorm(sqrt(nu[2]) * shape * sigma2^(-1/2) * (y - med)) + 
                 (1 - nu[1]) * dnorm(y, med, sqrt(sigma2)) * pnorm(shape * sigma2^(-1/2) * (y - med)))
  return(dens)
}

dSSC <- function(y, mu, sigma2, shape,nu){
  k1<- 2 * nu/(2 * nu-1)
  b<-  -sqrt(2/pi) * k1
  delta <- shape / (sqrt(1 + shape^2))
  Delta <- sqrt(sigma2) * delta
  med<- mu + b * Delta
  resp <- vector(mode = "numeric", length = length(y))
  for (i in 1:length(y)) {
    f <- function(u) 2 * nu * u^(nu - 1) * dnorm(y[i], med[i], sqrt(sigma2/u)) *
      pnorm(u^(1/2) * shape * (sigma2^(-1/2)) * (y[i] - med[i]))
    resp[i] <- integrate(f, 0, 0.9999)$value
  }
  return(resp)
}

d.mixedSTC <- function(y,x, beta, sigma2, shape, nu,r,alpha,class){
  g <- length(sigma2)
  dens = 0
  for(j in 1:g){
    Cent=y[class==j]-x[class==j,]%*%beta[,j]
    rnew=r[class==j,]
    for(i in 1:length(Cent)){
      if(j!=g) pI=c(exp(rnew[i,]%*%alpha[,j])/(1+sum(c(exp(t(alpha)%*%rnew[i,])))))
      if(j==g) pI=c(1/(1+sum(c(exp(t(alpha)%*%rnew[i,])))))
      dens <- dens + pI * dt.lsC( Cent[i],0, sigma2[j], shape[j], nu[j])
    }
  }
  return(dens)
}

d.mixedSNCC <- function(y,x, beta, sigma2, shape, nu,r,alpha,class){
  g <- length(sigma2)
  dens = 0
  for(j in 1:g){
    Cent=y[class==j]-x[class==j,]%*%beta[,j]
    rnew=r[class==j,]
    for(i in 1:length(Cent)){
      if(j!=g) pI=c(exp(rnew[i,]%*%alpha[,j])/(1+sum(c(exp(t(alpha)%*%rnew[i,])))))
      if(j==g) pI=c(1/(1+sum(c(exp(t(alpha)%*%rnew[i,])))))
      dens <- dens + pI * dSNCC(Cent[i],0, sigma2[j], shape[j], nu[(2*j-1) : (2*j)])
    }
  }
  return(dens)
}

d.mixedSSC <- function(y,x, beta, sigma2, shape, nu,r,alpha,class){
  g <- length(sigma2)
  dens = 0
  for(j in 1:g){
    Cent=y[class==j]-x[class==j,]%*%beta[,j]
    rnew=r[class==j,]
    for(i in 1:length(Cent)){
      if(j!=g) pI=c(exp(rnew[i,]%*%alpha[,j])/(1+sum(c(exp(t(alpha)%*%rnew[i,])))))
      if(j==g) pI=c(1/(1+sum(c(exp(t(alpha)%*%rnew[i,])))))
      dens <- dens + pI * dSSC(Cent[i],0, sigma2[j], shape[j], nu[j])
    }
  }
  return(dens)
}


smsn.mixReg.numesma <- function(y, 
                                x, 
                                nu, 
                                beta = NULL,  
                                sigma2 = NULL, 
                                shape = NULL, 
                                r,
                                alpha=NULL,
                                Class = NULL, 
                                g = NULL,
                                p=NULL,
                                family = "Skew.normal", 
                                error = 0.00001, 
                                iter.max = 100){
  
  
  if(is.null(sigma2) || is.null(alpha) || is.null(beta) || is.null(nu)||is.null(shape)){
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
    sigma2=rep(0,0)
    for(j in 1:g){
      beta[,j]=solve(t(x[class==j,])%*%x[class==j,])%*%t(x[class==j,])%*%y[class==j]
      sigma2[j]=var(y[class==j]-x[class==j,]%*%beta[,j])
    }
    pI = table(class)/n
    CC=1/pI[g]
    rstar=apply(r,2,mean)
    q=length(rstar)
    al<-function(par){
      ll=par[1]
      for(i in 2:(q-1)){
        ll=ll+par[i]*rstar[i]
      }
      ll=CC
    }
    opt<- nlminb(rep(1,q), al, lower = rep(0.01, q), 
                 upper = rep(20, q), control = list(iter.max = 2000))
    alpha=opt$par
    alpha=matrix(rep(alpha,g-1),ncol=(g-1))
    
  }
  g = length(pI)
  ################################################################################
  ##                                   Skew-t                                   ##
  ################################################################################
  if (family == "Skew.t"){
    shape=rep(2,g)
    nu=rep(10,g)
    k1 <- sqrt(nu/2) * gamma((nu - 1)/2) / gamma(nu/2)
    b <- -sqrt(2/pi) * k1
    
    n <- length(y)
    p <- ncol(x)
    
    delta <- Delta  <- Gama <- Gamaaux1 <- Gamaaux2 <- rep(0,g)
    
    mu <- matrix(0,n,g)
    media <- matrix(0,n,g)
    
    for (k in 1:g){
      delta[k] <- shape[k] / (sqrt(1 + shape[k]^2))
      Delta[k] <- sqrt(sigma2[k]) * delta[k]
      Gama[k] <- sigma2[k] - Delta[k]^2
      media[,k] <- x %*% beta[,k]
      mu[,k] <- media[, k] + b[k] * Delta[k]
    }
    
    teta <- c(as.vector(beta),as.vector(alpha), Delta, Gama,  nu)
    
    beta.old <- beta
    Delta.old <- Delta
    Gama.old <- Gama
    
    criterio <- 1
    count <- 0
    
    lk <- sum(log(d.mixedSTC(y,x, beta, sigma2, shape, nu,r,alpha,class)))
    
    while((criterio > error) && (count <= iter.max)){
      count <- count + 1
      tal <- matrix(0, n, g)
      for(j in 1:g){
        for(i in which(class==j)){
          e=x[i,]%*%beta[,j]
          pI=c(exp(t(alpha)%*%r[i,])/(1+sum(c(exp(t(alpha)%*%r[i,])))))
          pI=c(pI,c(1/(1+sum(c(exp(t(alpha)%*%r[i,]))))))
          tal[i, j] =pI[j] * dt.lsC(y[i]-e,0, sigma2[j], shape[j], nu[j])
        }
      }
      for(k in 1:n) if(all(tal[k,] == 0)) tal[k,] = .Machine$double.xmin
      tal = tal/rowSums(tal)
      class=apply(tal, 1, which.max)
      S1 <- matrix(0, n, g)
      S2 <- matrix(0, n, g)
      S3 <- matrix(0, n, g)
      for (j in 1:g){
        dj <- ((y[class==j] - mu[class==j,j]) / sqrt(sigma2[j]))^2
        Mtij2 <- 1/(1 + (Delta[j]^2) * (Gama[j]^(-1)))
        Mtij <- sqrt(Mtij2)
        mutij <- Mtij2 * Delta[j] * (Gama[j]^(-1))*(y[class==j] - mu[class==j,j]) + b[j]
        A <- (mutij - b[j]) / Mtij
        
        E = (2 * (nu[j])^(nu[j]/2) * gamma((2 + nu[j])/2) *
               ((dj + nu[j] + A^2))^(-(2 + nu[j])/2)) / (gamma(nu[j]/2) * 
                                                           pi * sqrt(sigma2[j]) * 
                                                           dt.lsC(y[class==j], media[class==j,j], sigma2[j], shape[j], nu[j]))
        u = ((4*(nu[j])^(nu[j]/2) * gamma((3 + nu[j])/2) * 
                (dj + nu[j])^(-(nu[j] + 3)/2)) / (gamma(nu[j]/2) * sqrt(pi) * 
                                                    sqrt(sigma2[j]) * 
                                                    dt.lsC(y[class==j], media[class==j,j], sigma2[j], shape[j], nu[j])) ) *
          pt(sqrt((3 + nu[j])/(dj + nu[j])) * A, 3 + nu[j])
        
        
        
        S1[class==j,j] <- tal[class==j,j] * u
        S2[class==j,j] <- tal[class==j,j] * (mutij * u + Mtij*E)
        S3[class==j,j] <- tal[class==j,j] * (mutij^2  *u + Mtij2 + Mtij * (mutij + b[j]) * E)
        
        S1aux<-as.vector(S1[class==j,j])
        S2aux<-as.vector(S2[class==j,j]/S1[class==j,j])
        
        ### M-step: 
        
        Delta[j] <- sum(S2[class==j, j] * (y[class==j] - mu[class==j,j] + b[j] * Delta[j])) / sum(S3[class==j,j])
        Gamaaux1[j] <- sum(S1[class==j,j]*(y[class==j] - mu[class==j,j] + b[j] * Delta[j])^2 - 2*(y[class==j] - mu[class==j,j] + b[j] * Delta[j])*
                             Delta[j] * S2[class==j,j] + Delta[j]^2 * S3[class==j,j]) / sum(tal[class==j,j])
        Gamaaux2[j]<-Gamaaux1[j] * sum(tal[,j])
        beta[,j] <- solve(t(x[class==j,]) %*% diag(S1aux) %*% x[class==j,]) %*% t(x[class==j,]) %*% diag(S1aux) %*% (y[class==j] - Delta[j] * S2aux)
        media[class==j,j] <- x[class==j,] %*% beta[,j]
        mu[class==j] <- media[class==j,j] + b[j]*Delta[j]
        
        ss=0
        for(i in 1:n) ss=ss+(tal[i, j]-c(exp(t(alpha)%*%r[i,])/(1+sum(c(exp(t(alpha)%*%r[i,])))))[j])*r[i,]
        if(j!=g) alpha[,j]=((4/sum(diag(r%*%t(r))))*ss)+alpha[,j]
      }
      
      Gamafim <- (1/n) * sum(Gamaaux2)
      
      for (jj in 1:g){
        Gama[jj]<- Gamafim
        sigma2[jj] <- Gama[jj] + Delta[jj]^2
        shape[jj] <- ((sigma2[jj]^(-1/2)) * Delta[jj] )/(sqrt(1 - (Delta[jj]^2) * (sigma2[jj]^(-1))))
      }
      
      logvero.ST <- function(nu) -sum(log( d.mixedSTC(y,x, beta, sigma2, shape, nu,r,alpha,class)))
      nu <- optim(nu, logvero.ST, method = "L-BFGS-B", lower = 2, upper = 50)$par
      k1 <- sqrt(nu/2) * gamma((nu - 1)/2) / gamma(nu/2)
      b <- -sqrt(2/pi) * k1
      lk1 <- sum(log( d.mixedSTC(y,x, beta, sigma2, shape, nu,r,alpha,class) ))
      
      
      param <- teta
      teta <- c(as.vector(beta),as.vector(alpha), Delta, Gama, nu)
      criterio <- abs(lk1/lk-1)
      
      beta.old <- beta
      Delta.old <- Delta
      Gama.old <- Gama
      lk <- lk1
    } 
    
    icl <- 0
    for(j in 1:g){
      Cent=y[class==j]-x[class==j,]%*%beta[,j]
      rnew=r[class==j,]
      for(i in 1:length(Cent)){
        if(j!=g) pI=c(exp(rnew[i,]%*%alpha[,j])/(1+sum(c(exp(t(alpha)%*%rnew[i,])))))
        if(j==g) pI=c(1/(1+sum(c(exp(t(alpha)%*%rnew[i,])))))
        icl <- icl + sum(log(pI * dt.lsC(Cent[i], 0,sigma2[j], shape[j], nu[j])))
      }
    }
  }
  
  if (family == "Skew.cn"){
    nu=cbind(rep(10,g),rep(10,g))
    shape=rep(2,g)
    k1 <- nu[1, ]/ nu[2, ]^(1/2) + 1 - nu[1, ]
    b <- -sqrt(2/pi) * k1
    n <- length(y)
    p <- ncol(x)
    
    delta <- Delta <- Gama <- rep(0,g)
    
    media <- mu <- matrix(0, n, g)
    
    for (k in 1:g){
      delta[k] <- shape[k] / (sqrt(1 + shape[k]^2))
      Delta[k] <- sqrt(sigma2[k]) * delta[k]
      Gama[k] <- sigma2[k] - Delta[k]^2
      media[, k] <- x %*% beta[,k]
      mu[, k] <- media[,k] + b[k] * Delta[k]
    }
    
    teta <- c(as.vector(beta),as.vector(alpha), Delta, Gama, nu)
    beta.old<- beta
    Delta.old <- Delta
    Gama.old <- Gama
    nu.old<-nu
    
    criterio <- 1
    count <- 0
    
    lk <- sum(log( d.mixedSNCC(y, x, beta, sigma2, shape, as.vector(nu),r,alpha,class)))
    
    while((criterio > error) && (count <= iter.max)){
      count <- count + 1
      tal <- matrix(0, n, g)
      for(j in 1:g){
        for(i in which(class==j)){
          e=x[i,]%*%beta[,j]
          pI=c(exp(t(alpha)%*%r[i,])/(1+sum(c(exp(t(alpha)%*%r[i,])))))
          pI=c(pI,c(1/(1+sum(c(exp(t(alpha)%*%r[i,]))))))
          tal[i, j] =pI[j] * dt.lsC(y[i]-e, 0,sigma2[j], shape[j], nu[j])
        }
      }
      for(k in 1:n) if(all(tal[k,] == 0)) tal[k,] = .Machine$double.xmin
      tal = tal/rowSums(tal)
      class=apply(tal, 1, which.max)
      S1 <- matrix(0, n, g)
      S2 <- matrix(0, n, g)
      S3 <- matrix(0, n, g)
      for (j in 1:g){
        ### E-step: 
        dj <- ((y[class==j] - mu[class==j,j])/sqrt(sigma2[j]))^2
        Mtij2 <- 1/(1 + (Delta[j]^2) * (Gama[j]^(-1)))
        Mtij <- sqrt(Mtij2)
        mutij <- Mtij2 * Delta[j] * (Gama[j]^(-1)) * (y[class==j] - mu[class==j,j]) + b[j]
        A <- (mutij - b[j]) / Mtij
        
        u=(2/dSNCC(y[class==j], media[class==j,j], sigma2[j], shape[j], nu[, j]))*
          (nu[1, j] * nu[2, j] * dnorm(y[class==j], mu[class==j,j], sqrt(sigma2[j] / nu[2, j])) *
             pnorm(sqrt(nu[2,j]) * A) + (1 - nu[1,j]) * dnorm(y[class==j], mu[class==j,j], sqrt(sigma2[j])) * pnorm(A))
        E=(2/dSNCC(y[class==j], media[class==j,j], sigma2[j], shape[j], nu[, j])) *
          (nu[1] * sqrt(nu[2, j]) * dnorm(y[class==j], mu[class==j, j], sqrt(sigma2[j] / nu[2, j])) *
             dnorm(sqrt(nu[2, j]) * A)+(1 - nu[1, j]) * dnorm(y[class==j], mu[class==j,j], sqrt(sigma2[j])) * dnorm(A))
        
        
        S1[,j] <- tal[class==j,j] * u
        S2[,j] <- tal[class==j,j] * (mutij*u + Mtij * E)
        S3[,j] <- tal[class==j,j] * (mutij^2*u + Mtij2 + Mtij*(mutij + b[j]) * E)
        
        S1aux<-as.vector(S1[class==j,j])
        S2aux<-as.vector(S2[class==j,j]/S1[class==j,j])
        
        ### M-step:
        Delta[j] <- sum(S2[class==j,j] * (y[class==j] - mu[class==j,j] + b[j] * Delta[j])) / sum(S3[class==j,j])
        Gama[j] <- sum(S1[class==j,j] * (y[class==j] - mu[class==j,j] + b[j] * Delta[j])^2 - 2*(y[class==j] - mu[class==j,j] + b[j]*Delta[j]) * 
                         Delta[j] * S2[,j] + Delta[j]^2 * S3[,j]) / sum(tal[,j])
        sigma2[j] <- Gama[j] + Delta[j]^2
        shape[j] <- ((sigma2[j]^(-1/2)) * Delta[j] )/(sqrt(1 - (Delta[j]^2) * (sigma2[j]^(-1))))
        beta[,j] <- solve(t(x[class==j,]) %*% diag(S1aux) %*% x[class==j,]) %*% t(x[class==j,]) %*% diag(S1aux) %*% (y[class==j] - Delta[j] * S2aux)
        media[class==j,j]<- x[class==j,] %*% beta[,j]
        mu[class==j,j]<- media[class==j,j] + b[j] * Delta[j]
      }
      
      logvero.SNC <- function(nu) -sum(log( d.mixedSNCC(y, x, beta, sigma2, shape, as.vector(nu),r,alpha,class) ))
      NUOPT = try(nlminb(as.vector(nu), logvero.SNC, lower = rep(0.01, 2), upper = rep(0.49, 0.99)))
      if (!('try-error' %in% class(NUOPT))) {
        nu = matrix(NUOPT$par, 2, g) 
      } else {nu = nu}
      lk1 <- sum(log(d.mixedSNCC(y, x, beta, sigma2, shape, as.vector(nu),r,alpha,class)))
      
      
      param <- teta
      teta <- c(as.vector(beta),as.vector(alpha), Delta, Gama, nu)
      
      criterio <- abs(lk1/lk-1)
      
      beta.old <- beta
      Delta.old <- Delta
      Gama.old <- Gama
      nu.old <- nu
      lk<-lk1
      
    }
    
    icl <- 0
    for(j in 1:g){
      Cent=y[class==j]-x[class==j,]%*%beta[,j]
      rnew=r[class==j,]
      for(i in 1:length(Cent)){
        if(j!=g) pI=c(exp(rnew[i,]%*%alpha[,j])/(1+sum(c(exp(t(alpha)%*%rnew[i,])))))
        if(j==g) pI=c(1/(1+sum(c(exp(t(alpha)%*%rnew[i,])))))
        icl <- icl + sum(log(pI * dSNCC(Cent[i], 0,sigma2[j], shape[j], nu[,j])))
      }
    }
  }
  
  if (family == "Skew.slash"){
    nu=rep(10,g)
    shape=rep(2,g)
    k1<- 2 * nu/(2 * nu - 1)
    b<- -sqrt(2/pi) * k1
    n <- length(y)
    p<- ncol(x)
    
    delta <- Delta  <- Gama<- Gamaaux1 <- Gamaaux2<-rep(0,g)
    
    mu<-matrix(0,n,g)
    media<-matrix(0,n,g)
    
    for (k in 1:g){
      delta[k] <- shape[k] / (sqrt(1 + shape[k]^2))
      Delta[k] <- sqrt(sigma2[k]) * delta[k]
      Gama[k] <- sigma2[k] - Delta[k]^2
      media[, k] <- x %*% Abetas[, k]
      mu[, k] <- media[,k] + b[k] * Delta[k]
    }
    
    teta <- c(as.vector(beta),as.vector(alpha), Delta, Gama, nu)
    
    beta.old<- beta
    Delta.old <- Delta
    Gama.old <- Gama
    nu.old<-nu
    criterio <- 1
    count <- 0
    
    lk <- sum(log( d.mixedSSC(y,x, beta, sigma2, shape, nu,r,alpha,class) ))
    
    while((criterio > error) && (count <= iter.max)){
      count <- count + 1
      tal <- matrix(0, n, g)
      for(j in 1:g){
        for(i in which(class==j)){
          e=x[i,]%*%beta[,j]
          pI=c(exp(t(alpha)%*%r[i,])/(1+sum(c(exp(t(alpha)%*%r[i,])))))
          pI=c(pI,c(1/(1+sum(c(exp(t(alpha)%*%r[i,]))))))
          tal[i, j] =pI[j] * dSSC(y[i]-e,0, sigma2[j], shape[j], nu[j])
        }
      }
      for(k in 1:n) if(all(tal[k,] == 0)) tal[k,] = .Machine$double.xmin
      tal = tal/rowSums(tal)
      class=apply(tal, 1, which.max)
      S1 <- matrix(0, n, g)
      S2 <- matrix(0, n, g)
      S3 <- matrix(0, n, g)
      for (j in 1:g){
        dj <- ((y[class==j] - mu[class==j, j])/sqrt(sigma2[j]))^2
        Mtij2 <- 1/(1 + (Delta[j]^2) * (Gama[j]^(-1)))
        Mtij <- sqrt(Mtij2)
        mutij <- Mtij2 * Delta[j] * (Gama[j]^(-1)) * (y[class==j] - mu[class==j,j]) + b[j]
        A <- (mutij - b[j]) / Mtij
        u <- vector(mode="numeric", length = n)
        E <- vector(mode="numeric", length = n)
        for(i in which(class==j)){
          E[i] <- (((2^(nu[j] + 1))*nu[j]*gamma(nu[j] + 1))/(dSSC(y[i], media[i,j], sigma2[j], shape[j], nu[j]) *
                                                               pi * sqrt(sigma2[j]))) * 
            ((dj[i] + A[i]^2)^(-nu[j]-1)) * pgamma(1, nu[j] + 1,(dj[i] + A[i]^2)/2)
          faux <- function(u) u^(nu[j]+0.5) * exp(-u * dj[i]/2) * pnorm(u^(1/2) * A[i])
          aux22 <- integrate(faux,0 , 1)$value
          u[i] <- ((sqrt(2) * nu[j]) / (dSSC(y[i], media[i,j], sigma2[j], shape[j], nu[j]) * 
                                          sqrt(pi) * sqrt(sigma2[j]))) * aux22
        }
        
        S1[class==j,j] <- tal[class==j,j] * u
        S2[class==j,j] <- tal[class==j,j] * (mutij * u + Mtij * E)
        S3[class==j,j] <- tal[class==j,j] * (mutij^2 * u + Mtij2 + Mtij * (mutij + b[j]) * E)
        
        S1aux<-as.vector(S1[class==j,j])
        S2aux<-as.vector(S2[class==j,j]/S1[class==j,j])
        
        ### M-step: 
        
        Delta[j] <- sum(S2[class==j,j] * (y[class==j] - mu[class==j,j] + b[j] * Delta[j])) / sum(S3[class==j,j])
        Gamaaux1[j] <- sum(S1[class==j,j]*(y[class==j] - mu[,j] + b[j] * Delta[j])^2 - 
                             2 * (y[class==j] - mu[class==j,j] + b[j] * Delta[j]) * Delta[j] * S2[class==j,j] + Delta[j]^2 * S3[class==j,j]) / sum(tal[,j])
        Gamaaux2[j]<-Gamaaux1[j] * sum(tal[,j])
        Abetas[,j] <- solve(t(x[class==j,]) %*% diag(S1aux) %*% x[class==j,]) %*% t(x[class==j,]) %*% diag(S1aux) %*% (y[class==j] - Delta[j] * S2aux)
        media[class==j,j]<-x[class==j,] %*% beta[,j]
        mu[class==j,j]<- media[class==j,j] + b[j] * Delta[j]
      }
      
      Gamafim<-(1/n) * sum(Gamaaux2)
      
      for (jj in 1:g){
        Gama[jj]<- Gamafim
        sigma2[jj] <- Gama[jj] + Delta[jj]^2
        shape[jj] <- ((sigma2[jj]^(-1/2)) * Delta[jj] )/(sqrt(1 - (Delta[jj]^2) * (sigma2[jj]^(-1))))
      }
      
      logvero.SS <- function(nu) -sum(log( d.mixedSSC(y,x, beta, sigma2, shape, nu,r,alpha,class) ))
      nu <- optim(nu, logvero.SS, method = "L-BFGS-B", lower = 1, upper = 50)$par
      
      lk1 <- sum(log( d.mixedSSC(y,x, beta, sigma2, shape, nu,r,alpha,class) ))
      
      
      
      param <- teta
      teta <- c(as,vector(beta),as.vector(alpha), Delta, Gama,  nu)
      criterio<-abs(lk1/lk-1)
      
      beta.old <- beta
      Delta.old <- Delta
      Gama.old <- Gama
      nu.old<-nu
      lk<-lk1
    }
    
    icl <- 0
    for(j in 1:g){
      Cent=y[class==j]-x[class==j,]%*%beta[,j]
      rnew=r[class==j,]
      for(i in 1:length(Cent)){
        if(j!=g) pI=c(exp(rnew[i,]%*%alpha[,j])/(1+sum(c(exp(t(alpha)%*%rnew[i,])))))
        if(j==g) pI=c(1/(1+sum(c(exp(t(alpha)%*%rnew[i,])))))
        icl <- icl + sum(log(pI * dSSC(Cent[i], 0,sigma2[j], shape[j], nu[j])))
      }
    }
  }
  
  
  d <- g * (nrow(beta) + 3) +(g-1)*nrow(alpha)
  if(family == "Skew.cn") d = d + g
  aic <- -2*lk + 2*d
  bic <- -2*lk + log(n)*d
  edc <- -2*lk + 0.2*sqrt(n)*d
  icl <- -2*icl + log(n)*d
  obj.out <- list(beta=beta, sigma2 = sigma2, shape = shape, alpha = alpha, 
                  nu = nu, aic = aic, bic = bic, edc = edc, icl = icl,
                  iter = count, n = length(y), group = Class)
  
  class(obj.out) <- family
  obj.out
}

