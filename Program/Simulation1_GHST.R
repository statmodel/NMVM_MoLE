rm(list=ls())
start.time <- Sys.time()
ptm<-proc.time()
#######------------import function  -------------#######
WD.PATH = paste(getwd(),"/Functions", sep = "")
source(paste(WD.PATH, '/Additional.r', sep = ""))
source(paste(WD.PATH, '/GHST-mix.r', sep = ""))
#######------------parameter -------------#######
nrun = 500
nsize = c(50,500, 1000, 2000)
g = 2 
q=3
p=4
beta0=matrix(c(0.5,-1.5,-1,-3,-1,2,1.5,3),ncol=g)
sigma20 = c(0.5,1.5)
lambda0 = c(3,-2)
alpha0 = matrix(c(0,2,4),ncol=1)
nu0 = cbind(2, 4)
for (m in 1:length(nsize)) {
  n = nsize[m]
  result1=matrix(NA,nrow=(3+p)*g+q*(g-1),ncol=nrun)
  result2=matrix(NA,nrow=nsize[m],ncol=nrun*2)
  result3=matrix(NA,nrow=8,ncol=nrun)
  mis = 0
  l = 0
  i=0
  while(l <= (nrun - 1)){
    l = l + 1
    x=cbind(rep(1,n),runif(n,-1,4),runif(n,-3,3),runif(n,1,4))
    r=cbind(rep(1,n),runif(n,-2,1),runif(n,-1,1))
    out.N = r.mix.reg.NMV(n,x,beta0, sigma20, lambda0, nu0, alpha0,r, family = "GHST")
    Class = out.N$class
    y = out.N$y
    Out.NMVBS = mix.reg.ST.EM (y, x,beta=NULL, sigma2 = NULL,
                               lambda = NULL, nu = NULL,  r, alpha=NULL,
                               g = g, Class = NULL, p=p, error = 0.00001,
                               iter.max = 1000, Stp.rule = "Log.like", 
                               error.est = F, per = 1, print = F, fix.sigma = F)
    result1[,l] = t(as.matrix(c(as.vector(Out.NMVBS$sigma2), as.vector(Out.NMVBS$lambda), as.vector(Out.NMVBS$alpha),
                                as.vector(Out.NMVBS$beta), as.vector(Out.NMVBS$nu))))
    result2[,c(2*i+1,2*i+2)]=cbind(Out.NMVBS$group,  Class)
    i=i+1
    result3[,l]=c(Out.NMVBS$loglike,Out.NMVBS$aic,Out.NMVBS$bic,Out.NMVBS$edc,Out.NMVBS$aic_c,Out.NMVBS$abic,Out.NMVBS$iter,Out.NMVBS$crite)
    print(l)
  }
  write.csv(result1, paste("Sim1/GHST_result1_",nsize[m],per[k],".csv"), append = TRUE,row.names = F,col.names = F)
  write.csv(result2, paste("Sim1/GHST_result2_",nsize[m],per[k],".csv"), append = TRUE,row.names = F,col.names = F)
  write.csv(result3, paste("Sim1/GHST_result3_",nsize[m],per[k],".csv"), append = TRUE,row.names = F,col.names = F)
  print(n)
}
proc.time()-ptm
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken