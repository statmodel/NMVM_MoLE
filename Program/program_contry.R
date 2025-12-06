rm(list = ls())

WD.PATH = paste(getwd(),"/Functions", sep = "")

source(paste(WD.PATH, '/Additional.r', sep = ""))
source(paste(WD.PATH, '/GHST-mix-censored.r', sep = ""))
source(paste(WD.PATH, '/NIG-mix-censored.r', sep = ""))
source(paste(WD.PATH, '/NMVBS-mix-censored.r', sep = ""))
source(paste(WD.PATH, '/NMVL-mix-censored.r', sep = ""))
source(paste(WD.PATH, '/SL-mix-censored.r', sep = ""))
source(paste(WD.PATH, '/VG-mix-censored.r', sep = ""))
source(paste(WD.PATH, '/normal-mix-censored.r', sep = ""))

#=================Data===============

df=county_dataset <- read.csv("C:/Users/ASUS/Dropbox/Dr Setudeh/Paper 3/without CR/Naderi/GitHub/Real data/county_population_poverty.csv")
head(county_dataset)
summary(county_dataset)
colnames(county_dataset)

n=nrow(county_dataset)
n

#==============preprocessing===================

x=cbind(rep(1,n),county_dataset$log_population)
r=cbind(rep(1,n),county_dataset$log_population)
y=county_dataset$pct_pov_2021


#================QQplot==================
library(ggplot2)
library(gridExtra)
library(grid)
library(gtable)
library(dplyr)

dens <- density(y)
dens_df <- data.frame(
  x = dens$x,             
  y = dens$y
)

p1=ggplot() +
  geom_histogram(aes(x = y, y = ..density..), 
                 bins = 10, 
                 fill = rgb(0, 0.5, 0.5, 0.7), 
                 color = "black", 
                 closed = "left") +
  labs(x = "y", y = "Density", title = "") +
  theme_minimal()
p1
dens <- density(x[,2])
dens_df <- data.frame(
  x = dens$x,             
  y = dens$y
)
p2=ggplot() +
  geom_histogram(aes(x = x[,2], y = ..density..), 
                 bins = 10, 
                 fill = rgb(0, 0.5, 0.5, 0.7), 
                 color = "black", 
                 closed = "left") +
  labs(x = "x", y = "Density", title = "") +
  theme_minimal()
p2
Histogram_xy <- arrangeGrob(p1,p2, ncol = 2)

grid.newpage()
grid.draw(Histogram_xy)

cairo_pdf("Histogram_y.pdf", width =10, height = 5)
grid.draw(p1)
dev.off()
####################################
df <- data.frame(y = y)

library(ggplot2)
library(cowplot)

p1 <- ggplot(df, aes(x = y)) +
  geom_boxplot(fill = "orange", width = 0.3) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = "Boxplot", title = NULL)

p2 <- ggplot(df, aes(x = y)) +
  geom_histogram(aes(y = ..density..), 
                 bins = 15, 
                 fill = rgb(0, 0.5, 0.5, 0.6), 
                 color = "black") +
  theme_minimal() +
  labs(x = "Poverty", y = "Density")

combined_plot1 <- plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(1, 3))
print(combined_plot1)

df <- data.frame(x = x[,2])

p3 <- ggplot(df, aes(x = x)) +
  geom_boxplot(fill = "orange", width = 0.3) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = "Boxplot", title = NULL)

p4 <- ggplot(df, aes(x = x)) +
  geom_histogram(aes(y = ..density..), 
                 bins = 15, 
                 fill = rgb(0, 0.5, 0.5, 0.6), 
                 color = "black") +
  theme_minimal() +
  labs(x = "Logarithm of population", y = "Density")

combined_plot2 <- plot_grid(p3, p4, ncol = 1, align = "v", rel_heights = c(1, 3))
print(combined_plot2)

Histogram_xy <- arrangeGrob(combined_plot1,combined_plot2, ncol = 2)

grid.newpage()
grid.draw(Histogram_xy)

cairo_pdf("Histogram_y_county.pdf", width =10, height = 6)
grid.draw(Histogram_xy)
dev.off()
ggsave("Histogram_y_county.png", plot = Histogram_xy, width =10, height = 6, dpi = 300)

#==================splite data=================
df=county_dataset

set.seed(36) 
train_indices <- sample(seq_len(n), size = 0.7 * n)
train_data <- df[train_indices, ]
test_data  <- df[-train_indices, ]

n <- nrow(train_data)
x=as.matrix(cbind(cons=rep(1,n),train_data[, c("log_population")]))
r=as.matrix(cbind(cons=rep(1,n),train_data[, c("log_population")]))
y=train_data$pct_pov_2021


d1=d2=0
q=2
p=2


g = 2
Out.NMVBS.ST2 = mix.arch.ST.EM.censor (y, x,d1,d2,typecensor="none",beta=NULL, sigma2 = NULL,
                                       lambda = NULL, nu = NULL,  r, alpha=NULL,
                                       g = g, Class = NULL, p=p, error = 0.0001,
                                       iter.max = 1000, Stp.rule = "Log.like", 
                                       error.est = F, per = 1, print = F, fix.sigma = F)
Out.NMVBS.NIG2 = mix.arch.nig.EM.censor (y, x,d1,d2,typecensor="none",beta=NULL, sigma2 = NULL,
                                         lambda = NULL, nu = NULL,  r, alpha=NULL,
                                         g = g, Class = NULL, p=p, error = 0.0001,
                                         iter.max = 1000, Stp.rule = "Log.like", 
                                         error.est = F, per = 1, print = F, fix.sigma = F)
Out.NMVBS.SL2 = mix.arch.SLap.EM.censor (y, x,d1,d2,typecensor="none",beta=NULL, sigma2 = NULL,
                                         lambda = NULL,  r, alpha=NULL,
                                         g = g, Class = NULL, p=p, error = 0.0001,
                                         iter.max = 50, Stp.rule = "Log.like", 
                                         error.est = F, per = 1, print = F, fix.sigma = F)
Out.NMVBS.NMVBS2 = mix.arch.NMVBS.EM.censor (y, x,d1,d2,typecensor="none",beta=NULL, sigma2 = NULL,
                                             lambda = NULL,  r, alpha=NULL,alpha2=NULL,
                                             g = g, Class = NULL, p=p, error = 0.0001,
                                             iter.max = 1000, Stp.rule = "Log.like", 
                                             error.est = F, per = 1, print = F, fix.sigma = F)
Out.NMVBS.NMVL2 = mix.arch.NMVL.EM.censor(y, x,d1,d2,typecensor="none",beta=NULL, sigma2 = NULL,
                                          lambda = NULL,  r, alpha=NULL, alpha2=NULL,
                                          g = g, Class = NULL, p=p, error = 0.0001,
                                          iter.max = 1000, Stp.rule = "Log.like", 
                                          error.est = F, per = 1, print = F, fix.sigma = F)
Out.NMVBS.VG2 = mix.arch.VG.EM.censor(y, x,d1,d2,typecensor="none",beta=NULL, sigma2 = NULL,
                                      lambda = NULL, nu=NULL, r, alpha=NULL, 
                                      g = g, Class = NULL, p=p, error = 0.0001,
                                      iter.max = 1000, Stp.rule = "Log.like", 
                                      error.est = F, per = 1, print = F, fix.sigma = F)

Out.NMVBS.N2 = mix.arch.norm.EM.censor(y, x, d1,d2,typecensor="none", beta=NULL, sigma2 = NULL, 
                                       g = g, r, alpha=NULL, Class = NULL,
                                       error = 0.0001,p=p, iter.max = 100, 
                                       Stp.rule = "Log.like", 
                                       error.est = F, per = 1, print = F, fix.sigma = F)


rbind(c(Out.NMVBS.ST2$loglike,Out.NMVBS.ST2$aic,Out.NMVBS.ST2$bic,Out.NMVBS.ST2$edc,Out.NMVBS.ST2$aic_c,Out.NMVBS.ST2$abic,Out.NMVBS.ST2$iter,Out.NMVBS.ST2$crite),
      c(Out.NMVBS.NIG2$loglike,Out.NMVBS.NIG2$aic,Out.NMVBS.NIG2$bic,Out.NMVBS.NIG2$edc,Out.NMVBS.NIG2$aic_c,Out.NMVBS.NIG2$abic,Out.NMVBS.NIG2$iter,Out.NMVBS.NIG2$crite),
      c(Out.NMVBS.SL2$loglike,Out.NMVBS.SL2$aic,Out.NMVBS.SL2$bic,Out.NMVBS.SL2$edc,Out.NMVBS.SL2$aic_c,Out.NMVBS.SL2$abic,Out.NMVBS.SL2$iter,Out.NMVBS.SL2$crite),
      c(Out.NMVBS.NMVL2$loglike,Out.NMVBS.NMVL2$aic,Out.NMVBS.NMVL2$bic,Out.NMVBS.NMVL2$edc,Out.NMVBS.NMVL2$aic_c,Out.NMVBS.NMVL2$abic,Out.NMVBS.NMVL2$iter,Out.NMVBS.NMVL2$crite),
      c(Out.NMVBS.NMVBS2$loglike,Out.NMVBS.NMVBS2$aic,Out.NMVBS.NMVBS2$bic,Out.NMVBS.NMVBS2$edc,Out.NMVBS.NMVBS2$aic_c,Out.NMVBS.NMVBS2$abic,Out.NMVBS.NMVBS2$iter,Out.NMVBS.NMVBS2$crite),
      c(Out.NMVBS.VG2$loglike,Out.NMVBS.VG2$aic,Out.NMVBS.VG2$bic,Out.NMVBS.VG2$edc,Out.NMVBS.VG2$aic_c,Out.NMVBS.VG2$abic,Out.NMVBS.VG2$iter,Out.NMVBS.VG2$crite),
      c(Out.NMVBS.N2$loglike,Out.NMVBS.N2$aic,Out.NMVBS.N2$bic,Out.NMVBS.N2$edc,Out.NMVBS.N2$aic_c,Out.NMVBS.N2$abic,Out.NMVBS.N2$iter,Out.NMVBS.N2$crite)
)


set.seed(36)
pred_GHST2=predict.mix.reg.NMV(n,x,Out.NMVBS.ST2$beta, Out.NMVBS.ST2$sigma2, Out.NMVBS.ST2$lambda, rbind(Out.NMVBS.ST2$nu), Out.NMVBS.ST2$alpha,r,family = "GHST",SN = FALSE)
pred_NIG2=predict.mix.reg.NMV(n,x,Out.NMVBS.NIG2$beta, Out.NMVBS.NIG2$sigma2, Out.NMVBS.NIG2$lambda, rbind(Out.NMVBS.NIG2$nu), Out.NMVBS.NIG2$alpha,r,family ="NIG",SN = FALSE)
pred_SL2=predict.mix.reg.NMV(n,x,Out.NMVBS.SL2$beta, Out.NMVBS.SL2$sigma2, Out.NMVBS.SL2$lambda, rbind(Out.NMVBS.SL2$nu), Out.NMVBS.SL2$alpha,r,family ="SL",SN = FALSE)
pred_NMVBS2=predict.mix.reg.NMV(n,x,Out.NMVBS.NMVBS2$beta, Out.NMVBS.NMVBS2$sigma2, Out.NMVBS.NMVBS2$lambda, rbind(Out.NMVBS.NMVBS2$alpha2), Out.NMVBS.NMVBS2$alpha,r,family ="NMVBS",SN = FALSE)
pred_NMVL2=predict.mix.reg.NMV(n,x,Out.NMVBS.NMVL2$beta, Out.NMVBS.NMVL2$sigma2, Out.NMVBS.NMVL2$lambda, rbind(Out.NMVBS.NMVL2$alpha2), Out.NMVBS.NMVL2$alpha,r,family ="NMVL",SN = FALSE)
pred_VG2=predict.mix.reg.NMV(n,x,Out.NMVBS.VG2$beta, Out.NMVBS.VG2$sigma2, Out.NMVBS.VG2$lambda, rbind(Out.NMVBS.VG2$nu), Out.NMVBS.VG2$alpha,r,family ="VG",SN = FALSE)
pred_N2=predict.mix.reg.NMV(n,x,Out.NMVBS.N2$beta, Out.NMVBS.N2$sigma2, Out.NMVBS.N2$lambda, rbind(Out.NMVBS.N2$alpha2), Out.NMVBS.N2$alpha,r,family ="normal",SN = FALSE)


result2=rbind(c(sqrt(mean((pred_GHST2$y-y)^2)),sqrt(mean((pred_NIG2$y-y)^2)),sqrt(mean((pred_SL2$y-y)^2)),sqrt(mean((pred_NMVBS2$y-y)^2)),sqrt(mean((pred_NMVL2$y-y)^2)),sqrt(mean((pred_VG2$y-y)^2)),sqrt(mean((pred_N2$y-y)^2))))

result2

rbind(c((mean(abs(pred_GHST2$y-y)/y)),(mean(abs(pred_NIG2$y-y)/y)),(mean(abs(pred_SL2$y-y)/y)),(mean(abs(pred_NMVBS2$y-y)/y)),(mean(abs(pred_NMVL2$y-y)/y)),(mean(abs(pred_VG2$y-y)/y)),(mean(abs(pred_N2$y-y)/y))))

#================parameter================================
cbind(c(Out.NMVBS.ST2$beta, Out.NMVBS.ST2$sigma2, Out.NMVBS.ST2$lambda, rbind(Out.NMVBS.ST2$nu),c(0,0), Out.NMVBS.ST2$alpha),
      c(Out.NMVBS.NIG2$beta, Out.NMVBS.NIG2$sigma2, Out.NMVBS.NIG2$lambda, rbind(Out.NMVBS.NIG2$nu),c(0,0), Out.NMVBS.NIG2$alpha),
      c(Out.NMVBS.SL2$beta, Out.NMVBS.SL2$sigma2, Out.NMVBS.SL2$lambda, c(0,0,0,0), Out.NMVBS.SL2$alpha),
      c(Out.NMVBS.NMVBS2$beta, Out.NMVBS.NMVBS2$sigma2, Out.NMVBS.NMVBS2$lambda, rbind(Out.NMVBS.NMVBS2$alpha2),c(0,0), Out.NMVBS.NMVBS2$alpha),
      c(Out.NMVBS.NMVL2$beta, Out.NMVBS.NMVL2$sigma2, Out.NMVBS.NMVL2$lambda, rbind(Out.NMVBS.NMVL2$alpha2),c(0,0), Out.NMVBS.NMVL2$alpha),
      c(Out.NMVBS.VG2$beta, Out.NMVBS.VG2$sigma2, Out.NMVBS.VG2$lambda, rbind(Out.NMVBS.VG2$nu), Out.NMVBS.VG2$alpha),
      c(Out.NMVBS.N2$beta, Out.NMVBS.N2$sigma2, c(0,0), c(0,0,0,0), Out.NMVBS.N2$alpha))


ntest <- nrow(test_data)
xtest=as.matrix(cbind(cons=rep(1,ntest),test_data[, c("log_population")]))
rtest=as.matrix(cbind(cons=rep(1,ntest),test_data[, c("log_population")]))
ytest=test_data$pct_pov_2021

set.seed(36)
pred_GHST2test=predict.mix.reg.NMV(ntest,xtest,Out.NMVBS.ST2$beta, Out.NMVBS.ST2$sigma2, Out.NMVBS.ST2$lambda, rbind(Out.NMVBS.ST2$nu), Out.NMVBS.ST2$alpha,rtest,family = "GHST",SN = FALSE)
pred_NIG2test=predict.mix.reg.NMV(ntest,xtest,Out.NMVBS.NIG2$beta, Out.NMVBS.NIG2$sigma2, Out.NMVBS.NIG2$lambda, rbind(Out.NMVBS.NIG2$nu), Out.NMVBS.NIG2$alpha,rtest,family ="NIG",SN = FALSE)
pred_SL2test=predict.mix.reg.NMV(ntest,xtest,Out.NMVBS.SL2$beta, Out.NMVBS.SL2$sigma2, Out.NMVBS.SL2$lambda, rbind(Out.NMVBS.SL2$nu), Out.NMVBS.SL2$alpha,rtest,family ="SL",SN = FALSE)
pred_NMVBS2test=predict.mix.reg.NMV(ntest,xtest,Out.NMVBS.NMVBS2$beta, Out.NMVBS.NMVBS2$sigma2, Out.NMVBS.NMVBS2$lambda, rbind(Out.NMVBS.NMVBS2$alpha2), Out.NMVBS.NMVBS2$alpha,rtest,family ="NMVBS",SN = FALSE)
pred_NMVL2test=predict.mix.reg.NMV(ntest,xtest,Out.NMVBS.NMVL2$beta, Out.NMVBS.NMVL2$sigma2, Out.NMVBS.NMVL2$lambda, rbind(Out.NMVBS.NMVL2$alpha2), Out.NMVBS.NMVL2$alpha,rtest,family ="NMVL",SN = FALSE)
pred_VG2test=predict.mix.reg.NMV(ntest,xtest,Out.NMVBS.VG2$beta, Out.NMVBS.VG2$sigma2, Out.NMVBS.VG2$lambda, rbind(Out.NMVBS.VG2$nu), Out.NMVBS.VG2$alpha,rtest,family ="VG",SN = FALSE)
pred_N2test=predict.mix.reg.NMV(ntest,xtest,Out.NMVBS.N2$beta, Out.NMVBS.N2$sigma2, Out.NMVBS.N2$lambda, rbind(Out.NMVBS.N2$alpha2), Out.NMVBS.N2$alpha,rtest,family ="normal",SN = FALSE)

c(sqrt(mean((pred_GHST2test$y-ytest)^2)),sqrt(mean((pred_NIG2test$y-ytest)^2)),
  sqrt(mean((pred_SL2test$y-ytest)^2)),sqrt(mean((pred_NMVBS2test$y-ytest)^2)),
  sqrt(mean((pred_NMVL2test$y-ytest)^2)),sqrt(mean((pred_VG2test$y-ytest)^2)),
  sqrt(mean((pred_N2test$y-ytest)^2)))

rbind(c((mean(abs(pred_GHST2test$y-ytest)/ytest)),(mean(abs(pred_NIG2test$y-ytest)/ytest)),
        (mean(abs(pred_SL2test$y-ytest)/ytest)),(mean(abs(pred_NMVBS2test$y-ytest)/ytest)),
        (mean(abs(pred_NMVL2test$y-ytest)/ytest)),(mean(abs(pred_VG2test$y-ytest)/ytest)),
        (mean(abs(pred_N2test$y-ytest)/ytest))))



###############################
df_new=rbind(train_data,test_data)
df_new$index=c(rep(1,2193),rep(2,941))
df_new$index=factor(df_new$index,labels =c("Train","Test"))
df_new$class=as.factor(c(Out.NMVBS.NMVL2$group,pred_NMVL2test$class))
df_new[df_new$pct_pov_2021>20,6]=2
df_new$res=(df_new$pct_pov_2021-c(pred_NMVL2$yhat,pred_NMVL2test$yhat))/3
df_new$id <- c(train_indices,setdiff(1:3134,train_indices) )


colnames(df_new)

p1=ggplot(df_new, aes(x = log_population, y = pct_pov_2021, color = class)) +
  geom_point(alpha = 0.8) +
  geom_abline(intercept = 12.8548762, slope =-0.1236657 , color = "#8A0000") + 
  geom_abline(intercept = 26.2643890, slope = -0.5648443,  color = "#004D8A") +
  scale_color_manual(values = c( "#8A0000",  "#004D8A")) +
  labs(x = "Logarithm of population", y = "Poverty", title = "") +
  theme_minimal()+theme(legend.position = "top")



p2=ggplot(df_new, aes(x = id, y = res, color = index)) +
  geom_point(size = 1) +
  geom_hline(yintercept = 0, color = "#212121", linetype = "dashed") +
  scale_color_manual(values = c("#8A0000", "#004D8A")) +
  labs(x = "Index", y = "Residual", title = "") +
  theme_minimal() +
  theme(legend.position = "top")

grid2<- arrangeGrob(p1, p2,
                    ncol = 2, widths = c(1, 1))

grid.newpage()
grid.draw(grid2)

cairo_pdf("train_test_res_county.pdf", width = 10, height = 5)
grid.draw(grid2)
dev.off()
ggsave("train_test-res_county.png", plot = grid2, width =10, height = 5, dpi = 300)


#######################
result2=rbind(c(sqrt(mean(((pred_GHST2$y-y)/3)^2)),sqrt(mean(((pred_NIG2$y-y)/3)^2)),sqrt(mean(((pred_SL2$y-y)/3)^2)),sqrt(mean(((pred_NMVBS2$y-y)/3)^2)),sqrt(mean(((pred_NMVL2$y-y)/3)^2)),sqrt(mean(((pred_VG2$y-y)/3)^2)),sqrt(mean(((pred_N2$y-y)/3)^2))))
result2
rbind(c((mean(abs((pred_GHST2$y-y)/3)/y)),(mean(abs((pred_NIG2$y-y)/3)/y)),(mean(abs((pred_SL2$y-y)/3)/y)),(mean(abs((pred_NMVBS2$y-y)/3)/y)),(mean(abs((pred_NMVL2$y-y)/3)/y)),(mean(abs((pred_VG2$y-y)/3)/y)),(mean(abs((pred_N2$y-y)/3)/y))))



set.seed(36)
pred_GHST2test=predict.mix.reg.NMV(ntest,xtest,Out.NMVBS.ST2$beta, Out.NMVBS.ST2$sigma2, Out.NMVBS.ST2$lambda, rbind(Out.NMVBS.ST2$nu), Out.NMVBS.ST2$alpha,rtest,family = "GHST",SN = FALSE)
pred_NIG2test=predict.mix.reg.NMV(ntest,xtest,Out.NMVBS.NIG2$beta, Out.NMVBS.NIG2$sigma2, Out.NMVBS.NIG2$lambda, rbind(Out.NMVBS.NIG2$nu), Out.NMVBS.NIG2$alpha,rtest,family ="NIG",SN = FALSE)
pred_SL2test=predict.mix.reg.NMV(ntest,xtest,Out.NMVBS.SL2$beta, Out.NMVBS.SL2$sigma2, Out.NMVBS.SL2$lambda, rbind(Out.NMVBS.SL2$nu), Out.NMVBS.SL2$alpha,rtest,family ="SL",SN = FALSE)
pred_NMVBS2test=predict.mix.reg.NMV(ntest,xtest,Out.NMVBS.NMVBS2$beta, Out.NMVBS.NMVBS2$sigma2, Out.NMVBS.NMVBS2$lambda, rbind(Out.NMVBS.NMVBS2$alpha2), Out.NMVBS.NMVBS2$alpha,rtest,family ="NMVBS",SN = FALSE)
pred_NMVL2test=predict.mix.reg.NMV(ntest,xtest,Out.NMVBS.NMVL2$beta, Out.NMVBS.NMVL2$sigma2, Out.NMVBS.NMVL2$lambda, rbind(Out.NMVBS.NMVL2$alpha2), Out.NMVBS.NMVL2$alpha,rtest,family ="NMVL",SN = FALSE)
pred_VG2test=predict.mix.reg.NMV(ntest,xtest,Out.NMVBS.VG2$beta, Out.NMVBS.VG2$sigma2, Out.NMVBS.VG2$lambda, rbind(Out.NMVBS.VG2$nu), Out.NMVBS.VG2$alpha,rtest,family ="VG",SN = FALSE)
pred_N2test=predict.mix.reg.NMV(ntest,xtest,Out.NMVBS.N2$beta, Out.NMVBS.N2$sigma2, Out.NMVBS.N2$lambda, rbind(Out.NMVBS.N2$alpha2), Out.NMVBS.N2$alpha,rtest,family ="normal",SN = FALSE)

c(sqrt(mean(((pred_GHST2test$y-ytest)/3)^2)),sqrt(mean(((pred_NIG2test$y-ytest)/3)^2)),
  sqrt(mean(((pred_SL2test$y-ytest)/3)^2)),sqrt(mean(((pred_NMVBS2test$y-ytest)/3)^2)),
  sqrt(mean(((pred_NMVL2test$y-ytest)/3)^2)),sqrt(mean(((pred_VG2test$y-ytest)/3)^2)),
  sqrt(mean(((pred_N2test$y-ytest)/3)^2)))

rbind(c((mean(abs((pred_GHST2test$y-ytest)/3)/ytest)),(mean(abs((pred_NIG2test$y-ytest)/3)/ytest)),
        (mean(abs((pred_SL2test$y-ytest)/3)/ytest)),(mean(abs((pred_NMVBS2test$y-ytest)/3)/ytest)),
        (mean(abs((pred_NMVL2test$y-ytest)/3)/ytest)),(mean(abs((pred_VG2test$y-ytest)/3)/ytest)),
        (mean(abs((pred_N2test$y-ytest)/3)/ytest))))

#############
df_new
sum(df_new$county=="hawaii")
df_new1=df_new[-which(df_new$county=="hawaii")]
library(ggplot2)
library(sf)
library(tigris)
library(dplyr)
library(gridExtra)
library(grid)
library(gtable)

options(tigris_use_cache = TRUE)
counties_map <- counties(cb = TRUE, resolution = "20m", year = 2022, class = "sf")

df_new$county <- iconv(df_new$county, from = "", to = "UTF-8", sub = "byte")
df_new$county <- gsub(" County", "", df_new$county)
df_new$county <- tolower(df_new$county)
counties_map$NAME <- tolower(counties_map$NAME)

map_data <- counties_map %>%
  left_join(df_new, by = c("NAME" = "county"))

map_data <- map_data %>% filter(!is.na(class))

p <- ggplot(map_data) +
  geom_sf(aes(fill = as.factor(class)), color = "white", size = 0.1) +
  scale_fill_manual(
    values = c("1" = "#8A0000", "2" = "#004D8A"),
    name = "Class",
    labels = c("1", "2")
  ) +
  labs(title = "County Clustering Map") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  ) +
  coord_sf(
    xlim = range(st_coordinates(map_data)[,1], na.rm = TRUE),
    ylim = range(st_coordinates(map_data)[,2], na.rm = TRUE),
    expand = FALSE
  )

cairo_pdf("county_map.pdf", width = 12, height = 6)
grid.draw(p)
dev.off()
ggsave("county_map.png", plot = p, width =12, height = 6, dpi = 300)


write.csv(df_new, "my_data_county.csv", row.names = FALSE)
