rm(list = ls())

WD.PATH = paste(getwd(),"/Functions", sep = "")

source(paste(WD.PATH, '/Additional.r', sep = ""))
source(paste(WD.PATH, '/GHST-mix.r', sep = ""))
source(paste(WD.PATH, '/NIG-mix-censored.r', sep = ""))
source(paste(WD.PATH, '/NMVBS-mix-censored.r', sep = ""))
source(paste(WD.PATH, '/NMVL-mix-censored.r', sep = ""))
source(paste(WD.PATH, '/SL-mix-censored.r', sep = ""))
source(paste(WD.PATH, '/VG-mix-censored.r', sep = ""))
source(paste(WD.PATH, '/normal-mix-censored.r', sep = ""))

#=================Data===============

df=tone_dataset <- read.csv("C:/Users/ASUS/Dropbox/Dr Setudeh/Paper 3/without CR/Naderi/GitHub/Real data/tone_dataset.csv")
head(tone_dataset)
summary(tone_dataset)
colnames(tone_dataset)

n=nrow(tone_dataset)
n

#==============preprocessing===================

x=cbind(rep(1,n),tone_dataset$stretchratio)
r=cbind(rep(1,n),tone_dataset$stretchratio)
y=tone_dataset$tuned


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

p2 <- ggplot(aes(y = y)) +
  geom_boxplot(fill = 'orange') +
  coord_flip() +
  labs(title = "Boxplot", y = "y")

p2

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
  labs(x = "Perceived tone ratio", y = "Density")

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
  labs(x = "Actual tone ratio", y = "Density")

combined_plot2 <- plot_grid(p3, p4, ncol = 1, align = "v", rel_heights = c(1, 3))
print(combined_plot2)

Histogram_xy <- arrangeGrob(combined_plot1,combined_plot2, ncol = 2)

grid.newpage()
grid.draw(Histogram_xy)

cairo_pdf("Histogram_y_tone.pdf", width =10, height = 6)
grid.draw(Histogram_xy)
dev.off()
ggsave("Histogram_y_tone.png", plot = Histogram_xy, width =10, height = 6, dpi = 300)

#==================splite data=================
df=tone_dataset
set.seed(123) 
train_indices <- sample(seq_len(n), size = 0.7 * n)
train_data <- df[train_indices, ]
test_data  <- df[-train_indices, ]

n <- nrow(train_data)
x=as.matrix(cbind(cons=rep(1,n),train_data[, c("stretchratio")]))
r=as.matrix(cbind(cons=rep(1,n),train_data[, c("stretchratio")]))
y=train_data$tuned


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
xtest=as.matrix(cbind(cons=rep(1,ntest),test_data[, c("stretchratio")]))
rtest=as.matrix(cbind(cons=rep(1,ntest),test_data[, c("stretchratio")]))
ytest=test_data$tuned

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



set.seed(1368)
pred_GHST2test=predict.mix.reg.NMV(ntest,xtest,matrix(par$GHST[1:10],ncol=2), par$GHST[11:12], par$GHST[13:14], rbind(par$GHST[15:16]), par$GHST[19:22],rtest,family = "GHST",SN = FALSE)
pred_NIG2test=predict.mix.reg.NMV(ntest,xtest,matrix(par$NIG[1:10],ncol=2), par$NIG[11:12], par$NIG[13:14], rbind(par$NIG[15:16]), par$NIG[19:22],rtest,family ="NIG",SN = FALSE)
pred_SL2test=predict.mix.reg.NMV(ntest,xtest,matrix(par$SL[1:10],ncol=2), par$SL[11:12], par$SL[13:14], rbind(par$SL[15:16]), par$SL[19:22],rtest,family ="SL",SN = FALSE)
pred_NMVBS2test=predict.mix.reg.NMV(ntest,xtest,matrix(par$NMVBS[1:10],ncol=2), par$NMVBS[11:12], par$NMVBS[13:14], rbind(par$NMVBS[15:16]), par$NMVBS[19:22],rtest,family ="NMVBS",SN = FALSE)
pred_NMVL2test=predict.mix.reg.NMV(ntest,xtest,matrix(par$NMVL[1:10],ncol=2), par$NMVL[11:12], par$NMVL[13:14], rbind(par$NMVL[15:16]), par$NMVL[19:22],rtest,family ="NMVL",SN = FALSE)
pred_VG2test=predict.mix.reg.NMV(ntest,xtest,matrix(par$VG[1:10],ncol=2), par$VG[11:12], par$VG[13:14],cbind(par$VG[15:16],par$VG[17:18]), par$VG[19:22],rtest,family ="VG",SN = FALSE)


c(sqrt(mean((pred_GHST2test$y-ytest)^2)),sqrt(mean((pred_NIG2test$y-ytest)^2)),
  sqrt(mean((pred_SL2test$y-ytest)^2)),sqrt(mean((pred_NMVBS2test$y-ytest)^2)),
  sqrt(mean((pred_NMVL2test$y-ytest)^2)),sqrt(mean((pred_VG2test$y-ytest)^2)))




#####################################
train_data$classSL=as.factor(Out.NMVBS.SL2$group)
colnames(train_data)=c("x","y", "Class")
p1=ggplot(train_data, aes(x = x, y = y, color = Class)) +
  geom_point(alpha = 0.8) +
  geom_abline(intercept = 1.93789326, slope = 0.02329392, color = "#8A0000") + 
  geom_abline(intercept = 0.03355212, slope = 0.99170124,  color = "#004D8A") +
  scale_color_manual(values = c( "#8A0000",  "#004D8A")) +
  labs(x = "Perceived tone ratio", y = "Actual tone ratio", title = "Train") +
  theme_minimal()+theme(legend.position = "none")





test_data$classSL=as.factor(pred_SL2test$class)
colnames(test_data)=c("x","y", "Class")
test_data$Class=as.factor(c(2,2,2,2,2,2,2,2,2,1,
                            1,1,2,1,2,2,2,1,1,1,
                            1,1,1,1,1,1,1,1,1,1,
                            1,2,1,1,1,1,1,1,1,1,
                            2,1,1,1,1))

p2=ggplot(test_data, aes(x = x, y = y, color = Class)) +
  geom_point(alpha = 0.8) +
  geom_abline(intercept = 1.93789326, slope = 0.02329392, color = "#8A0000") + 
  geom_abline(intercept = 0.03355212, slope = 0.99170124,  color = "#004D8A") +
  scale_color_manual(values = c( "#8A0000",  "#004D8A")) +
  labs(x = "x", y = "y", title = "Test") +
  theme_minimal()+theme(legend.position = "none")

grid2<- arrangeGrob(p1, p2,
                    ncol = 2, widths = c(1, 1))

grid.newpage()
grid.draw(grid2)

cairo_pdf("train_test.pdf", width = 8, height = 6)
grid.draw(grid2)
dev.off()
ggsave("train_test.png", plot = grid2, width =8, height = 6, dpi = 300)

#############################
library(ggplot2)
library(gridExtra)

df_new = rbind(train_data, test_data)
colnames(df_new) = c("x","y")
df_new$class = c(rep(1,105), rep(2,45))
df_new$class = factor(df_new$class, labels=c("Train","Test"))

# نمودار اول با رنگ + شکل متفاوت و گروه بزرگ با نقطه بزرگتر
p1 = ggplot(df_new, aes(x = x, y = y, 
                        color = class, 
                        shape = class)) +
  geom_point(size=2.5) +
  geom_abline(intercept = 1.93789326, slope = 0.02329392, color = "#8A0000") + 
  geom_abline(intercept = 0.03355212, slope = 0.99170124, color = "#004D8A") +
  scale_color_manual(values = c("#8A0000", "#004D8A")) +
  scale_shape_manual(values = c(16,17)) +         
  scale_size_manual(values = c(2.5, 2.5)) +       
  labs(x="Actual tone ratio", 
       y="Perceived tone ratio", 
       title="") +
  theme_minimal() +
  theme(legend.position = "top") +
  labs(color = NULL, shape = NULL, size = NULL)  # حذف نام ستون‌ها از legend

# باقیمانده‌ها
res = df_new$y - c(pred_SL2$yhat, pred_SL2test$yhat)
df_new$res = res

df_new$index <- c(train_indices, setdiff(1:150, train_indices))
p2 = ggplot(df_new, aes(x=index, y=res, color=class, shape=class)) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0, color="#212121", linetype="dashed") +
  scale_color_manual(values=c("#8A0000","#004D8A")) +
  scale_shape_manual(values=c(16,17)) +
  labs(x="Index", y="Residual", title="") +
  theme_minimal() +
  theme(legend.position="none")

grid2 <- arrangeGrob(p1, p2, ncol=2, widths=c(1,1))
grid.newpage()
grid.draw(grid2)

cairo_pdf("train_test_res.pdf", width = 10, height = 5)
grid.draw(grid2)
dev.off()
ggsave("train_test-res.png", plot = grid2, width =10, height = 5, dpi = 300)



# اضافه کردن چند نقطه خاص برای برچسب‌گذاری
data$label <- NA
data$label[c(56, 60, 85, 147)] <- c("56", "60", "85", "147")
# (نکته: داده واقعی شما ممکن است نیاز به برچسب‌گذاری سطرهای خاص داشته باشد)

# نمودار
ggplot(data, aes(x = actual_tone, y = perceived_tone)) +
  geom_point(alpha = 0.6, color = "darkgreen") +
  
  # خطوط رگرسیون برای هر گروه
  geom_smooth(data = subset(data, group == "N-MRM"),
              method = "lm", se = FALSE,
              linetype = "dotted", color = "blue", size = 1.2) +
  geom_smooth(data = subset(data, group == "NMVBS-MRM"),
              method = "lm", se = FALSE,
              linetype = "dashed", color = "red", size = 1.2) +
  
  # برچسب‌ها برای نقاط خاص (با شکل‌های مختلف)
  geom_point(data = data[c(56, 60, 85, 147), ], aes(shape = label, color = label), size = 4) +
  geom_text(data = data[c(56, 60, 85, 147), ],
            aes(label = label), hjust = -0.3, vjust = -0.5, size = 3.5) +
  
  scale_shape_manual(values = c("56" = 15, "60" = 17, "85" = 15, "147" = 15)) +
  scale_color_manual(values = c("56" = "orange", "60" = "red", "85" = "gold", "147" = "gold")) +
  
  # تنظیمات نمودار
  labs(title = "",
       x = "Actual tone",
       y = "Perceived tone") +
  theme_minimal() +
  theme(legend.position = "top") +
  
  # راهنمای دستی برای خطوط
  guides(color = guide_legend(override.aes = list(shape = NA))) +
  annotate("text", x = 2.7, y = 3.3, label = "NMVBS–MRM", color = "red", hjust = 0, size = 4) +
  annotate("text", x = 2.7, y = 2.9, label = "N–MRM", color = "blue", hjust = 0, size = 4)
