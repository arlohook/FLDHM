

library(tidyverse)
library(fda)
library(deSolve)
library(reshape)
library(ggpubr)
library(patchwork)
library(refund)
library(expm)
library(parallel)
library(foreach)
library(mgcv)
library(doParallel)
library(doRNG)  

# Read in data

Metrics = readRDS("FLDHM - Synthetic Gait.rds")

# make smoothing function
smooth.gcv = function(X){
  
  t = seq(0,1 , 0.01)
  basis = create.fourier.basis(c(0,1), 101)
  lambdas = c(1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5)
  gcv = numeric(length(lambdas))
  
  for(i in lambdas){
    
    par = fdPar(basis, 4, i)
    fit = smooth.basis(t, X, par, dfscale = 1.2)
    
    gcv[which(lambdas == i)] = mean(fit$gcv)
    
  }
  
  bestpar = fdPar(basis, 4, lambdas[which.min(gcv)])
  bestfd = smooth.basis(t, X, bestpar)$fd
  print(lambdas[which.min(gcv)[1]])
  return(bestfd)  
  
}

# create functional data
Metrics.fd = lapply(Metrics, function(X){smooth.gcv(X)})


# estimate derivatives
hip.fd = Metrics.fd$Hip
dhip.fd = deriv.fd(hip.fd, 1)
d2hip.fd = deriv.fd(hip.fd, 2)

knee.fd = Metrics.fd$Knee
dknee.fd = deriv.fd(knee.fd, 1)
d2knee.fd = deriv.fd(knee.fd, 2)


ankle.fd = Metrics.fd$Ankle
dankle.fd = deriv.fd(ankle.fd, 1)
d2ankle.fd = deriv.fd(ankle.fd, 2)

par(mfrow = c(1,3))
plot(hip.fd, main = "Hip")
plot(knee.fd, main = "Knee")
plot(ankle.fd, main = "Ankle")


# center data

hip.fd.c = center.fd(hip.fd)
dhip.fd.c = center.fd(dhip.fd)
d2hip.fd.c = center.fd(d2hip.fd)

knee.fd.c = center.fd(knee.fd)
dknee.fd.c = center.fd(dknee.fd)
d2knee.fd.c = center.fd(d2knee.fd)


ankle.fd.c = center.fd(ankle.fd)
dankle.fd.c = center.fd(dankle.fd)
d2ankle.fd.c = center.fd(d2ankle.fd)


# set up data frames for modelling
tfine = seq(0,1, 0.01)
n = ncol(Metrics$Hip)

set.seed(1)
Session = sample(1:5, n, T)

pffr.data.h = data.frame("Session" = Session,
                         "H" = I(-1*t(eval.fd(tfine, hip.fd.c))),
                         "K" = I(t(eval.fd(tfine, knee.fd.c))),
                         "dH" = I(-1*t(eval.fd(tfine, dhip.fd.c))),
                         "dK" = I(t(eval.fd(tfine, dknee.fd.c))),
                         "d2H" = I(t(eval.fd(tfine, d2hip.fd.c))),
                         "d2K" = I(t(eval.fd(tfine, d2knee.fd.c))))




pffr.data.k = data.frame("Session" = Session,
                         "H" = I(t(eval.fd(tfine, hip.fd.c))),
                         "K" = I(-1*t(eval.fd(tfine, knee.fd.c))),
                         "A" = I(t(eval.fd(tfine, ankle.fd.c))),
                         "dH" = I(t(eval.fd(tfine, dhip.fd.c))),
                         "dK" = I(-1*t(eval.fd(tfine, dknee.fd.c))),
                         "dA" = I(t(eval.fd(tfine, dankle.fd.c))),
                         "d2H" = I(t(eval.fd(tfine, d2hip.fd.c))),
                         "d2K" = I(t(eval.fd(tfine, d2knee.fd.c))),
                         "d2A" = I(t(eval.fd(tfine, d2ankle.fd.c))))


pffr.data.a = data.frame("Session" = Session,
                         "K" = I(t(eval.fd(tfine, knee.fd.c))),
                         "A" = I(-1*t(eval.fd(tfine, ankle.fd.c))),
                         "dK" = I(t(eval.fd(tfine, dknee.fd.c))),
                         "dA" = I(-1*t(eval.fd(tfine, dankle.fd.c))),
                         "d2K" = I(t(eval.fd(tfine, d2knee.fd.c))),
                         "d2A" = I(t(eval.fd(tfine, d2ankle.fd.c))))



# look at manual selection of smpar with fixed rich basis (n = 60)


lambdas <- 10^(-10:10)



# Folds by session/date
folds <- split(seq_len(n), pffr.data.h$Session)
kfold <- length(folds)

#set up function to get out of sample R-squared
rsq_for_lambda <- function(lambda) {
  # prediction containers
  h.pred <- matrix(NA_real_, n, 101)
  k.pred <- matrix(NA_real_, n, 101)
  a.pred <- matrix(NA_real_, n, 101)
  
  for (b in seq_len(kfold)) {
    idx_test <- folds[[b]]
    
    train.h <- pffr.data.h[-idx_test, ]
    train.k <- pffr.data.k[-idx_test, ]
    train.a <- pffr.data.a[-idx_test, ]
    
    test.h  <- pffr.data.h[idx_test, ]
    test.k  <- pffr.data.k[idx_test, ]
    test.a  <- pffr.data.a[idx_test, ]
    
    # NOTE: 'sp = lambda' is a starting value for the y-index margin.
    hip.pda <- pffr(
      d2H ~ H + dH + K + dK,
      data = train.h, yind = tfine,
      bs.yindex = list(k = 60, m = c(3, 2), bs = "cc", sp = lambda),
      bs.int = list(k = 60, m = c(3, 2), bs = "cc", sp = lambda),
      method = "REML"
    )
    
    knee.pda <- pffr(
      d2K ~ H + dH + K + dK + A + dA,
      data = train.k, yind = tfine,
      bs.yindex = list(k = 60, m = c(3, 2), bs = "cc", sp = lambda),
      bs.int = list(k = 60, m = c(3, 2), bs = "cc", sp = lambda),
      method = "REML"
    )
    
    ankle.pda <- pffr(
      d2A ~ K + dK + A + dA,
      data = train.a, yind = tfine,
      bs.yindex = list(k = 60, m = c(3, 2), bs = "cc", sp = lambda),
      bs.int = list(k = 60, m = c(3, 2), bs = "cc", sp = lambda),
      method = "REML"
    )
    
    h.pred[idx_test, ] <- predict(hip.pda,  test.h)
    k.pred[idx_test, ] <- predict(knee.pda, test.k)
    a.pred[idx_test, ] <- predict(ankle.pda, test.a)
  }
  
  # True matrices
  true.h <- pffr.data.h$d2H
  true.k <- pffr.data.k$d2K
  true.a <- pffr.data.a$d2A
  
  
  TSS <- function(M) {
    cm <- colMeans(M)
    sum(colSums((M - rep(cm, each = nrow(M)))^2))
  }
  RSS <- function(M, P) sum(colSums((M - P)^2))
  
  rsq.h <- 1 - RSS(true.h, h.pred) / TSS(true.h)
  rsq.k <- 1 - RSS(true.k, k.pred) / TSS(true.k)
  rsq.a <- 1 - RSS(true.a, a.pred) / TSS(true.a)

  
  c("hip" = rsq.h, "knee" = rsq.k, "ankle" = rsq.a)
  
  
}

# set up cluster
n_cores <- 20
cl <- makeCluster(n_cores)
registerDoParallel(cl)


clusterExport(
  cl,
  c("pffr.data.h", "pffr.data.k", "pffr.data.a", "tfine", "folds", "kfold", "n", "rsq_for_lambda"),
  envir = environment()
)

# find lambda with best out of sample r-squared
res_mat <- foreach(L = lambdas, .combine = rbind,
                   .packages = c("refund", "mgcv")) %dopar% {
                     out <- rsq_for_lambda(L)
                     c(lambda = L, out)
                   }





#saveRDS(res_mat, "Lambda CV BAM 101 Gait.rds")

#res_mat = readRDS("Lambda CV CC Gait.rds")

stopCluster(cl)

lam.h = res_mat[which.max(res_mat[,2]),1]
lam.k = res_mat[which.max(res_mat[,3]),1]
lam.a = res_mat[which.max(res_mat[,4]),1]
###############################################################################

hip.pda = pffr(d2H ~ H + dH + K + dK, data = pffr.data.h, yind = tfine, 
               bs.yindex = list(k = 60, m = c(3,2), bs = 'cc', sp = lam.h),
               bs.int = list(k = 60, m = c(3,2), bs = 'cc', sp = lam.h))
summary(hip.pda)


knee.pda = pffr(d2K ~ H + dH + K + dK + A + dA, data = pffr.data.k, yind = tfine, 
                bs.yindex = list(k = 60, m = c(3,2), bs = 'cc', sp = lam.k),
                bs.int = list(k = 60, m = c(3,2), bs = 'cc', sp = lam.k))

summary(knee.pda)

ankle.pda = pffr(d2A ~ K + dK + A + dA, data = pffr.data.a, yind = tfine, 
                 bs.yindex = list(k = 60, m = c(3,2), bs = 'cc', sp = lam.a),
                 bs.int = list(k = 60, m = c(3,2), bs = 'cc', sp = lam.a))
summary(ankle.pda)

# calculate integrated r-squared
R2.h = 1 - (sum(apply(matrix(hip.pda$residuals,n,101,T), 2, function(x){sum(x^2)})))/(sum(apply(matrix(hip.pda$y, n, 101, T), 2, function(x){sum((x-mean(x))^2)})))
R2.k = 1 - (sum(apply(matrix(knee.pda$residuals,n,101,T), 2, function(x){sum(x^2)})))/(sum(apply(matrix(knee.pda$y, n, 101, T), 2, function(x){sum((x-mean(x))^2)})))
R2.a = 1 - (sum(apply(matrix(ankle.pda$residuals,n,101,T), 2, function(x){sum(x^2)})))/(sum(apply(matrix(ankle.pda$y, n, 101, T), 2, function(x){sum((x-mean(x))^2)})))

# quick plot of functional coefficients
par(mfrow= c(2,3))
plot(hip.pda, scale = 0)

par(mfrow = c(3,3))
plot(knee.pda, scale = 0)

par(mfrow = c(2,3))
plot(ankle.pda, scale = 0)


# maek function to extract functional coefficients
extract.coefs = function(mod){
  
  
  coefs = coef(mod)
  smooth = lapply(coefs[[2]], function(x){x$coef})
  smooth[[1]]$value = smooth[[1]]$value+coefs$pterms[1,1]
  smooth[[1]]$se = smooth[[1]]$se + coefs$pterms[1,2]
  
  return(smooth)
}

# extract coefficients
hip.co = extract.coefs(hip.pda)
knee.co = extract.coefs(knee.pda)
ankle.co = extract.coefs(ankle.pda)


# make coeficients into functions
hfun = lapply(hip.co, function(X){splinefun(x = X$tfine.vec, y = X$value)})
kfun = lapply(knee.co, function(X){splinefun(x = X$tfine.vec, y = X$value)})
afun = lapply(ankle.co, function(X){splinefun(x = X$tfine.vec, y = X$value)})


# create state transition matrix
A_t = function(t){
  
  A = matrix(c(0,               1,                  0,                  0,                  0,                  0,   
               -hfun[[2]](t),  -hfun[[3]](t),      hfun[[4]](t),     hfun[[5]](t),             0,                  0,
               0,               0,                  0,                  1,                  0,                  0,
               kfun[[2]](t),  kfun[[3]](t),       -kfun[[4]](t),  -kfun[[5]](t),       kfun[[6]](t),      kfun[[7]](t),
               0,               0,                  0,                  0,                  0,                  1,
               0,               0,                afun[[2]](t),          afun[[3]](t),   -afun[[4]](t),      -afun[[5]](t)),
             ncol = 6, nrow = 6, byrow = T)
  
  return(A)
  
}

# make Bu matrix
B = matrix(c(0,1,0,1,0,1), nrow = 6)

Bu_t = function(t){
  
  U = matrix(
    c(0, hfun[[1]](t), 0, kfun[[1]](t), 0, afun[[1]](t)),
    ncol = 1)
  
  Bu = B*U
  return(Bu)  
}

# define the ode
LL = function(t,y,params){
  
  
  dy = A_t(t)%*%y + Bu_t(t)
  
  return(list(dy))
}

# quick test
initial = c(0,0,0,0,0,0)

solution <- ode(y = initial, times = tfine, func = LL, parms = NULL, method = "lsoda")

plot(solution)





# look at reconstruction of data

time = tfine

inits = cbind(t(eval.fd(0, hip.fd.c)),
              t(eval.fd(0, dhip.fd.c)), 
              t(eval.fd(0, knee.fd.c)),
              t(eval.fd(0, dknee.fd.c)),
              t(eval.fd(0, ankle.fd.c)),
              t(eval.fd(0, dankle.fd.c)))
colnames(inits) = c("H", "dH", "K", "dK", "A", "dA")
rownames(inits) = 1:n


h.mat = c()
dh.mat = c()
k.mat = c()
dk.mat = c()
a.mat = c()
da.mat = c()


for(i in 1:nrow(inits)){
  
  initial = inits[i,]
  
  solution <- ode(y = initial, times = time, func = LL, parms = NULL, method = "lsoda")
  
  h.mat = cbind(h.mat, solution[,2])
  dh.mat = cbind(dh.mat, solution[,3])
  k.mat = cbind(k.mat, solution[,4])
  dk.mat = cbind(dk.mat, solution[,5])
  a.mat = cbind(a.mat, solution[,6])
  da.mat = cbind(da.mat, solution[,7])

}


# add back in the mean

h.hat = apply(h.mat, 2, function(x){x+eval.fd(time, mean.fd(hip.fd))})
dh.hat = apply(dh.mat, 2, function(x){x+eval.fd(time, mean.fd(dhip.fd))})
k.hat = apply(k.mat, 2, function(x){x+eval.fd(time, mean.fd(knee.fd))})
dk.hat = apply(dk.mat, 2, function(x){x+eval.fd(time, mean.fd(dknee.fd))})
a.hat = apply(a.mat, 2, function(x){x+eval.fd(time, mean.fd(ankle.fd))})
da.hat = apply(da.mat, 2, function(x){x+eval.fd(time, mean.fd(dankle.fd))})

par(mfrow = c(1,3))
matplot(h.hat, type = 'l', main = "Esimated Hip", ylim = c(-20, 60))
matplot(eval.fd(time, hip.fd), type = 'l', main = "True Hip", ylim = c(-20, 60))
matplot(eval.fd(time, hip.fd)-h.hat, type = 'l', main = "Hip Residuals")

par(mfrow = c(1,3))
matplot(k.hat, type = 'l', main = "Esimated Knee", ylim = c(-20, 120))
matplot(eval.fd(time, knee.fd), type = 'l', main = "True Knee", ylim = c(-20, 120))
matplot(eval.fd(time, knee.fd) - k.hat, type = 'l', main = "Knee Residuals")

par(mfrow = c(1,3))
matplot(a.hat, type = 'l', main = "Esimated Ankle", ylim = c(-40,40))
matplot(eval.fd(time, ankle.fd), type = 'l', main = "True Ankle", ylim = c(-40,40))
matplot(eval.fd(time, ankle.fd)-a.hat, type = 'l', main = "Ankle Residuals")




# look at eigen values

eigmat = matrix(NA, nrow = 101, ncol = 6)
eigvecs = array(NA, c(101,6,6))

for(i in 1:101){
  
  A = A_t(time[i])
  
  EIG = eigen(A)
  
  eigmat[i,] = EIG$values
  
  eigvecs[i,,1] = t(EIG$vectors[,1])
  eigvecs[i,,2] = t(EIG$vectors[,2])
  eigvecs[i,,3] = t(EIG$vectors[,3])
  eigvecs[i,,4] = t(EIG$vectors[,4])
  eigvecs[i,,5] = t(EIG$vectors[,5])
  eigvecs[i,,6] = t(EIG$vectors[,6])
  
 
  
}

# take the real component for plotting
R.eigmat = apply(eigmat, 2, Re) 

real.plot = melt(R.eigmat)

TVES = ggplot(real.plot, aes(x = X1, y = value, colour = as.factor(X2)))+
  geom_line()+
  geom_hline(yintercept = 0, lty = 2)+
  labs(x = " ", y  ="Eigen Value \n (Real Part)")+
  theme_light()+
  facet_wrap( ~as.factor(X2),ncol = 6)+
  theme(legend.position = 'none', text = element_text(size = 14), axis.text.x = element_text(size = 10))
TVES

# look at eigen vectors

evec1 = Re(eigvecs[,,1])
colnames(evec1) =  c("H", "dH", "K", "dK", "A", "dA")
evc1.plot = melt(evec1)

evec2 = Re(eigvecs[,,2])
colnames(evec2) =  c("H", "dH", "K", "dK", "A", "dA")
evc2.plot = melt(evec2)

evec3 = Re(eigvecs[,,3])
colnames(evec3) =  c("H", "dH", "K", "dK", "A", "dA")
evc3.plot = melt(evec3)

evec4 = Re(eigvecs[,,4])
colnames(evec4) =  c("H", "dH", "K", "dK", "A", "dA")
evc4.plot = melt(evec4)

evec5 = Re(eigvecs[,,5])
colnames(evec5) =  c("H", "dH", "K", "dK", "A", "dA")
evc5.plot = melt(evec5)

evec6 = Re(eigvecs[,,6])
colnames(evec6) =  c("H", "dH", "K", "dK", "A", "dA")
evc6.plot = melt(evec6)


EVPLOT = data.frame(rbind(evc1.plot, evc2.plot, evc3.plot, evc4.plot, evc5.plot, evc6.plot))
EVPLOT$Eigenvector = as.factor(c(rep(1:6, each = 606)))


my_palette <- c(
  "#0072B2", # Blue
  "#D55E00", # Vermillion
  "#F0E442", # Yellow
  "#009E73", # Green
  "#CC79A7", # Reddish purple
  "#E69F00")  # Orange


TVCS = ggplot(EVPLOT, aes(x = X1, y = abs(value), colour = X2))+
  geom_line()+
  labs(x = "% Gait Cycle", y = "Eigen Vector Value \n (Real Part)", colour = "Variable")+
  theme_light()+
  facet_wrap(~Eigenvector, ncol = 6)+
  scale_colour_manual(values = my_palette)+
  theme(legend.position = 'bottom', text = element_text(size = 14), axis.text.x = element_text(size = 10))

TVCS
layout = "AAA
          AAA
          AAA
          AAA
          AAA
          AAA
          AAA
          BBB
          BBB
          BBB
          BBB
          BBB
          BBB
          BBB
          CCC"

EIGVALVEC = wrap_plots(TVES, TVCS, guide_area())+ plot_layout(design = layout)
EIGVALVEC
#ggsave("FLDHM Figure 5.3 Sensitivity.jpg", EIGVALVEC, dpi = 300, height = 16, width = 20, units = 'cm')




