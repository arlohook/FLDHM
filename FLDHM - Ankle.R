
library(tidyverse)
library(fda)
library(deSolve)
library(reshape)
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

# create smoothing function
smooth.gcv = function(X){
  
  t = seq(0, 1, 0.01)
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
ankle.fd = Metrics.fd$Ankle
dankle.fd = deriv.fd(ankle.fd, 1)
d2ankle.fd = deriv.fd(ankle.fd, 2)

par(mfrow = c(1,3))
plot(ankle.fd, main = "0")
plot(dankle.fd, main = "1")
plot(d2ankle.fd, main = "2")


# center data

ankle.fd.c = center.fd(ankle.fd)
dankle.fd.c = center.fd(dankle.fd)
d2ankle.fd.c = center.fd(d2ankle.fd)


# set up data for modelling
tfine = seq(0,1, 0.01)
n <- ncol(Metrics$Hip) 

set.seed(1)

pffr.data.a = data.frame("Session" = sample(1:5, n, T),     # sessions randomly allocated here as not preserved in data generation
                         "A" = I(-1*t(eval.fd(tfine, ankle.fd.c))),
                         "dA" = I(-1*t(eval.fd(tfine, dankle.fd.c))),
                         "d2A" = I(t(eval.fd(tfine, d2ankle.fd.c))))



# look at manual selection of smpar with fixed rich basis (n = 60)

lambdas <- 10^(-10:15)

        

# Folds by session
folds <- split(seq_len(n), pffr.data.a$Session)
kfold <- length(folds)


# set up k-fold validation function
rsq_for_lambda <- function(lambda) {
  # prediction containers
 
  a.pred <- matrix(NA, n, 101)
  
  for (b in seq_len(kfold)) {
    idx_test <- folds[[b]]
    
    
    train.a <- pffr.data.a[-idx_test, ]
    
    
    test.a  <- pffr.data.a[idx_test, ]

    ankle.pda <- pffr(
      d2A ~ A + dA,
      data = train.a, yind = tfine,
      bs.yindex = list(k = 60, m = c(3, 2), bs = "cc", sp = lambda),
      bs.int = list(k = 60, m = c(3, 2), bs = "cc", sp = lambda),
      method = "REML"
    )
    

    a.pred[idx_test, ] <- predict(ankle.pda, test.a)
  }
  
  # True matrices

  true.a <- pffr.data.a$d2A
  
  
  TSS <- function(M) {
    cm <- colMeans(M)
    sum(colSums((M - rep(cm, each = nrow(M)))^2))
  }
  RSS <- function(M, P) sum(colSums((M - P)^2))
  

  rsq.a <- 1 - RSS(true.a, a.pred) / TSS(true.a)
  
  c("ankle" = rsq.a)
  
  
}

# set up cluster
n_cores <- 20
cl <- makeCluster(n_cores)
registerDoParallel(cl)


clusterExport(
  cl,
  c("pffr.data.a", "tfine", "folds", "kfold", "n", "rsq_for_lambda"),
  envir = environment()
)


# calculate out of sample R-squared for each lambda
res_mat <- foreach(L = lambdas, .combine = rbind,
                   .packages = c("refund", "mgcv")) %dopar% {
                     out <- rsq_for_lambda(L)
                     c(lambda = L, out)
                   }




stopCluster(cl)

# take best lambda
lam = res_mat[which.max(res_mat[,2]),1]


# fit the model
ankle.pda = pffr(d2A ~ A + dA, data = pffr.data.a, yind = tfine, 
                 bs.yindex = list(k = 60, m = c(3,2), bs = 'cc', sp = lam),
                 bs.int = list(k = 60, m = c(3,2), bs = 'cc', sp = lam))

summary(ankle.pda)

# calculate integrated r-squared
R2 = 1 - (sum(apply(matrix(ankle.pda$residuals,n,101,T), 2, function(x){sum(x^2)})))/(sum(apply(matrix(ankle.pda$y, n, 101, T), 2, function(x){sum((x-mean(x))^2)})))


# plot parameters 
par(mfrow = c(1,3))
plot(ankle.pda, scale = 0)


# creaet function to extrcat functional coefs
extract.coefs = function(mod){
  
  
  coefs = coef(mod)
  smooth = lapply(coefs[[2]], function(x){x$coef})
  smooth[[1]]$value = smooth[[1]]$value+coefs$pterms[1,1]
  smooth[[1]]$se = smooth[[1]]$se + coefs$pterms[1,2]
  
  return(smooth)
}



ankle.co = extract.coefs(ankle.pda)



# set up bootstrap samples for confidence bands

B = 1000

set.seed(1)
boot.ind = matrix(sample(1:n, n*B, T), n, B)

# write boootstrap function
boot.fun = function(i){
  
  pda.data = pffr.data.a[boot.ind[,i],]
  
  fit = pffr(d2A ~ A + dA, data = pda.data, yind = tfine, 
             bs.yindex = list(k = 101, m = c(3,2), bs = 'cc', sp = 1e8),
             bs.int = list(k = 101, m = c(3,2), bs = 'cc', sp = 1e8))
  
  coefs = extract.coefs(fit)
  
  raw = lapply(coefs, function(c){c$value})
  
  return(raw)
  
}


# set up cluster
n_cores <- 10
cl <- makeCluster(n_cores)
registerDoParallel(cl)


clusterExport(
  cl,
  c("pffr.data.a", "boot.ind", "extract.coefs", "boot.fun"),
  envir = environment()
)


# run bootstrap
boot.mats <- foreach(L = 1:B, #.combine = 'c',
                   .packages = c("refund", "mgcv"),
                   .combine = function(x, y) {
                                        list(
                                          int = cbind(x[[1]], y[[1]]),
                                          beta_A = cbind(x[[2]], y[[2]]),
                                          beta_dA = cbind(x[[3]], y[[3]]))},
                                      .init = list(
                                        int = NULL,
                                        beta_A = NULL,
                                        beta_dA = NULL
                                      )) %dopar% {boot.fun(L)}
                     
                     
                     
                   



stopCluster(cl)



# update se

ankle.co[[1]]$seb = apply(boot.mats$int, 1, sd)
ankle.co[[2]]$seb = apply(boot.mats$beta_A, 1, sd)
ankle.co[[3]]$seb = apply(boot.mats$beta_dA, 1, sd)


# plot coefficients

Stiff = ggplot(ankle.co[[2]], aes(x = tfine.vec*100, y = value, 
                                 ymin = value-2*seb, ymax = value+2*seb))+
  geom_ribbon(alpha = 0.4, fill = "orange")+
  geom_line()+
  labs(x = "% Gait Cycle", y = "Coefficient Value", title = "A")+
  geom_hline(yintercept = 0, lty = 2)+
  theme_light()+
  theme(text = element_text(size = 14))

Damp = ggplot(ankle.co[[3]], aes(x = tfine.vec*100, y = value, 
                                ymin = value-2*seb, ymax = value+2*seb))+
  geom_ribbon(alpha = 0.4, fill = "blue")+
  geom_line()+
  labs(x = "% Gait Cycle", y = " ", title = "B")+
  geom_hline(yintercept = 0, lty = 2)+
  theme_light()+
  theme(text = element_text(size = 14))


COEFPLOT = Stiff + Damp
COEFPLOT

#ggsave("FLDHM Figure 3.2.Ankle.jpg", COEFPLOT, dpi = 300, height = 10, width = 18, units = 'cm')




# create a function for each coefficient
afun = lapply(ankle.co, function(X){splinefun(x = X$tfine.vec, y = X$value)})


# construct the state transition matrix
A_t = function(t){
  
  A = matrix(c(0,               1, 
               -afun[[2]](t),  -afun[[3]](t)),
             ncol = 2, nrow = 2, byrow = T)
  
  return(A)
  
}

# construct the Bu components
B = matrix(c(0,1), nrow = 2)

Bu_t = function(t){
  
  U = matrix(
    c(0, afun[[1]](t)),
    ncol = 1)
  
  Bu = B*U
  return(Bu)  
}

# define the ODE
LL = function(t,y,params){
  
  
  dy = A_t(t)%*%y + Bu_t(t)
  
  return(list(dy))
}

# quick test, this shold return functions very close to zero
initial = c(0,0)

solution <- ode(y = initial, times = tfine, func = LL, parms = NULL, method = "lsoda")

plot(solution)





# look at reconstruction of data

time = tfine

inits = cbind(t(eval.fd(0, ankle.fd.c)),
              t(eval.fd(0, dankle.fd.c)))

colnames(inits) = c("A", "dA")
rownames(inits) = 1:n

a.mat = c()
da.mat = c()


for(i in 1:nrow(inits)){
  
  initial = inits[i,]
  
  solution <- ode(y = initial, times = time, func = LL, parms = NULL, method = "lsoda")
  
  a.mat = cbind(a.mat, solution[,2])
  da.mat = cbind(da.mat, solution[,3])
}


# add back in the mean

a.hat = apply(a.mat, 2, function(x){x+eval.fd(time, mean.fd(ankle.fd))})
da.hat = apply(da.mat, 2, function(x){x+eval.fd(time, mean.fd(dankle.fd))})


par(mfrow = c(1,3))
matplot(a.hat, type = 'l', main = "Esimated Ankle", ylim = c(-40,40))
matplot(eval.fd(time, ankle.fd), type = 'l', main = "True Ankle", ylim = c(-40,40))
matplot(eval.fd(time, ankle.fd)-a.hat, type = 'l', main = "Ankle Residuals")




# look at eigen values of A(t)

eigmat = matrix(NA, nrow = 101, ncol = 2)
eigvecs = array(NA, c(101,2,2))

for(i in 1:101){
  
  A = A_t(time[i])
  
  EIG = eigen(A)
  
  eigmat[i,] = EIG$values
  
  eigvecs[i,,1] = t(EIG$vectors[,1])
  eigvecs[i,,2] = t(EIG$vectors[,2])
  
  
  
}

# extract just the real component and plot
R.eigmat = apply(eigmat, 2, Re) 

real.plot = melt(R.eigmat)

TVES = ggplot(real.plot, aes(x = X1, y = value, colour = as.factor(X2)))+
  geom_line()+
  geom_hline(yintercept = 0, lty = 2)+
  labs(x = " ", y  ="Eigen Value \n (Real Part)")+
  theme_light()+
  facet_wrap( ~as.factor(X2),ncol = 6, scales = 'free')+
  theme(legend.position = 'none',
        text = element_text(size = 16))
TVES

# look at eigen vectors

evec1 = Re(eigvecs[,,1])
colnames(evec1) =  c("A", "dA")
evc1.plot = melt(evec1)

evec2 = Re(eigvecs[,,2])
colnames(evec2) =  c("A", "dA")
evc2.plot = melt(evec2)



EVPLOT = data.frame(rbind(evc1.plot, evc2.plot))
EVPLOT$Eigenvector = as.factor(c(rep(1:2, each = 202)))


my_palette <- c(
  "#0072B2", # Blue
  "#D55E00")  # Orange


TVCS = ggplot(EVPLOT, aes(x = X1, y = abs(value), colour = X2))+
  geom_line()+
  labs(x = "% Gait Cycle", y = "Eigen Vector \n (Real Part)", colour = "Variable")+
  theme_light()+
  facet_wrap(~Eigenvector, ncol = 6)+
  scale_colour_manual(values = my_palette)+
  theme(legend.position = 'bottom', text = element_text(size = 16))

TVCS
layout = "AAA
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
          CCC"

EIGVALVEC = wrap_plots(TVES, TVCS, guide_area())+ plot_layout(design = layout)
EIGVALVEC
#ggsave("FLDHM Figure 5.3 - Ankle Sensitivity.jpg", EIGVALVEC, dpi = 300, height = 16, width = 15, units = 'cm')






