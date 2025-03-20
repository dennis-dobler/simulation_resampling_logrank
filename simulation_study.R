library(survival)

lr.test <- function(dN1, dN2, Y1, Y2){
  test.stat.num <- as.numeric((Y2/(Y1+Y2)) %*% dN1 - (Y1/(Y1+Y2)) %*% dN2)
  test.stat.denom <- as.numeric((Y2^2/(Y1+Y2)^2 ) %*% dN1 + (Y1^2/(Y1+Y2)^2) %*% dN2)
  test.stat <- test.stat.num^2 / test.stat.denom
  return(test.stat)
}

lr.test.perm <- function(datf, n){
  datf.perm <- sample(datf[,3])
  dN1 <- datf[,2] * (2 - datf.perm)
  dN2 <- datf[,2] * (datf.perm - 1)
  Y1 <- c(n*40/100, n*40/100 - cumsum(datf.perm==1))[-(n+1)]
  Y2 <- c(n*60/100, n*60/100 - cumsum(datf.perm==2))[-(n+1)]
  return(lr.test(dN1, dN2, Y1, Y2))
}

lr.test.pool.bs <- function(datf, n){
  bs.ind <- sample(1:n, replace=T)
  datf.pool.bs <- datf
  datf.pool.bs[,1] <- datf[bs.ind,1]
  datf.pool.bs[,2] <- datf[bs.ind,2]
  datf.pool.bs <- datf.pool.bs[order(datf.pool.bs[,1]),]
  dN1 <- datf.pool.bs[,2] * (2 - datf.pool.bs[,3])
  dN2 <- datf.pool.bs[,2] * (datf.pool.bs[,3] - 1)
  # As if ties were broken
  Y1 <- c(n*40/100, n*40/100 - cumsum(datf.pool.bs[,3]==1))[-(n+1)]
  Y2 <- c(n*60/100, n*60/100 - cumsum(datf.pool.bs[,3]==2))[-(n+1)]
  return(lr.test(dN1, dN2, Y1, Y2))
}


lr.test.wild <- function(dN1, dN2, Y1, Y2, n){
  G <- (rpois(n,lambda=1)-1)
  test.stat.num <- as.numeric((Y2/(Y1+Y2)) %*% (G * dN1) - (Y1/(Y1+Y2)) %*% (G * dN2))
  test.stat.denom <- sqrt(as.numeric((Y2^2/(Y1+Y2)^2 ) %*% (G^2 * dN1) + (Y1^2/(Y1+Y2)^2) %*% (G^2 * dN2)))
  test.stat <- (test.stat.num / test.stat.denom)^2
  if(is.nan(test.stat)) test.stat <- 0
  return(test.stat)
}


one.test <- function(itnr, n, distr.ind, c.s){
  times1 <- distr.vec[[distr.ind]](n)
  groups <- c(rep(1, n*40/100), rep(2, n*60/100))
  times <- pmin(times1, c(runif(n*40/100, min=0.1*c.s, max=4*c.s), runif(n*60/100, min=1*c.s, max=3*c.s)))
  deltas <- times==times1
  
  times.sort <- sort(times)
  
  datf <- data.frame(times, deltas, groups)
  datf <- datf[order(times),]
  
  dN1 <- datf[,2] * (2 - datf[,3])
  dN2 <- datf[,2] * (datf[,3] - 1)
  
  Y1 <- c(n*40/100, n*40/100 - cumsum(datf[,3]==1))[-(n+1)]
  Y2 <- c(n*60/100, n*60/100 - cumsum(datf[,3]==2))[-(n+1)]
  
  lr.asy <- lr.test(dN1=dN1, dN2=dN2, Y1=Y1, Y2=Y2)
  lr.asy.p <- 1 - pchisq(lr.asy, df=1)
  
  lr.perm <- replicate(B, lr.test.perm(datf, n))
  lr.perm <- lr.perm[which(!is.na(lr.perm))]
  lr.perm.p <- mean(lr.perm >= lr.asy)
  
  lr.pool.bs <- replicate(B, lr.test.pool.bs(datf, n))
  lr.pool.bs <- lr.pool.bs[which(!is.na(lr.pool.bs))]
  lr.pool.bs.p <- mean(lr.pool.bs >= lr.asy)
  
  lr.wild <- replicate(B, lr.test.wild(dN1, dN2, Y1, Y2, n))
  lr.wild.p <- mean(lr.wild >= lr.asy)
  
  return(c(lr.asy.p <= 0.05, lr.perm.p <= 0.05, lr.pool.bs.p <= 0.05, lr.wild.p <= 0.05))
}


setwd("/home/dobler/Dokumente/Konferenzen und Dienstreisen/2025/SAfJR/talk/R")

set.seed(210325)
K <- 2000
B <- 1000

n.vec <- c(25,50,100,200)
distr.vec <- list(function(n) rexp(n, rate=1), function(n) rweibull(n, shape=2, scale=0.5), function(n) rlnorm(n, meanlog=-0.5, sdlog=0.5))
cens.scales <- c(1,1.5)

results <- array(numeric(length(n.vec)* length(distr.vec)* length(cens.scales)*4), dim = c(length(n.vec), length(distr.vec), length(cens.scales), 4))

for(n in n.vec){
  for(distr.ind in 1:length(distr.vec)){
    for(c.s in cens.scales){
      results[which(n.vec==n), distr.ind, which(cens.scales==c.s), ] <- rowMeans(sapply(1:K, one.test, n, distr.ind, c.s))
      print(c(n, distr.ind, c.s))
      print(results[which(n.vec==n), distr.ind, which(cens.scales==c.s), ])
      if(T){
        L <- list(res = results[,,,1])
        save(L, file="res.asy.rdata")
        L <- list(res = results[,,,2])
        save(L, file="res.perm.rdata")
        L <- list(res = results[,,,3])
        save(L, file="res.pool.rdata")
        L <- list(res = results[,,,4])
        save(L, file="res.wild.rdata")
      }
    }
  }
}

if(F){
  L <- list(res = results[,,,1])
  save(L, file="res.asy.rdata")
  L <- list(res = results[,,,2])
  save(L, file="res.perm.rdata")
  L <- list(res = results[,,,3])
  save(L, file="res.pool.rdata")
  L <- list(res = results[,,,4])
  save(L, file="res.wild.rdata")
}


load("res.asy.rdata")
asy <- L$res
load("res.perm.rdata")
perm <- L$res
load("res.pool.rdata")
pool <- L$res
load("res.wild.rdata")
wild <- L$res
boxplot(c(asy), c(perm), c(pool), c(wild), names=c("asy", "perm", "pooled bs", "wild"), ylim=c(0.035,0.065))
abline(h=0.05, col=2, lwd=3)


load("res.asy.rdata")
asy <- L$res[1,,]
load("res.perm.rdata")
perm <- L$res[1,,]
load("res.pool.rdata")
pool <- L$res[1,,]
load("res.wild.rdata")
wild <- L$res[1,,]
boxplot(c(asy), c(perm), c(pool), c(wild), names=c("asy", "perm", "pooled bs", "wild"), ylim=c(0.035,0.065))
abline(h=0.05, col=2, lwd=3)



load("res.asy.rdata")
asy <- L$res[2,,]
load("res.perm.rdata")
perm <- L$res[2,,]
load("res.pool.rdata")
pool <- L$res[2,,]
load("res.wild.rdata")
wild <- L$res[2,,]
boxplot(c(asy), c(perm), c(pool), c(wild), names=c("asy", "perm", "pooled bs", "wild"), ylim=c(0.035,0.065))
abline(h=0.05, col=2, lwd=3)



load("res.asy.rdata")
asy <- L$res[3,,]
load("res.perm.rdata")
perm <- L$res[3,,]
load("res.pool.rdata")
pool <- L$res[3,,]
load("res.wild.rdata")
wild <- L$res[3,,]
boxplot(c(asy), c(perm), c(pool), c(wild), names=c("asy", "perm", "pooled bs", "wild"), ylim=c(0.035,0.065))
abline(h=0.05, col=2, lwd=3)



load("res.asy.rdata")
asy <- L$res[4,,]
load("res.perm.rdata")
perm <- L$res[4,,]
load("res.pool.rdata")
pool <- L$res[4,,]
load("res.wild.rdata")
wild <- L$res[4,,]
boxplot(c(asy), c(perm), c(pool), c(wild), names=c("asy", "perm", "pooled bs", "wild"), ylim=c(0.035,0.065))
abline(h=0.05, col=2, lwd=3)