rm(list=ls())
library(MASS)
library(invgamma)
library(survival)
library(parallel)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)

######### data generation ##################
data_gene = function(dataid,m=10,m0=6,rhoxu=0,rhoyu=0,effect_true=2){
  #data_gene = function(dataid,m=500,rhoxu=0,rhoyu=0,effect_true=2){
  set.seed(dataid)
  n = as.matrix(sample(10:50,m,replace=T,prob=rep(1/41,41)))
  #n = as.matrix(sample(100:200,m,replace=T,prob=rep(1/101,101)))
  id = unlist(apply(n,1,seq,from=1,by=1))
  clu = rep(1:m,n)
  data = as.data.frame(cbind(clu,id))
  N = dim(data)[1]
  data$x1 = rnorm(N,0,1)
  data$x2 = sample(c(-1,0,1),N,replace=T,prob=rep(1/3,3))
  x1mean = rep(-99,m) 
  x2mean = rep(-99,m)
  for (i in 1:m){
    x1mean[i] = mean(data[data$clu==i,]$x1)
    x2mean[i] = mean(data[data$clu==i,]$x2)
  }
  u = apply(as.matrix(-rhoxu*(x1mean+x2mean)),1,rnorm,n=1,sd=1)
  u0 = rnorm(m0,0,0.01)
  u[1:m0]=u0
  data$U = rep(u,n)
  psfunc=function(x,para){
    1/(1+exp(-x%*%para))
  }
  data$pstrue = apply(as.matrix(data[,3:5]),1,psfunc,para=c(1,1,1))
  
  ### remove too small/large pscore
  data = data[data$pstrue>=0.05 & data$pstrue<=0.95,]
  m = length(unique(data$clu))
  n = as.numeric(table(data$clu)) 
  n_ind = which(n>1)
  m_vec = c(1:m)[n_ind]
  data = data[data$clu %in%m_vec,]
  N = dim(data)[1]
  u = unique(data$U)
  m = length(unique(data$clu))
  n = as.numeric(table(data$clu)) 
  data$clu = rep(1:m,n)
  data$id = unlist(apply(matrix(n),1,seq,from=1,by=1))
  
  zij = apply(matrix(data$pstrue),1,rbinom,n=1,size=1)
  data$z = zij
  suff_T = rep(-99,m)
  for (i in 1:m){
    suff_T[i] = sum(data[data$clu==i,]$z)
  }
  ##### check each cluster, must have treatment and control both ####
  while (suff_T[1]==n[1] | suff_T[1]==0){
    sap = apply(as.matrix(data$pstrue[1:n[1]]),1,rbinom,n=1,size=1)
    data$z[1:n[1]] = sap
    suff_T[1] = sum(sap)
  }
  begin = 1
  total = n[1]
  for (i in 2:m){
    begin = begin + n[i-1]
    total = total+ n[i]
    while (suff_T[i]==n[i] | suff_T[i]==0){
      sap = apply(as.matrix(data$pstrue[begin:total]),1,rbinom,n=1,size=1) 
      data$z[begin:total] = sap
      suff_T[i] = sum(sap)
    }
  }
  data$sufft = rep(suff_T,n)
  data$y0 = 1 + data$x1 + data$x2 + rnorm(N,0,1)
  data$y1 = 1 + data$x1 + data$x2 + effect_true + rhoyu*data$U + rnorm(N,0,1)
  data$y_obs = data$y1*data$z + data$y0*(1-data$z)
  
  data_obs = data[,c(1:4,7,8,11)] # include int
  print(dataid)
  return (list(data=data,data_obs=data_obs,u=u,m=m,n=n,suff_T=suff_T))
}
simdatam20= mclapply(1:2,data_gene,m=10,m0=6,rhoxu=0,rhoyu=0,effect_true=2,mc.cores=1)

########## sufficient formula ###############
psfunc2_suffi = function(dataid,simudata){
  ptm = proc.time()
  set.seed(dataid)
  data = simudata[[dataid]]$data
  data_obs = simudata[[dataid]]$data_obs
  N = simudata[[dataid]]$N
  m = simudata[[dataid]]$m
  n = simudata[[dataid]]$n
  suff_T = simudata[[dataid]]$suff_T
  
  clog_summ = clogit(z~x1+x2+strata(clu), data=data)
  beta_est = clog_summ$coefficients
  
  binary_suffi_func=function(cls,data){
    prop_ps = c()
    dem=rep(-99,m)
    ## denominator part
    data_clu = as.matrix(data[data$clu==cls,3:4])
    ni=n[cls]
    suffT =suff_T[cls]
    rdenom = matrix(-99,nrow=suffT+1,ncol=ni+1)
    rdenom[lower.tri(rdenom)]=0 ## lower part =0
    rdenom[1,]=1
    for (rowid in 2:(suffT+1)){
      for (colid in 2:(ni+1)){
        rdenom[rowid,colid]=rdenom[rowid,colid-1] + exp(data_clu[colid-1,]%*%beta_est)*rdenom[rowid-1,colid-1]
      }
    }
    dem[cls] = rdenom[suffT+1,ni+1]
    for (sub in 1:n[cls]){
      ## numerator part
      rnum = matrix(-99,nrow=suffT+1,ncol=ni+1)
      rnum[lower.tri(rnum)]=0 ## lower part =0
      rnum[1,]=1
      for (rowid in 2:(suffT+1)){
        for (colid in 2:(ni+1)){
          rnum[rowid,colid]=rnum[rowid,colid-1] + exp(data_clu[colid-1,]%*%beta_est)*(1-as.numeric(colid-1==sub))*rnum[rowid-1,colid-1]
        }
      }
      if (data[data$clu==cls,]$z[sub]==0){
        nume= rnum[suffT+1,ni+1]
      }else{
        nume= dem[cls]-rnum[suffT+1,ni+1]
      }
      prob = nume/dem[cls]
      prop_ps = c(prop_ps,prob)
    }
    return(prop_ps)
  }
  data$ps_weight = unlist(lapply(c(1:m),binary_suffi_func,data_obs))
  effect_suffi = mean(data$z*data$y_obs/data$ps_weight-(1-data$z)*data$y_obs/(data$ps_weight))
  print(dataid)
  time_span = proc.time()-ptm
  return(list(effect_suffi=effect_suffi,time_span=time_span))
}
effect_suffim20 = mclapply(1:2,psfunc2_suffi,simdatam20,mc.cores=1)
