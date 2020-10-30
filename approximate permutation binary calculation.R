rm(list=ls())
library(MASS)
library(hier.part)
library(invgamma)
library(survival)
library(tictoc)
library(parallel)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)
########################Part 1: Data generation #######################
data_gene = function(dataid,rhoxu = 0,rhoyu = 0,m = 50,effect_true=2){
  n = as.matrix(sample(22:40,m,replace=T,prob=rep(1/19,19)))
  id = (unlist(apply(n,1,seq,from=1,by=1)))
  clu = rep(1:m,n)
  data = as.data.frame(cbind(clu,id))
  N = dim(data)[1]
  inter = rep(1,N)
  data$Int = inter
  data$x1 = rnorm(N,0,1)
  data$x2 = sample(c(-1,0,1),N,replace=T,prob=rep(1/3,3))
  x1mean = rep(-99,m) 
  x2mean = rep(-99,m)
  for (i in 1:m){
    x1mean[i] = mean(data[data$clu==i,]$x1)
    x2mean[i] = mean(data[data$clu==i,]$x2)
  }
  u = apply(as.matrix(-rhoxu*(x1mean+x2mean)),1,rnorm,n=1,sd=1)
  data$U = rep(u,n)
  ps=function(para){# x1,x2,U
    exp(para[1]+para[2]+para[3])/(1+exp(para[1]+para[2]+para[3]))
  }
  ps_true = apply(data[,4:6],1,ps)
  aij = apply(as.matrix(ps_true),1,rbinom,n=1,size=1)
  data$A = aij
  suff_T = rep(-99,m)
  for (i in 1:m){
    suff_T[i] = sum(data[data$clu==i,]$A)
  }
  while (suff_T[1]==n[1] | suff_T[1]==0){
    sap = apply(as.matrix(ps_true[1:n[1]]),1,rbinom,n=1,size=1)
    data$A[1:n[1]] = sap
    suff_T[1] = sum(sap)
  }
  begin = 1
  total = n[1]
  for (i in 2:m){
    begin = begin + n[i-1]
    total = total+ n[i]
    while (suff_T[i]==n[i] | suff_T[i]==0){
      sap = apply(as.matrix(ps_true[begin:total]),1,rbinom,n=1,size=1) 
      data$A[begin:total] = sap
      suff_T[i] = sum(sap)
    }
  }
  data$sufft = rep(suff_T,n)
  data$pstrue = ps_true
  data$y0 = 1 + data$x1 + data$x2 + rhoyu*data$U+ rnorm(N,0,1)
  data$y1 = 1 + data$x1 + data$x2 + effect_true + rhoyu*data$U + rnorm(N,0,1)
  data$y_obs = data$y1*data$A + data$y0*(1-data$A)
  effec_simu = mean(data$y1-data$y0)
  effec_naiv = mean(data$A*data$y_obs-(1-data$A)*data$y_obs)
  data_obs = data[,c(1:5,7,8,12)]
  return(list(data=data,data_obs=data_obs,N=N,m=m,n=n,u=u,suff_T = suff_T,effec_simu=effec_simu,effec_naiv=effec_naiv))
}
simudata = mclapply(1:3,data_gene,rhoxu = 0,rhoyu = 0,m = 50,effect_true=,mc.cores=1)

####################### formula based #################
t1 = proc.time()
psfunc_suffi = function(dataid,simudata){
  data = simudata[[dataid]]$data
  data_obs = simudata[[dataid]]$data_obs
  N = simudata[[dataid]]$N
  m = simudata[[dataid]]$m
  n = simudata[[dataid]]$n
  u = simudata[[dataid]]$u
  suff_T = simudata[[dataid]]$suff_T
  
  clog_summ = clogit(A~x1+x2+strata(clu), data=data)
  beta_est = clog_summ$coefficients
  dem=rep(-99,m)
  prop_ps = 0
  for (cls in 1:m){
    ## denominator part
    data_clu = as.matrix(data[data$clu==cls,4:5])
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
      if (data[data$clu==cls,]$A[sub]==0){
        nume= rnum[suffT+1,ni+1]
      }else{
        nume= dem[cls]-rnum[suffT+1,ni+1]
      }
      
      prob = nume/dem[cls]
      prop_ps = rbind(prop_ps,prob)
    }
  }
  ps_suffi=prop_ps[-1]

  print(dataid)
  return(ps_suffi)
}
suffieffect = mclapply(1:3,psfunc_suffi,simudata,mc.cores=1)
t11 = proc.time()-t1

t2 = proc.time()
psfunc2_suffi = function(dataid,simudata){
   data = simudata[[dataid]]$data
   data_obs = simudata[[dataid]]$data_obs
   N = simudata[[dataid]]$N
   m = simudata[[dataid]]$m
   n = simudata[[dataid]]$n
   u = simudata[[dataid]]$u
   suff_T = simudata[[dataid]]$suff_T
   
   clog_summ = clogit(A~x1+x2+strata(clu), data=data)
   beta_est = clog_summ$coefficients

   binary_suffi_func=function(cls,data){
     prop_ps = c()
     dem=rep(-99,m)
     ## denominator part
     data_clu = as.matrix(data[data$clu==cls,4:5])
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
       if (data[data$clu==cls,]$A[sub]==0){
         nume= rnum[suffT+1,ni+1]
       }else{
         nume= dem[cls]-rnum[suffT+1,ni+1]
       }
       prob = nume/dem[cls]
       prop_ps = c(prop_ps,prob)
     }
     return(prop_ps)
   }
   ps_suffi2_list = lapply(c(1:m),binary_suffi_func,data_obs)
   ps_suffi2 = unlist(ps_suffi2_list)
   return(ps_suffi2)
}
suffieffect2 = mclapply(1:3,psfunc2_suffi,simudata,mc.cores=1)
t22 = proc.time()-t2

############# resampling method ########################
a = proc.time()
func_suffi = function(dataid,simdata,resamp=2000){
  set.seed(dataid)
  data = simdata[[dataid]]$data
  n = simdata[[dataid]]$n
  m = simdata[[dataid]]$m
  suff_T = simdata[[dataid]]$suff_T
  
  clog_summ = clogit(A~x1+x2+strata(clu), data=data)
  beta_est = clog_summ$coefficients
  ####### construct Vij
  prop_ps = 0
  for (i in 1:m){
    trt_com = data[data$clu==i,]$A # trt assignments within cluster
    bbb = t(replicate(resamp,sample(trt_com,n[i],FALSE)))
    denomi = sum(exp(bbb%*%as.matrix(data[data$clu==i,4:5])%*%beta_est))
    for (j in 1:n[i]){
      aaa = bbb[bbb[,j]==data[data$clu==i,]$A[j],]
      nume = exp(data[data$clu==i,]$A[j]*as.matrix(data[data$clu==i,4:5][j,])%*%beta_est)*sum(exp(aaa[,-j]%*%(as.matrix(data[data$clu==i,4:5])[-j,]%*%beta_est)))
      
      prob = nume/denomi
      prop_ps = rbind(prop_ps,prob)
    }
    print(i)
  }
  propps = prop_ps[-1,]
  return(propps)
}
suffieffect3 = mclapply(1:3,func_suffi,simudata,resamp=2000,mc.cores=1)
b = proc.time()-a
print(b)

aa = proc.time()
resamp_suffi = function(dataid,simdata,resamp=1000){
  set.seed(dataid)
  data = simdata[[dataid]]$data
  n = simdata[[dataid]]$n
  m = simdata[[dataid]]$m
  suff_T = simdata[[dataid]]$suff_T
  
  clog_summ = clogit(A~x1+x2+strata(clu), data=data)
  beta_est = clog_summ$coefficients
  binary_suffi_samp = function(cls,data){
    prop_ps = c()
    trt_com = data[data$clu==cls,]$A # trt assignments within cluster
    bbb = t(replicate(resamp,sample(trt_com,n[cls],FALSE)))
    denomi = sum(exp(bbb%*%as.matrix(data[data$clu==cls,4:5])%*%beta_est))
    for (j in 1:n[cls]){
      aaa = bbb[bbb[,j]==data[data$clu==cls,]$A[j],]
      nume = exp(data[data$clu==cls,]$A[j]*as.matrix(data[data$clu==cls,4:5][j,])%*%beta_est)*sum(exp(aaa[,-j]%*%(as.matrix(data[data$clu==cls,4:5])[-j,]%*%beta_est)))
      
      prob = nume/denomi
      prop_ps = c(prop_ps,prob)
    }
    return(prop_ps)
  }
  ps_suffi3 = unlist(lapply(c(1:m),binary_suffi_samp,data))
  return(ps_suffi3)
}
suffieffect4 = mclapply(1:3,resamp_suffi,simudata,resamp=2000,mc.cores=1)
bb = proc.time()-aa
print(bb)


plot(suffieffect[[1]],suffieffect2[[1]])
plot(suffieffect[[1]],suffieffect5[[1]])
plot(suffieffect[[1]],suffieffect3[[1]])
plot(suffieffect5[[1]],suffieffect4[[1]])

########### estimate results ################
data = simudata[[1]]$data
ps1 = suffieffect[[1]]
ps2 = suffieffect2[[1]]
ps3 = suffieffect3[[1]]
ps4 = suffieffect4[[1]]
ps5 = suffieffect5[[1]]
data$w_suffi1 = ifelse(data$A==1,ps1,1-ps1)
data$w_suffi2 = ifelse(data$A==1,ps2,1-ps2)
data$w_suffi3 = ifelse(data$A==1,ps3,1-ps3)
data$w_suffi4 = ifelse(data$A==1,ps4,1-ps4)
data$w_suffi5 = ifelse(data$A==1,ps5,1-ps5)
effect_suffi1 = mean(data$A*data$y_obs/data$w_suffi1-(1-data$A)*data$y_obs/(1-data$w_suffi1))
effect_suffi2 = mean(data$A*data$y_obs/data$w_suffi2-(1-data$A)*data$y_obs/(1-data$w_suffi2))
effect_suffi3 = mean(data$A*data$y_obs/data$w_suffi3-(1-data$A)*data$y_obs/(1-data$w_suffi3))
effect_suffi4 = mean(data$A*data$y_obs/data$w_suffi4-(1-data$A)*data$y_obs/(1-data$w_suffi4))
effect_suffi5 = mean(data$A*data$y_obs/data$w_suffi5-(1-data$A)*data$y_obs/(1-data$w_suffi5))
