rm(list=ls())
library(parallel)
library(nnet)
library(Matrix)
library(mclogit)
library(survival)
library(gtools)
########################Part 1: Data generation #######################
data_gene = function(dataid,rhoxu1 = 0,rhoxu2 = 0,rhoyu1 = 0,rhoyu2 = 0,m = 100){
  n = as.matrix(sample(3:5,m,replace=T,prob=rep(1/3,3)))
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
  u1 = apply(as.matrix(-rhoxu1*(x1mean+x2mean)),1,rnorm,n=1,sd=1)
  u2 = apply(as.matrix(-rhoxu2*(x1mean+x2mean)),1,rnorm,n=1,sd=1)
  data$U1 = rep(u1,n)
  data$U2 = rep(u2,n)
  ps0 = function(para){
    1/(exp(para[1]+para[2]+para[3])+exp(para[1]+para[2]+para[4])+1)
  }
  ps1 = function(para){
    exp(para[1]+para[2]+para[3])/(exp(para[1]+para[2]+para[3])+exp(para[1]+para[2]+para[4])+1)
  }
  pstrue = matrix(-99,ncol=3,nrow=N)
  pstrue[,1] = apply(data[,4:7],1,ps0)
  pstrue[,2] = apply(data[,4:7],1,ps1)
  pstrue[,3] = 1-pstrue[,1]- pstrue[,2]
  aij = apply(pstrue,1,sample,x=c(0,1,2),size=1,replace=F)
  data$A = aij
  suffi = matrix(-99,nrow=m,ncol=2)
  for (i in 1:m){
    cluA = data[data$clu==i,]$A
    suffi[i,1] = length(cluA[cluA==1])
    suffi[i,2] = length(cluA[cluA==2])
  }
  while (suffi[1,1]==0|suffi[1,2]==0|sum(suffi[1,])==n[1] | sum(suffi[1,])==0){
    sap = apply(as.matrix(pstrue[1:n[1],]),1,sample,x=c(0,1,2),size=1,replace=F)
    data$A[1:n[1]] = sap
    suffi[1,] = c(length(sap[sap==1]),length(sap[sap==2]))
  }
  begin = 1
  total = n[1]
  for (i in 2:m){
    begin = begin + n[i-1]
    total = total+ n[i]
    while (suffi[i,1]==0|suffi[i,2]==0|sum(suffi[i,])==n[i] | sum(suffi[i,])==0){
      sap = apply(as.matrix(pstrue[begin:total,]),1,sample,x=c(0,1,2),size=1,replace=F) 
      data$A[begin:total] = sap
      suffi[i,] = c(length(sap[sap==1]),length(sap[sap==2]))
    }
  }
  data$sufft1 = rep(suffi[,1],n)
  data$sufft2 = rep(suffi[,2],n)
  data$y0 = 1 + data$x1 + data$x2 + rhoyu1*data$U1 +rhoyu2*data$U2+rnorm(N,0,1)
  data$y1 = 1 + data$x1 + data$x2 + 2 + rhoyu1*data$U1 +rhoyu2*data$U2 + rnorm(N,0,1)
  data$y2 = 1 + data$x1 + data$x2 + 4 + rhoyu1*data$U1 +rhoyu2*data$U2 + rnorm(N,0,1)
  data$y_obs = ifelse(data$A==0,data$y0,ifelse(data$A==1,data$y1,data$y2))
  effec_true = c(2,4,2) #effect 01,02,12
  effec_simu = c(mean(data$y1-data$y0),mean(data$y2-data$y0),mean(data$y2-data$y1))
  print(dataid)
  return(list(data=data,m=m,n=n,N=N,suffi=suffi,pstrue=pstrue,effec_simu=effec_simu))
}
simudata = mclapply(1:1,data_gene,rhoxu1 = 5,rhoxu2 = 5,rhoyu1 = 5,rhoyu2 = 5,m = 500,mc.cores=1)
################### Part 2: PS estimation  ###################
## naive way 
effectnaivefunc = function(dataid){
  data=simudata[[dataid]]$data
  N = simudata[[dataid]]$N
  n = simudata[[dataid]]$n
  m = simudata[[dataid]]$m
  suffi = simudata[[dataid]]$suffi
  y0_naive = mean(data[data$A==0,]$y_obs)
  y1_naive = mean(data[data$A==1,]$y_obs)
  y2_naive = mean(data[data$A==2,]$y_obs)
  effect_naive = c(y1_naive-y0_naive,y2_naive-y0_naive,y2_naive-y1_naive)
  print(dataid)
  return(effect_naive=effect_naive)
}
res_naive = mclapply(1:1,effectnaivefunc,mc.cores=1)

## fix effect 
effectfixfunc = function(dataid){
  data=simudata[[dataid]]$data
  N = simudata[[dataid]]$N
  n = simudata[[dataid]]$n
  m = simudata[[dataid]]$m
  suffi = simudata[[dataid]]$suffi
  ####################### pscore estimation ########################
  data$A_recode = relevel(factor(data$A), ref = "0")
  logfix = multinom(A_recode~x1+x2+factor(clu),data=data,MaxNWts =10000000)  # if factor clu, error for too many weights, so give large wts
  w_fix = predict(logfix,type="probs") ## or logfix[["fitted.values"]]
  data$w_fix0 = w_fix[,1]
  data$w_fix1 = w_fix[,2]
  data$w_fix2 = w_fix[,3]
  data$w_fix = ifelse(data$A==0,data$w_fix0,ifelse(data$A==1,data$w_fix1,data$w_fix2))
  data0 = data[data$A==0,]
  data1 = data[data$A==1,]
  data2 = data[data$A==2,]
  # y0_fix = mean(data0$y_obs/data0$w_fix)
  y0_fix = sum(data0$y_obs/data0$w_fix)/sum(1/data0$w_fix)
  y1_fix = sum(data1$y_obs/data1$w_fix)/sum(1/data1$w_fix)
  y2_fix = sum(data2$y_obs/data2$w_fix)/sum(1/data2$w_fix)
  effect_fix = c(y1_fix-y0_fix,y2_fix-y0_fix,y2_fix-y1_fix)
  print(dataid)
  return (effect_fix=effect_fix)
}
res_fix = mclapply(1:1,effectfixfunc,mc.cores=1)
## method 1 
# logfix = mblogit(A_recode~x1+x2+factor(clu),random=NULL, data=data)
# w_fix = matrix(logfix$fitted.values,byrow=T,nrow=N,ncol=3)
## method 2 

## random model 
# library(lme4)
# logran = lmer(A ~ x1 + x2 + 1|clu,data=data)
# predict(logran)
# library(mlogit)
effectranfunc = function(dataid){
  data=simudata[[dataid]]$data
  N = simudata[[dataid]]$N
  n = simudata[[dataid]]$n
  m = simudata[[dataid]]$m
  suffi = simudata[[dataid]]$suffi
  ####################### pscore estimation ########################
  data$A_recode = relevel(factor(data$A), ref = "0")
  ran_summ = mblogit(A_recode~x1+x2,random=~1|clu, data=data)
  w_ran = matrix(ran_summ$fitted.values,byrow=T,nrow=N,ncol=3)
  data$w_ran0 = w_ran[,1]
  data$w_ran1 = w_ran[,2]
  data$w_ran2 = w_ran[,3]
  data$w_ran = ifelse(data$A==0,data$w_ran0,ifelse(data$A==1,data$w_ran1,data$w_ran2))
  data0 = data[data$A==0,]
  data1 = data[data$A==1,]
  data2 = data[data$A==2,]
  y0_ran = sum(data0$y_obs/data0$w_ran)/sum(1/data0$w_ran)
  y1_ran = sum(data1$y_obs/data1$w_ran)/sum(1/data1$w_ran)
  y2_ran = sum(data2$y_obs/data2$w_ran)/sum(1/data2$w_ran)
  effect_ran = c(y1_ran-y0_ran,y2_ran-y0_ran,y2_ran-y1_ran)
  print(dataid)
  return(effect_ran=effect_ran)
}
res_ran = mclapply(1:1,effectranfunc,mc.cores=1)

# data_ran = data[,c(1,4,5,15)]
# data_ran$A_recode = ifelse(data_ran$A_recode=="0","g0",ifelse(data_ran$A_recode=="1","g1","g2"))
# data_ran$choice = ifelse(data_ran$A_recode=="g0","yes","no")
# dr = mlogit.data(data_ran,choice="choice",shape="wide")
# logran = mlogit(choice ~ x1 + x2,rpar=(clu="n"),data=dr)
# data$w_ran = predict(logran)
# effec_ran = mean(data$A*data$y_obs/data$w_ran-(1-data$A)*data$y_obs/(1-data$w_ran))


## sufficient statistics 
effectsuffifunc = function(dataid){
  data=simudata[[dataid]]$data
  N = simudata[[dataid]]$N
  n = simudata[[dataid]]$n
  m = simudata[[dataid]]$m
  suffi = simudata[[dataid]]$suffi
  ####################### pscore estimation ########################
  #data$choice = ifelse(data$A==0,1,0)
  #clog_summ = mclogit(cbind(A,choice)~x1+x2,random=~1|clu, data=data)
  #beta_mle = (clog_summ$coefficients)[1:2]
  beta1_mle = c(1,1,1) ## beta_0 as reference = 0
  beta2_mle = c(1,1,1)
  prop_ps = 0
  for (cl in 1:m){
    data_c = as.matrix(data[data$clu==cl,c(3,4:5,8)]) ## intercept, x1,x2,A
    for (sub in 1:n[cl]){
      data_num = data_c[-sub,1:3]
      data_den = data_c[,1:3]
      allposs = permutations(n=3,r=n[cl],v=c(0,1,2),repeats.allowed=T)
      Vbar = allposs
      Vij = allposs[allposs[,sub]==data_c[sub,4],]
      keep.crit = rep(-99,dim(Vij)[1])
      keep.crit2 = rep(-99,dim(Vbar)[1])
      for (i in 1:dim(Vij)[1]){
        t = Vij[i,]
        keep.crit[i] = ifelse(length(t[t==1])==suffi[cl,1] & length(t[t==2])==suffi[cl,2],"keep","discard")
      }
      ind = which(keep.crit=="keep")
      Vij = Vij[ind,]
      for (i in 1:dim(Vbar)[1]){
        t = Vbar[i,]
        keep.crit2[i] = ifelse(length(t[t==1])==suffi[cl,1] & length(t[t==2])==suffi[cl,2],"keep","discard")
      }
      ind2 = which(keep.crit2=="keep")
      Vbar = Vbar[ind2,]
      ########## calculate numerator
      #xij = data_c[sub,1:2]
      #d_need = c(xij%*%beta1_mle,xij%*%beta2_mle,1)
      sum_num = rep(-99,dim(Vij)[1])
      sum_astar = rep (-99,n[cl]-1)
      for (k in 1:dim(Vij)[1]){
        astar_j = Vij[k,-sub]
        Iasta = matrix(-99,nrow=n[cl]-1,ncol=3)
        for (q in 1:(n[cl]-1)){
          Iasta[q,1] = ifelse(astar_j[q]==1,1,0)
          Iasta[q,2] = ifelse(astar_j[q]==2,1,0)
          Iasta[q,3] = ifelse(astar_j[q]==0,1,0)
        }
        # sum_astar = sum(Iasta%*%d_need)
        for (p in 1:(n[cl]-1)){
          sum_astar[p] = Iasta[p,]%*%c(data_num[p,]%*%beta1_mle,data_num[p,]%*%beta2_mle,1)
        }
        
        sum_num[k] = exp(sum(sum_astar))
      }
      Isub = rep(-99,3)
      Isub[1] = ifelse(data_c[sub,4]==1,1,0)
      Isub[2] = ifelse(data_c[sub,4]==2,1,0)
      Isub[3] = ifelse(data_c[sub,4]==0,1,0)
      # num_sub = Isub%*%d_need
      num_sub = Isub%*%c(data_c[sub,1:3]%*%beta1_mle,data_c[sub,1:3]%*%beta2_mle,1)
      num = exp(num_sub)*sum(sum_num)
      ########### calculate denominator
      sum_den = rep(-99,dim(Vbar)[1])
      sum_abar = rep(-99,n[cl])
      for (k in 1:dim(Vbar)[1]){
        abar = Vbar[k,]
        Iabar = matrix(-99,nrow=n[cl],ncol=3)
        for (q in 1:n[cl]){
          Iabar[q,1] = ifelse(abar[q]==1,1,0)
          Iabar[q,2] = ifelse(abar[q]==2,1,0)
          Iabar[q,3] = ifelse(abar[q]==0,1,0)
        }
        # sum_abar = sum(Iabar%*%d_need)
        for (p in 1:n[cl]){
          sum_abar[p] = Iabar[p,]%*%c(data_den[p,]%*%beta1_mle,data_den[p,]%*%beta2_mle,1)
        }
        sum_den[k] = exp(sum(sum_abar))
      }
      prob = num/sum(sum_den)
      prop_ps = rbind(prop_ps,prob)  
    }
  }
  propps = prop_ps[-1]
  data$ps_weight = propps
  data0 = data[data$A==0,]
  data1 = data[data$A==1,]
  data2 = data[data$A==2,]
  ybar_0 = sum(data0$y_obs/data0$ps_weight)/sum(1/data0$ps_weight)
  ybar_1 = sum(data1$y_obs/data1$ps_weight)/sum(1/data1$ps_weight)
  ybar_2 = sum(data2$y_obs/data2$ps_weight)/sum(1/data2$ps_weight)
  effect_suffi = c(ybar_1-ybar_0,ybar_2-ybar_0,ybar_2-ybar_1)
  print(dataid)
  return(effect_suffi = effect_suffi)
}
res_suffi = mclapply(1:1,effectsuffifunc,mc.cores=1)

### result summary
datagene=1
simu = matrix(-99,nrow=datagene,ncol=3)
for (i in 1:datagene){
  simu[i,1] = simudata[[i]]$effec_simu[1]
  simu[i,2] = simudata[[i]]$effec_simu[2]
  simu[i,3] = simudata[[i]]$effec_simu[3]
}
fix = matrix(unlist(res_fix),byrow=T,ncol=3,nrow=datagene)
ran = matrix(unlist(res_ran),byrow=T,ncol=3,nrow=datagene)
suffi = matrix(unlist(res_suffi),byrow=T,ncol=3,nrow=datagene)
naive = matrix(unlist(res_naive),byrow=T,ncol=3,nrow=datagene)

res_summary=as.matrix(rbind(apply(simu,2,mean),apply(naive,2,mean),apply(fix,2,mean),apply(ran,2,mean),apply(suffi,2,mean)))
rownames(res_summary)=c("simu","naive","fix","random","suffi")
colnames(res_summary) = c("effect10","effect20","effect21")
