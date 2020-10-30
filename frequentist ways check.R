rm(list=ls())
load("C:/Users/emmal/Desktop/simdata00m500_small.RData")
load("C:/Users/emmal/Desktop/simdata05m500_small.RData")
load("C:/Users/emmal/Desktop/simdata50m500_small.RData")
load("C:/Users/emmal/Desktop/simdata55m500_small.RData")

dataIndex = 100
library(parallel)
library(lme4)
library(survival)
library(hier.part)

func_naiv = function(dataid,simdata){
  set.seed(dataid)
  data = simdata[[dataid]]$data_obs
  effec_naiv = mean(data$z*data$y_obs-(1-data$z)*data$y_obs)
  return (effec_naiv=effec_naiv)
}
effect_naic_summ00 = unlist(mclapply(1:dataIndex,func_naiv,simdata=simdata00m500_small,mc.cores = 1))
effect_naic_summ05 = unlist(mclapply(1:dataIndex,func_naiv,simdata=simdata05m500_small,mc.cores = 1))
effect_naic_summ50 = unlist(mclapply(1:dataIndex,func_naiv,simdata=simdata50m500_small,mc.cores = 1))
effect_naic_summ55 = unlist(mclapply(1:dataIndex,func_naiv,simdata=simdata55m500_small,mc.cores = 1))
save(effect_naic_summ00,file="effect_naiv00_res_summary.RData")
save(effect_naic_summ05,file="effect_naiv05_res_summary.RData")
save(effect_naic_summ50,file="effect_naiv50_res_summary.RData")
save(effect_naic_summ55,file="effect_naiv55_res_summary.RData")

func_ran = function(simdata,dataid){
  set.seed(dataid)
  data = simdata[[dataid]]$data_obs
  control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))
  logran = glmer(z ~ -1 + x1 + x2+ (1|clu),data=data,family=binomial,control=control)
  beta_ran = summary(logran)$coefficients[,1]
  m_ran = data$x1*beta_ran[1]+data$x2*beta_ran[2]
  w_ran = exp(m_ran)/(1+exp(m_ran))
  #w_ran = predict(logran,type="response")
  data$w_ran = w_ran
  effec_ran = mean(data$z*data$y_obs/data$w_ran-(1-data$z)*data$y_obs/(1-data$w_ran))
  print(dataid)
  return (effec_ran = effec_ran)
}
effect_ran_summ00 = unlist(mclapply(1:dataIndex,func_ran,simdata=simdata00m500_small,mc.cores = 1))
effect_ran_summ05 = unlist(mclapply(1:dataIndex,func_ran,simdata=simdata05m500_small,mc.cores = 1))
effect_ran_summ50 = unlist(mclapply(1:dataIndex,func_ran,simdata=simdata50m500_small,mc.cores = 1))
effect_ran_summ55 = unlist(mclapply(1:dataIndex,func_ran,simdata=simdata55m500_small,mc.cores = 1))
save(effect_ran_summ00,file="effect_ran00_res_summary.RData")
save(effect_ran_summ05,file="effect_ran05_res_summary.RData")
save(effect_ran_summ50,file="effect_ran50_res_summary.RData")
save(effect_ran_summ55,file="effect_ran55_res_summary.RData")

func_fix = function(simdata,dataid){
  set.seed(dataid)
  data = simdata[[dataid]]$data_obs
  logfix = glm(z~-1+x1+x2+factor(clu),family = binomial(link = "logit"),data=data)
  data$w_fix = predict(logfix,type="response")
  effec_fix = mean(data$z*data$y_obs/data$w_fix-(1-data$z)*data$y_obs/(1-data$w_fix))
  print(dataid)
  return (effec_fix=effec_fix)
}
effect_fix_summ00 = unlist(mclapply(1:dataIndex,func_fix,simdata=simdata00m500_small,mc.cores = 1))
effect_fix_summ05 = unlist(mclapply(1:dataIndex,func_fix,simdata=simdata05m500_small,mc.cores = 1))
effect_fix_summ50 = unlist(mclapply(1:dataIndex,func_fix,simdata=simdata50m500_small,mc.cores = 1))
effect_fix_summ55 = unlist(mclapply(1:dataIndex,func_fix,simdata=simdata55m500_small,mc.cores = 1))
save(effect_fix_summ00,file="effect_fix00_res_summary.RData")
save(effect_fix_summ05,file="effect_fix05_res_summary.RData")
save(effect_fix_summ50,file="effect_fix50_res_summary.RData")
save(effect_fix_summ55,file="effect_fix55_res_summary.RData")



func_suffi = function(simdata,dataid){
  set.seed(dataid)
  data = simdata[[dataid]]$data
  data_obs = simdata[[dataid]]$data_obs
  m = simdata[[dataid]]$m
  n = simdata[[dataid]]$n
  suff_T = simdata[[dataid]]$suff_T
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
  print (dataid)
  return (effect_suffi=effect_suffi)
}
effect_suffi_summ00 = unlist(mclapply(1:dataIndex,func_suffi,simdata=simdata00m500_small,mc.cores = 1))
effect_suffi_summ05 = unlist(mclapply(1:dataIndex,func_suffi,simdata=simdata05m500_small,mc.cores = 1))
effect_suffi_summ50 = unlist(mclapply(1:dataIndex,func_suffi,simdata=simdata50m500_small,mc.cores = 1))
effect_suffi_summ55 = unlist(mclapply(1:dataIndex,func_suffi,simdata=simdata55m500_small,mc.cores = 1))
save(effect_suffi_summ00,file="effect_suffi00_res_summary.RData")
save(effect_suffi_summ05,file="effect_suffi05_res_summary.RData")
save(effect_suffi_summ50,file="effect_suffi50_res_summary.RData")
save(effect_suffi_summ55,file="effect_suffi55_res_summary.RData")


################## result 00 ######################
true00 = rep(-99,100)
true05 = rep(-99,100)
true50 = rep(-99,100)
true55 = rep(-99,100)
for (i in 1:100){
  true00[i] = mean(simdata00m500_small[[i]]$data$y1 - simdata00m500_small[[i]]$data$y0)
  true05[i] = mean(simdata05m500_small[[i]]$data$y1 - simdata05m500_small[[i]]$data$y0)
  true50[i] = mean(simdata50m500_small[[i]]$data$y1 - simdata50m500_small[[i]]$data$y0)
  true55[i] = mean(simdata55m500_small[[i]]$data$y1 - simdata55m500_small[[i]]$data$y0)
}
mean(true00)
sd(true00)/10
mean(true05)
sd(true05)/10
mean(true50)
sd(true50)/10
mean(true55)
sd(true55)/10


res_summ00 = as.data.frame(matrix(-99,nrow=4,ncol=4))
res_summ00[1,] = c(mean(effect_naic_summ00),mean(effect_naic_summ00)-mean(true00),sd(effect_naic_summ00)/10,mean((effect_naic_summ00-true00)^2))
res_summ00[2,] = c(mean(effect_fix_summ00),mean(effect_fix_summ00)-mean(true00),sd(effect_fix_summ00)/10,mean((effect_fix_summ00-true00)^2))
res_summ00[3,] = c(mean(effect_ran_summ00),mean(effect_ran_summ00)-mean(true00),sd(effect_ran_summ00)/10,mean((effect_ran_summ00-true00)^2))
res_summ00[4,] = c(mean(effect_suffi_summ00),mean(effect_suffi_summ00)-mean(true00),sd(effect_suffi_summ00)/10,mean((effect_suffi_summ00-true00)^2))
rownames(res_summ00)=c("naive","fix","ran","suffi")
colnames(res_summ00)=c("mean","bias","s.e","mse")

################## result 05 ######################
res_summ05 = as.data.frame(matrix(-99,nrow=4,ncol=4))
res_summ05[1,] = c(mean(effect_naic_summ05),mean(effect_naic_summ05)-mean(true05),sd(effect_naic_summ05)/10,mean((effect_naic_summ05-true05)^2))
res_summ05[2,] = c(mean(effect_fix_summ05),mean(effect_fix_summ05)-mean(true05),sd(effect_fix_summ05)/10,mean((effect_fix_summ05-true05)^2))
res_summ05[3,] = c(mean(effect_ran_summ05),mean(effect_ran_summ05)-mean(true05),sd(effect_ran_summ05)/10,mean((effect_ran_summ05-true05)^2))
res_summ05[4,] = c(mean(effect_suffi_summ05),mean(effect_suffi_summ05)-mean(true05),sd(effect_suffi_summ05)/10,mean((effect_suffi_summ05-true05)^2))
rownames(res_summ05)=c("naive","fix","ran","suffi")
colnames(res_summ05)=c("mean","bias","s.e","mse")

################# result 50 #######################
res_summ50 = as.data.frame(matrix(-99,nrow=4,ncol=4))
res_summ50[1,] = c(mean(effect_naic_summ50),mean(effect_naic_summ50)-mean(true50),sd(effect_naic_summ50)/10,mean((effect_naic_summ50-true50)^2))
res_summ50[2,] = c(mean(effect_fix_summ50),mean(effect_fix_summ50)-mean(true50),sd(effect_fix_summ50)/10,mean((effect_fix_summ50-true50)^2))
res_summ50[3,] = c(mean(effect_ran_summ50),mean(effect_ran_summ50)-mean(true50),sd(effect_ran_summ50)/10,mean((effect_ran_summ50-true50)^2))
res_summ50[4,] = c(mean(effect_suffi_summ50),mean(effect_suffi_summ50)-mean(true50),sd(effect_suffi_summ50)/10,mean((effect_suffi_summ50-true50)^2))
rownames(res_summ50)=c("naive","fix","ran","suffi")
colnames(res_summ50)=c("mean","bias","s.e","mse")



################## result 55 ######################
res_summ55 = as.data.frame(matrix(-99,nrow=4,ncol=4))
res_summ55[1,] = c(mean(effect_naic_summ55),mean(effect_naic_summ55)-mean(true55),sd(effect_naic_summ55)/10,mean((effect_naic_summ55-true55)^2))
res_summ55[2,] = c(mean(effect_fix_summ55),mean(effect_fix_summ55)-mean(true55),sd(effect_fix_summ55)/10,mean((effect_fix_summ55-true55)^2))
res_summ55[3,] = c(mean(effect_ran_summ55),mean(effect_ran_summ55)-mean(true55),sd(effect_ran_summ55)/10,mean((effect_ran_summ55-true55)^2))
res_summ55[4,] = c(mean(effect_suffi_summ55),mean(effect_suffi_summ55)-mean(true55),sd(effect_suffi_summ55)/10,mean((effect_suffi_summ55-true55)^2))
rownames(res_summ55)=c("naive","fix","ran","suffi")
colnames(res_summ55)=c("mean","bias","s.e","mse")
