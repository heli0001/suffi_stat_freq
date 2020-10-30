rm(list=ls())
########################Part 1: Data generation #######################
rhoxu = 0
rhoyu = 0
m = 500
n = as.matrix(sample(2:5,m,replace=T,prob=rep(0.25,4)))
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
ps=function(para){
  exp(para[1]+para[2]+para[3])/(1+exp(para[1]+para[2]+para[3]))
}
pstrue = apply(data[,4:6],1,ps)
aij = apply(as.matrix(pstrue),1,rbinom,n=1,size=1)
data$A = aij
data$y0 = 1 + data$x1 + data$x2 + rnorm(N,0,1)
data$y1 = 1 + data$x1 + data$x2 + 2 + rhoyu*data$U + rnorm(N,0,1)
data$y_obs = data$y1*data$A + data$y0*(1-data$A)
effec_true = 2
effec_simu = mean(data$y1-data$y0)

############################## Part 2: PS Estimation ###############################
## 2-1 naive estimator, no weight adjustment
############### 
# effec_naiv = mean(data[data$A==1,9])-mean(data[data$A==0,9])
effec_naiv = mean(data$A*data$y_obs-(1-data$A)*data$y_obs)
######### random model #########
library(Matrix)
library(lme4)
library(nlme)
logran = glmer(A ~ x1 + x2+ (1|clu),data=data,family=binomial,control = glmerControl(optimizer = "bobyqa"),
               nAGQ = 10)
data$w_ran = predict(logran)
effec_ran = mean(data$A*data$y_obs/data$w_ran-(1-data$A)*data$y_obs/(1-data$w_ran))

######### fixed model ##########
#library(fastDummies)
#data_new = dummy_cols(data, select_columns = "clu")
#bbb = glm(as.factor(A) ~ .,family = "binomial",data=data_new[,c(3,4,6,11:510)])
logfix = glm(factor(A)~1 + x1+x2+factor(clu),family = binomial(link = "logit"),data=data)
#summary(logfix)$coef
data$w_fix = predict(logfix,type="response")

effec_fix = mean(data$A*data$y_obs/data$w_fix-(1-data$A)*data$y_obs/(1-data$w_fix))
effec_fix_baye1 = mean(data$A*data$y_obs/data$ps-(1-data$A)*data$y_obs/(1-data$ps))
effec_fix_baye2 = mean(data$A*data$y_obs/data$ps2-(1-data$A)*data$y_obs/(1-data$ps2))
effec_ran_baye = mean(data$A*data$y_obs/data$psran-(1-data$A)*data$y_obs/(1-data$psran))
######### new proposed weight

##### calculate the conditional probability 
data = data[,1:9]
suff_T = rep(-99,m)
for (i in 1:m){
  suff_T[i] = sum(data[data$clu==i,]$A)
}
##### obtain MLE
data$sufft = rep(suff_T,n)
library(survival)
clog_summ = clogit(A~x1+x2+sufft+strata(clu), data=data)
beta_mle = exp(clog_summ$coefficients)[1:2]

###########
# if multinomial, could use mclogit
#library(mclogit)
#clog_summ = mclogit(cbind(clu, A)~x1+x2, data=data)
#beta_mle = exp(clog_summ$coefficients)


####### construct Vij
beta_est=beta_mle
prop_ps = 0
library(hier.part)
for (i in 1:m){
  for (j in 1:n[i]){
    astar = rbind(rep(0,n[i]),combos(n[i])$binary) ## all combinations of treatment
    astarj = astar[astar[,j]==data[data$clu==i,]$A[j],] ## jth component is equal to the observed
    sumrow = rowSums(astarj)
    ind = which(sumrow==suff_T[i])
    astar_f = matrix(astarj[ind,],ncol=n[i]) ## and only keep rowsum equal to sufficient statistics 
    nume = exp(data[data$clu==i,]$A[j]*as.matrix(data[data$clu==i,4:5][j,])%*%beta_est)*sum(exp(astar_f[,-j]%*%(as.matrix(data[data$clu==i,4:5])[-j,]%*%beta_est)))
    
    ind2 = which(rowSums(astar)==suff_T[i]) # find rows sum=sufficient statistics
    abar = matrix(astar[ind2,],ncol=n[i]) # only keep those
    denomi = sum(exp(abar%*%as.matrix(data[data$clu==i,4:5])%*%beta_est))
    
    prob = nume/denomi
    prop_ps = rbind(prop_ps,prob)
  }
}
propps = prop_ps[-1]
if (sum(data[data$clu==1,]$A)==0){
  propps[1:n[1]]=0
}
begin = 1
total = n[1]
for (i in 2:m){
  begin = begin + n[i-1]
  total = total+ n[i]
  if (sum(data[data$clu==i,]$A)==0){
    propps[begin:total]=0
  }
  else if (sum(data[data$clu==i,]$A)==n[i]){
    propps[begin:total]=1
  }
}
data$ps_weight = propps
effec_new = mean(data[data$A==1,]$y_obs/data[data$A==1,]$ps_weight) - mean(data[data$A==0,]$y_obs/(data[data$A==0,]$ps_weight))
mean(data$A*data$y_obs/data$ps_weight-(1-data$A)*data$y_obs/(data$ps_weight))














