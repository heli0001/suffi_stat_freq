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
suff_T = rep(-99,m)
for (i in 1:m){
  suff_T[i] = sum(data[data$clu==i,]$A)
}
# which(suff_T==n|suff_T==0)

##### check each cluster, must have treatment and control both ####
while (suff_T[1]==n[1] | suff_T[1]==0){
  sap = apply(as.matrix(pstrue[1:n[1]]),1,rbinom,n=1,size=1)
  data$A[1:n[1]] = sap
  suff_T[1] = sum(sap)
}
begin = 1
total = n[1]
for (i in 2:m){
  begin = begin + n[i-1]
  total = total+ n[i]
  while (suff_T[i]==n[i] | suff_T[i]==0){
    sap = apply(as.matrix(pstrue[begin:total]),1,rbinom,n=1,size=1) 
    data$A[begin:total] = sap
    suff_T[i] = sum(sap)
  }
}
data$sufft = rep(suff_T,n)
data$y0 = 1 + data$x1 + data$x2 + rnorm(N,0,1)
data$y1 = 1 + data$x1 + data$x2 + 2 + rhoyu*data$U + rnorm(N,0,1)
data$y_obs = data$y1*data$A + data$y0*(1-data$A)
effect_true = 2
effec_simu = mean(data$y1-data$y0)
############################## Part 2: PS Estimation ###############################

## 2-1 naive estimator, no weight adjustment
############### 
effec_naiv = mean(data$A*data$y_obs-(1-data$A)*data$y_obs)

## 2-2 random effect model
###############
library(lme4)
control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))
logran = glmer(A ~ x1 + x2+ (1|clu),data=data,family=binomial,control=control)
beta_ran = summary(logran)$coefficients[,1]
m_ran = beta_ran[1]+data$x1*beta_ran[2]+data$x2*beta_ran[3]
w_ran = exp(m_ran)/(1+exp(m_ran))
mean(data$A*data$y_obs/w_ran-(1-data$A)*data$y_obs/(1-w_ran))
data$w_ran = w_ran
effec_ran = mean(data$A*data$y_obs/data$w_ran-(1-data$A)*data$y_obs/(1-data$w_ran))

## 2-3 fixed effect model
###############
logfix = glm(A~x1+x2+factor(clu),family = binomial(link = "logit"),data=data)
data$w_fix = predict(logfix,type="response")
effec_fix = mean(data$A*data$y_obs/data$w_fix-(1-data$A)*data$y_obs/(1-data$w_fix))

## 2-4 new proposed weight
############### 
##### calculate the conditional probability 
library(survival)
clog_summ = clogit(A~x1+x2+strata(clu), data=data)
test = matrix(-99,nrow=10,ncol=2)
for (i in 1:10){
  clog_summ = clogit(A~x1+x2+strata(clu), data=data)
  test[i,] = (clog_summ$coefficients)
}
beta_mle = (clog_summ$coefficients)
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
propps = prop_ps[-1,]
# if (sum(data[data$clu==1,]$A)==0){
#   propps[1:n[1]]=0
# }
# begin = 1
# total = n[1]
# for (i in 2:m){
#   begin = begin + n[i-1]
#   total = total+ n[i]
#   if (sum(data[data$clu==i,]$A)==0){
#     propps[begin:total]=1
#   }
#   else if (sum(data[data$clu==i,]$A)==n[i]){
#     propps[begin:total]=1
#   }
# }
data$ps_weight = propps
effec_new = mean(data$A*data$y_obs/data$ps_weight-(1-data$A)*data$y_obs/(data$ps_weight))















