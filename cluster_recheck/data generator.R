library(parallel)
rhoxu = 0
rhoyu = 0
m = 500

data_gene = function(index,rhoxu=0,rhoyu=0,m=500,effec_true = 2){
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
  data$y1 = 1 + data$x1 + data$x2 + effec_true + rhoyu*data$U + rnorm(N,0,1)
  data$y_obs = data$y1*data$A + data$y0*(1-data$A)
  return (list(n=n,suff_T=suff_T,data=data))
}
simdata00 = mclapply(1:10,data_gene,rhoxu=0,rhoyu=0,m=500,effec_true = 2,mc.cores=1)
simdata05 = mclapply(1:10,data_gene,rhoxu=0,rhoyu=5,m=500,effec_true = 2,mc.cores=1)
simdata50 = mclapply(1:10,data_gene,rhoxu=5,rhoyu=0,m=500,effec_true = 2,mc.cores=1)
simdata55 = mclapply(1:10,data_gene,rhoxu=5,rhoyu=5,m=500,effec_true = 2,mc.cores=1)
save(simdata00,file="1000 simulated00 data.RData")
save(simdata05,file="1000 simulated05 data.RData")
save(simdata50,file="1000 simulated50 data.RData")
save(simdata55,file="1000 simulated55 data.RData")