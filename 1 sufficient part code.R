################### efficient way to calculate binary case ##############
SumAllComb <- function(m1, n, x1) {
  # x1 is a vector 
  # n = length(x1) 
  if (m1==0) {
    return(1)
  }else if (m1>n){
    return(0)
  }else if (m1<0) {
    return(0)
  }else if (n==0 & m1>0) {
    return(0)
  }else {
    return(SumAllComb(m1,n-1,x1) + SumAllComb(m1-1,n-1,x1)*x1[n])
  }
}

icpwfunc <- function(exb1, a1, ni){
  # calculate icpw in a culster
  # exb1 is a vector of (exp(X_i1%*%beta1),...,exp(X_ini%*%beta1)) w.r.t. the beta1 in category 1
  # ni is sample size of cluster i
  # a1 is vector of (A_i1,..,A_ini)
  t1 <- sum(a1)
  if(t1==ni | ni==1 | t1==0){
    return(rep(1,ni))
  }else{
    denom <- SumAllComb(m1=t1, n=ni, x1=exb1)
    numer <- sapply(c(1:ni), FUN=function(j){SumAllComb(x1=exb1[-j],m1=(t1-a1[j]),n=(ni-1))})
    numer[a1==1] = numer[a1==1]*exb1[a1==1]
    ps_suffi_cl = numer/denom
    return(ps_suffi_cl)
  }
}

icpw_all = function(mind,data,n){
  data_cl1 = data[data$clu==mind,]
  ps_suffi_cl1 = icpwfunc(data_cl1$exb1,data_cl1$a1,n[mind])
  print(mind)
  return(ps_suffi_cl1)
}

ps_est = function(x,para){
  1/(1+exp(-x[1:2]%*%para-x[3]))
}


iter=11000
burn=1000

effect_func = function(dataid,simdata,clu_res,iter=iter,burn=burn){
  ptm = proc.time()
  set.seed(dataid)
  data_obs = simdata[[dataid]]$data_obs
  m = simdata[[dataid]]$m
  n = simdata[[dataid]]$n
  gamma_pos = clu_res[[dataid]]$gamma_pos
  delta = clu_res[[dataid]]$delta
  zeta = clu_res[[dataid]]$zeta
  
  ## dummy variable/indicator matrix for cluster
  cluind = data_obs$clu
  Istar = matrix(0, nrow=length(cluind), ncol=length(unique(cluind)))
  Istar[cbind(seq_along(cluind), cluind)] = 1
  
  
  ###### based on NMIG for significant clusters, use sufficient method to calculate pscore, 
  ###### others still use the logistic one
  gamma_est = apply(gamma_pos[(burn+1):iter,],2,mean)
  exb1 = exp(as.matrix(data_obs[,3:4])%*%gamma_est)
  a1 = data_obs$z
  data_obs$exb1 = exb1
  data_obs$a1 = a1
  
  # nonsig clusters
  inclu_prob = apply(delta[(burn+1):iter,],2,mean)
  non_inclu = which(inclu_prob<0.5)
  data_ps = cbind(data_obs,Istar)
  data_nsig = data_ps[data_ps$clu %in% non_inclu,]
  zeta_est =  apply(zeta[(burn+1):iter,],2,mean)[non_inclu]
  data_need = cbind(data_nsig$x1,data_nsig$x2,as.matrix(data_nsig[,-c(1:9)][,non_inclu])%*%zeta_est)
  data_nsig$ps_est = apply(data_need,1,ps_est,para=gamma_est)
  
  # sig clusters
  inclu = which(inclu_prob>=0.5)
  data_sig = data_obs[data_obs$clu %in% inclu,]
  data_sig$ps_est = unlist(lapply(inclu,icpw_all,data_obs,n))
  
  #################### calculate using IPW ###########################
  data_new = rbind(data_nsig[,c(1:9,m+10)],data_sig)
  effect_nmig_suffi = mean(data_new$z*data_new$y_obs/data_new$ps_est-(1-data_new$z)*data_new$y_obs/(1-data_new$ps_est))

  print(dataid)
  time_span = proc.time() - ptm
  return(list(effect_nmig_suffi=effect_nmig_suffi,time_span=time_span))
}
effect_res = mclapply(1,effect_func,simdatam20,clu_res,iter=iter,burn=burn,mc.cores = 1)



