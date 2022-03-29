ddm_forward<-function(parameters,task){
  
  library('RWiener') #for pwiener
  
  likelihood=rep(0,nrow(task))

  for (t in 1:nrow(task)){
    probA<-pwiener(10, alpha=parameters$alpha, tau=parameters$tau, beta=parameters$beta, delta=parameters$delta,resp='lower')
    
    if (task$choice[t]==1){
      likelihood[t]=probA
    } else {
      likelihood[t]=1-probA
    }
    
  } # end of t loop
  return(likelihood)
}