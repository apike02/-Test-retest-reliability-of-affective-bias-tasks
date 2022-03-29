gng_m5_forward<-function(parameters,task){
  
  library('boot') #for inv.logit
  
  sv=rep(0,4)
  wv_g=rep(0,4)
  wv_ng=rep(0,4)
  qv_g=rep(0,4)
  qv_ng=rep(0,4)
  sv=rep(0,4)
  pGo=rep(0,4)
  
  likelihood=rep(0,nrow(task))
  
  for (t in 1:nrow(task)){
  
    if (sv[task$cue[t]]>0){
      pavlovian=parameters$app;
    } else {
      pavlovian=parameters$av;
    }
    
    wv_g[task$cue[t]]  = qv_g[task$cue[t]] + parameters$gobias + pavlovian * sv[task$cue[t]];
    wv_ng[task$cue[t]] = qv_ng[task$cue[t]];  # qv_ng is always equal to wv_ng (regardless of action)
    pGo[task$cue[t]]   = inv.logit(wv_g[task$cue[t]] - wv_ng[task$cue[t]])
    
    # noise
    pGo[task$cue[t]]   = pGo[task$cue[t]] * (1 - parameters$xi)
    pGo[task$cue[t]]   = pGo[task$cue[t]] + parameters$xi/2
    
    if (task$keyPressed[t]==1){
      likelihood[t]=pGo[task$cue[t]]
    } else {
      likelihood[t]=1-pGo[task$cue[t]]
    }
    
    # after receiving feedback, update sv[t + 1]
    if (task$outcome[t] >= 0) {
      sv[task$cue[t]] = sv[task$cue[t]] + parameters$ep * (parameters$rhoRew * task$outcome[t] - sv[task$cue[t]])
    } else {
      sv[task$cue[t]] = sv[task$cue[t]] + parameters$ep * (parameters$rhoPun * task$outcome[t] - sv[task$cue[t]])
    }
    
    # update action values
    if (task$keyPressed[t]==1) { # update go value
      if (task$outcome[t] >=0) {
        qv_g[task$cue[t]] = qv_g[task$cue[t]] + parameters$ep * (parameters$rhoRew * task$outcome[t] - qv_g[task$cue[t]]);
      } else {
        qv_g[task$cue[t]] = qv_g[task$cue[t]] + parameters$ep * (parameters$rhoPun * task$outcome[t] - qv_g[task$cue[t]]);
      }
    } else { # update no-go value
      if (task$outcome[t] >=0) {
        qv_ng[task$cue[t]] = qv_ng[task$cue[t]] + parameters$ep * (parameters$rhoRew * task$outcome[t] - qv_ng[task$cue[t]]);
      } else {
        qv_ng[task$cue[t]] = qv_ng[task$cue[t]] + parameters$ep * (parameters$rhoPun * task$outcome[t] - qv_ng[task$cue[t]]);
      }
    }
  } # end of t loop
  return(likelihood)
}