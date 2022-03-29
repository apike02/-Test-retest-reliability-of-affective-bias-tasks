model_selection<-function(fit_filename,ntrials,nparam,llname='loglik',workingdir) {
  name<-load(paste0(workingdir,fit_filename))
  assign('fit_dataframe',get(name)$fit)
  print('Make sure youre using something obtained by sampling!')
  library('loo')
  loo_data<-loo(fit_dataframe,llname)$estimates
  effective_parameters<-loo_data['p_loo','Estimate']
  looic<-loo_data['looic','Estimate']
  log_like<-extract_log_lik(fit_dataframe,llname)
  mean_log_like<-mean(log_like)
  ibic=-2*mean_log_like+effective_parameters*log(ntrials)
  bic=-2*mean_log_like+nparam*log(ntrials)
  return(ic=list(ibic=ibic,bic=bic,looic=looic))
}
