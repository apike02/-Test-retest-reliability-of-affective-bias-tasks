get_mean_trialwise_likelihood<-function(model,data,parameters,session_parameters,session_task,mean_parameters=0){
  
  if (model=='gng'){
    task_full<-read_tsv(paste0(workingdir,'modelling_gng_',session_task,'.txt'))
  } else if (model=='ezddm'|model=='ddm'){
    task_full<-read_tsv(paste0(workingdir,'modelling_circles_',session_task,'.txt'))
    task_full$subjID<-paste0('P',task_full$subjID)
  }
  likelihood<-list()
  
  if (session_parameters==1&model!='ezddm'){
    suffix=''
  } else {
    suffix=paste0('_',session_parameters)
  }
  
  if  (mean_parameters==1){
    for (participant in unique(data$id)){
      parameters<-data%>%
        select(all_of(paste0(params,suffix)))%>%
        summarise(across(everything(),~mean(.x)))%>%
        rename_with(.cols=everything(),.fn = ~sub("[_]2$", "", .)) #removes _2 from end
      task<-task_full%>%
        filter(subjID==participant)
      if (model=='gng'){
        likelihood[[participant]]<-gng_m5_forward(parameters,task)
      } else if (model=='ezddm'){
        likelihood[[participant]]<-ezddm_forward(parameters,task)
      } else if (model=='ddm'){
        likelihood[[participant]]<-ddm_forward(parameters,task)
      }
    }
  } else if (mean_parameters==0){
  
    for (participant in unique(data$id)){
      parameters<-data%>%
        filter(id==participant)%>%
        select(all_of(paste0(params,suffix)))%>%
        rename_with(.cols=everything(),.fn = ~sub("[_][1-9]$", "", .)) #removes _1 or _2 from end
      task<-task_full%>%
        filter(subjID==participant)
      if (model=='gng'){
        likelihood[[participant]]<-gng_m5_forward(parameters,task)
      } else if (model=='ezddm'){
        likelihood[[participant]]<-ezddm_forward(parameters,task)
      } else if (model=='ddm'){
        likelihood[[participant]]<-ddm_forward(parameters,task)
      }
    }
  }
  
  likelihood<-unlist(lapply(likelihood,mean))
  
  return(likelihood)

}