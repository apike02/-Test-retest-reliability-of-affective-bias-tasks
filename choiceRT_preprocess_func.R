choiceRT_preprocess_func_testretest <- function(raw_data_t1, raw_data_t2, RTbound = 0.1) {
  #from hbayesdm but modified for testretest
  
  # Use raw_data as a data.frame
  raw_data_t1 <- as.data.frame(raw_data_t1)
  raw_data_t2 <- as.data.frame(raw_data_t2)
  
  # Use general_info of raw_data
  subjs   <- unique(raw_data_t1$subjID)
  n_subj  <- length(subjs)
  
  # Number of upper and lower boundary responses
  Nu <- cbind(with(raw_data_t1, aggregate(choice == 2, by = list(y = subjID), FUN = sum)[["x"]]),
              with(raw_data_t2, aggregate(choice == 2, by = list(y = subjID), FUN = sum)[["x"]]))
  Nl <- cbind(with(raw_data_t1, aggregate(choice == 1, by = list(y = subjID), FUN = sum)[["x"]]),
              with(raw_data_t2, aggregate(choice == 1, by = list(y = subjID), FUN = sum)[["x"]]))
  
  # Reaction times for upper and lower boundary responses
  RTu <- array(-1, c(n_subj, 2, max(Nu)))
  RTl <- array(-1, c(n_subj, 2, max(Nl)))
  
  raw_data<-list(raw_data_t1,raw_data_t2)
  
  for (time in 1:2){
    for (i in 1:n_subj) {
      subj <- subjs[i]
      subj_data <- subset(raw_data[[time]], raw_data[[time]]$subjID == subj)
      
      if (Nu[i] > 0)
        RTu[i, time, 1:Nu[i,time]] <- subj_data$RT[subj_data$choice == 2]  # (Nu/Nl[i]+1):Nu/Nl_max will be padded with 0's
      if (Nl[i] > 0)
        RTl[i, time, 1:Nl[i,time]] <- subj_data$RT[subj_data$choice == 1]  # 0 padding is skipped in likelihood calculation
    }
  }
  
  # Minimum reaction time
  minRT <- cbind(with(raw_data_t1, aggregate(RT, by = list(y = subjID), FUN = min)[["x"]]),
                 with(raw_data_t2, aggregate(RT, by = list(y = subjID), FUN = min)[["x"]]))
  
  # Wrap into a list for Stan
  data_list <- list(
    N       = n_subj,   # Number of subjects
    N_time  = 2,        # Number of timepoints
    Nu_max  = max(Nu),  # Max (across subjects) number of upper boundary responses
    Nl_max  = max(Nl),  # Max (across subjects) number of lower boundary responses
    Nu      = Nu,       # Number of upper boundary responses for each subject
    Nl      = Nl,       # Number of lower boundary responses for each subject
    RTu     = RTu,      # Upper boundary response times
    RTl     = RTl,      # Lower boundary response times
    minRT   = minRT,    # Minimum RT for each subject
    RTbound = RTbound   # Lower bound of RT across all subjects (e.g., 0.1 second)
  )
  
  # Returned data_list will directly be passed to Stan
  return(data_list)
}