## Functions for obtaining transformed effects

transformed_effects <- function(data, model, re_pass = NA){
  if(is.na(re_pass)){
    data$Correct_hat <- rvar(posterior_linpred(model,
                                               re_formula = NA,
                                               resp = "Correct",
                                               transform = TRUE,
                                               newdata = data))
    
    data %>%
      group_by(Feedback) -> data
    
  }else{
    data$Correct_hat <- rvar(posterior_linpred(model,
                                               re_formula = re_pass,
                                               resp = "Correct",
                                               transform = TRUE,
                                               newdata = data))
    
    gr_vars <- c("Cue", "Sub", "Experiment", "TrialGr")
    
    gr_vars <- gr_vars[which(!is.na((str_extract(deparse(re_pass), gr_vars))))]
    
    data %>%
      group_by(Feedback, !!as.symbol(gr_vars)) -> data
  }
  
  if(isTRUE(is.null(model$data$test_effect))){
    data %>%
      reframe(te = Correct_hat[Condition == "Ret"]-
                Correct_hat[Condition == "Study"]) -> data
    
    bind_cols(data, summarize_draws(data$te,
                                    mean,
                                    ~quantile2(.x, probs = c(.025, .975)))) %>%
      select(!c(te, variable)) -> data
  }else{
    bind_cols(data, summarize_draws(data$Correct_hat,
                                    mean,
                                    ~quantile2(.x, probs = c(.025, .975)))) %>%
      summarize(across(mean:q97.5, ~mean(.x))) -> data
  }
  
  return(data)
}

transformed_groups <- function(data, gr_var = NULL){
  
  if(is.null(gr_var)){
    cols <- "Feedback"
  }else{
    cols <- c(gr_var, "Feedback")
  }
  
  data %>%
    distinct(Condition, !!!rlang::syms(cols)) %>%
    mutate(set2 = TRUE,
           set3 = TRUE,
           InitRet_z = 0,
           PropPracCorrect_z = 0) -> pred_dat
  
  data %>%
    group_by(!!!rlang::syms(cols)) %>%
    summarize(test_effect = mean(Correct[Condition == "Ret"] -
                                   mean(Correct[Condition == "Study"]))) %>%
    ungroup() -> dat_temp
  
  if(is.null(gr_var)){
    max_formula <- NULL
  }else if(nrow(distinct(dat_temp, !!as.symbol(gr_var), Feedback)) ==
           nrow(distinct(dat_temp, !!as.symbol(gr_var)))){
    max_formula <- paste0("~(Condition|gr(",
                          gr_var,
                          ", by = Feedback))")
    
    int_formula <- paste0("~(1|gr(",
                          gr_var,
                          ", by = Feedback))")
  }else{
    max_formula <- paste0("~(Condition*Feedback|",
                          gr_var,
                          ")")
    
    int_formula <- paste0("~(1|",
                          gr_var,
                          ")")
  }
  
  # Maximal model  
  dat_temp %>% 
    left_join(transformed_effects(pred_dat, 
                                  fit_Maximal,
                                  ifelse(is.null(gr_var),
                                         NA,
                                         max_formula)),
              by = cols) %>%
    mutate(model = "Maximal") -> temp_Maximal
  
  # Intercept model
  dat_temp %>% 
    left_join(transformed_effects(pred_dat, 
                                  fit_Intercept, 
                                  ifelse(is.null(gr_var),
                                         NA,
                                         int_formula)),
              by = cols) %>%
    mutate(model = "Intercept") -> temp_Intercept
  
  # Aggregate model
  dat_temp %>% 
    left_join(transformed_effects(pred_dat, 
                                  fit_Aggregate),
              by = "Feedback") %>%
    mutate(model = "Aggregate") -> temp_Aggregate
  
  temp_Maximal %>% bind_rows(temp_Intercept, temp_Aggregate) -> dat_temp
  
  if(!is.null(gr_var)){
    dat_temp %>%
      group_by(model) %>%
      mutate(se = (mean-test_effect)^2) %>%
      summarize(p_capture = length(which(test_effect > q2.5 & 
                                           test_effect < q97.5))/
                  length(test_effect),
                mse = mean(se),
                rmse = sqrt(mse)) -> dat_temp
    
    temp_Maximal %>% 
      left_join(temp_Intercept, by = c(gr_var, "Feedback", "test_effect")) %>% 
      left_join(temp_Aggregate, by = c(gr_var, "Feedback", "test_effect")) %>%
      rename("Aggregate" = "mean",
             "Intercept" = "mean.y",
             "Maximal" = "mean.x") -> temp_wide
    
    fit_r <- brm(bf(mvbind(Aggregate, Intercept, Maximal, test_effect)~1,
                    family = student) + set_rescor(TRUE),
                 data = temp_wide,
                 threads = threading(floor(ncores/4)))
    
    as_draws_rvars(fit_r) -> d
    
    summarize_draws(d, mean, ~quantile2(.x, c(.025, .975))) %>%
      filter(str_detect(variable, "rescor") & str_detect(variable, "test")) %>%
      rename("r" = "mean") -> temp_cor
    
    dat_temp %>% bind_cols(temp_cor[,2:4]) -> dat_temp
  }
  
  dat_temp$group <- gr_var
  
  return(distinct(dat_temp))
  # return(dat_temp)
}
