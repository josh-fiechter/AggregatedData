
sim_data <- function(n_subject, 
                     n_item, 
                     group_sd, 
                     fit_update = FALSE, 
                     test_fb_correct = .43, 
                     other_correct = .30,
                     data_return = FALSE){
  
  ##################
  ## Generate data
  ##################
  
  subject <- c(replicate(n_item, c(1:n_subject)))
  
  tibble(subject) %>%
    group_by(subject) %>%
    mutate(item = seq_along(subject),
           Feedback = sample(c("no", "yes"), 1)) %>%
    group_by(subject, item) %>%
    mutate(
      Condition = sample(c("Ret","Study"), 1)) %>%
    ungroup() %>%
    mutate(p_population = case_when(
      Condition == "Ret" & Feedback == "yes" ~ logit_scaled(test_fb_correct),
      TRUE ~ logit_scaled(other_correct)
    )) %>%
    group_by(subject) %>%
    mutate(p_sub = case_when(
      Condition == "Ret" & Feedback == "no" ~ rnorm(1, 0 , group_sd),
      Condition == "Study" & Feedback == "no" ~ rnorm(1, 0 , group_sd),
      Condition == "Ret" & Feedback == "yes" ~ rnorm(1, 0 , group_sd),
      Condition == "Study" & Feedback == "yes" ~ rnorm(1, 0 , group_sd)
    )) %>%
    group_by(item) %>%
    mutate(p_item = case_when(
      Condition == "Ret" & Feedback == "no" ~ rnorm(1, 0 , group_sd),
      Condition == "Study" & Feedback == "no" ~ rnorm(1, 0 , group_sd),
      Condition == "Ret" & Feedback == "yes" ~ rnorm(1, 0 , group_sd),
      Condition == "Study" & Feedback == "yes" ~ rnorm(1, 0 , group_sd)
    )) %>% 
    rowwise() %>%
    mutate(Correct = rbinom(1, 1, inv_logit_scaled(sum(c(p_population, p_sub, p_item))))) %>%
    arrange(subject, item) -> dat
  
  dat <- dat[sample(nrow(dat)),]
  
  if(data_return == TRUE){
    return(dat)
  } 
  else{
    ###################
    ## Model the data
    ###################
    
    ## generate aggregate data
    
    dat %>% 
      group_by(subject, Condition, Feedback) %>% 
      summarize(Correct = mean(Correct)) %>%
      group_by(subject, Feedback) %>%
      reframe(test_effect = Correct[Condition == "Ret"]-
                Correct[Condition == "Study"]) -> dat_agg
    
    if(fit_update == FALSE){
      ## Aggregate model ##
      
      bf_1 <- bf(test_effect ~ Feedback,
                 family = gaussian)
      
      prior_1 <- prior(cauchy(0,1), class = b)
      
      fit_Aggregate <- brm(bf_1,
                           prior = prior_1,
                           sample_prior = TRUE,
                           data = dat_agg,
                           threads = threading(floor(ncores/4)))
      
      ## Intercept model ##
      
      bf_1 <- bf(Correct ~ Condition*Feedback +
                   (1|gr(subject, by = Feedback)) +
                   (1|item),
                 family = bernoulli)
      
      prior_1 <- prior(cauchy(0,1), class = b)
      
      fit_Intercept <- brm(bf_1,
                           prior = prior_1,
                           sample_prior = TRUE,
                           data = dat,
                           threads = threading(floor(ncores/4)))
      
      ## Maximal model ##
      
      bf_1 <- bf(Correct ~ Condition*Feedback +
                   (Condition|gr(subject, by = Feedback)) +
                   (Condition*Feedback|item),
                 family = bernoulli)
      
      prior_1 <- prior(cauchy(0,1), class = b)
      
      fit_Maximal <- brm(bf_1,
                         prior = prior_1,
                         sample_prior = TRUE,
                         data = dat,
                         threads = threading(floor(ncores/4)))
    }
    else{
      fit_Aggregate <- update(fit_Aggregate, newdata = dat_agg)
      fit_Intercept <- update(fit_Intercept, newdata = dat)
      fit_Maximal <- update(fit_Maximal, newdata = dat)
    }
    
    return(list(fit_Aggregate = fit_Aggregate,
                fit_Intercept = fit_Intercept,
                fit_Maximal = fit_Maximal)) 
  }
}

########################################################
########################################################

## Run sims

sim_fits <- function(n_sims, 
                     n_subject, 
                     n_item, 
                     group_sd){
  
  dat_sims <- NULL
  
  for (i in 1:n_sims){
    if(exists("fit_Aggregate")){
      temp <- sim_data(n_subject = n_subject, n_item = n_item, group_sd = group_sd, fit_update = TRUE)
    }else{
      temp <- sim_data(n_subject = n_subject, n_item = n_item, group_sd = group_sd, fit_update = FALSE)
    }  
    
    fit_Aggregate <<- temp$fit_Aggregate
    fit_Intercept <<- temp$fit_Intercept
    fit_Maximal <<- temp$fit_Maximal
    
    dat_population <- transformed_groups(fit_Maximal$data)
    dat_population$sim.num <- i
    dat_sims %>% bind_rows(dat_population) -> dat_sims
  }
  
  ## Get average rate of replication
  dat_sims %>% 
    mutate(Aggregate_add = ifelse(model == "Aggregate", -1, 0),
           Intercept_add = ifelse(model == "Intercept", -1, 0),
           Maximal_add = ifelse(model == "Maximal", -1, 0)) %>%
    rowwise() %>%
    mutate(Aggregate_rep = (length(
      which(
        temp$mean[temp$model == "Aggregate" & 
                    temp$Feedback == Feedback] >= q2.5 & 
          temp$mean[temp$model == "Aggregate" & 
                      temp$Feedback == Feedback] <= q97.5)) + Aggregate_add)/(max(temp$sim.num) + Aggregate_add),
      Intercept_rep = (length(
        which(
          temp$mean[temp$model == "Intercept" & 
                      temp$Feedback == Feedback] >= q2.5 & 
            temp$mean[temp$model == "Intercept" & 
                        temp$Feedback == Feedback] <= q97.5)) + Intercept_add)/(max(temp$sim.num) + Intercept_add),
      Maximal_rep = (length(
        which(
          temp$mean[temp$model == "Maximal" & 
                      temp$Feedback == Feedback] >= q2.5 & 
            temp$mean[temp$model == "Maximal" & 
                        temp$Feedback == Feedback] <= q97.5)) + Maximal_add)/(max(temp$sim.num) + Maximal_add),
      n_items = n_item,
      group_sd = group_sd) %>%
    select(Feedback, test_effect, mean, q2.5, q97.5, model, sim.num, 
           Aggregate_rep, Intercept_rep, Maximal_rep, n_items, group_sd)-> dat_sims
  
  return(dat_sims)
}
