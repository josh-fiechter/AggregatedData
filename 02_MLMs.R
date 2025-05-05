
################################
## STEP 1: Manipulate the data
################################

dat %>%
  mutate(set1 = ifelse(Condition == "Ret",
                       TRUE,
                       FALSE),
         set2 = ifelse(Condition == "Ret" & NumPracTrial > 1, 
                       TRUE, 
                       FALSE)) -> dat

# Shuffle the data
dat[sample(nrow(dat)),] -> dat

###########################
## STEP 3: Model the data
###########################

## Intercept model

bf_1 <- bf(Correct ~ Condition*Feedback +
             (1|r1|gr(Sub, by = Feedback)) +
             (1|r2|Cue) +
             (1|r3|TrialGr) +
             (1|r4|gr(Experiment, by = Feedback)),
           family = bernoulli)

bf_2 <- bf(InitRet|subset(set1) ~ Feedback +
             (1|r1|gr(Sub, by = Feedback)) +
             (1|r2|Cue) +
             (1|r3|TrialGr) +
             (1|r4|gr(Experiment, by = Feedback)),
           family = bernoulli)

bf_3 <- bf(NumPracCorrect|subset(set2) + trials(NumPracTrial) ~ Feedback +
             (1|r1|gr(Sub, by = Feedback)) +
             (1|r2|Cue) +
             (1|r3|TrialGr) +
             (1|r4|gr(Experiment, by = Feedback)),
           family = binomial)

prior_mv <- prior(cauchy(0,1), class = b, resp = Correct) +
  prior(cauchy(0,1), class = b, resp = InitRet) +
  prior(cauchy(0,1), class = b, resp = NumPracCorrect)

fit_Intercept <- brm(bf_1 + bf_2 + bf_3 + set_rescor(FALSE),
                     data = dat,
                     prior = prior_mv,
                     sample_prior = TRUE,
                     iter = 5000,
                     seed = 123,
                     control = list(adapt_delta = .999),
                     save_pars = save_pars(all = TRUE),
                     threads = threading(floor(ncores/4)))

# Fit time: 4590 s
# 5k transitions exceeded max tree depth

rm(bf_1, bf_2, bf_3, prior_mv)

# Testing effect--no feedback
# Testing effect--feedback

bf_Intercept <- hypothesis(fit_Intercept,
                           c("Correct_ConditionStudy = 0",
                             "Correct_FeedbackTRUE = 
                        Correct_FeedbackTRUE + 
                        Correct_ConditionStudy + 
                        Correct_ConditionStudy:FeedbackTRUE"))

1/bf_Intercept$hypothesis$Evid.Ratio

#########################################################################
#########################################################################

## Maximal model
bf_1 <- bf(Correct ~ Condition*Feedback +
             (Condition|r1|gr(Sub, by = Feedback)) +
             (Condition*Feedback|r2|Cue) +
             (Condition*Feedback|r3|TrialGr) +
             (Condition|r4|gr(Experiment, by = Feedback)),
           family = bernoulli)

bf_2 <- bf(InitRet|subset(set1) ~ Feedback +
             (1|r1|gr(Sub, by = Feedback)) +
             (Feedback|r2|Cue) +
             (Feedback|r3|TrialGr) +
             (1|r4|gr(Experiment, by = Feedback)),
           family = bernoulli)

bf_3 <- bf(NumPracCorrect|subset(set2) + trials(NumPracTrial) ~ Feedback +
             (1|r1|gr(Sub, by = Feedback)) +
             (Feedback|r2|Cue) +
             (Feedback|r3|TrialGr) +
             (1|r4|gr(Experiment, by = Feedback)),
           family = binomial)

prior_mv <- prior(cauchy(0,1), class = b, resp = Correct) +
  prior(cauchy(0,1), class = b, resp = InitRet) +
  prior(cauchy(0,1), class = b, resp = NumPracCorrect)

fit_Maximal <- brm(bf_1 + bf_2 + bf_3 + set_rescor(FALSE),
               data = dat,
               prior = prior_mv,
               sample_prior = TRUE,
               iter = 5000,
               seed = 123,
               control = list(adapt_delta = .99),
               save_pars = save_pars(all = TRUE),
               threads = threading(floor(ncores/4)))

# Fit time: 8098 s
# Max treedepth on 2500 iterations

# This model outperforms the one with only 2 DVs:
#                elpd_diff se_diff
# fit_Maximal_3dv   0.0       0.0  
# fit_Maximal     -38.0       6.9  

rm(bf_1, bf_2, bf_3, prior_mv)

# Testing effect--no feedback
# Testing effect--feedback

bf_Maximal <- hypothesis(fit_Maximal,
                      c("Correct_ConditionStudy = 0",
                        "Correct_FeedbackTRUE = 
                        Correct_FeedbackTRUE + 
                        Correct_ConditionStudy + 
                        Correct_ConditionStudy:FeedbackTRUE"))

1/bf_Maximal$hypothesis$Evid.Ratio

# Differences in `Sub`-level correlations between `InitRet` and `Study` 
# (i.e., testing effect) as a function of feedback
hypothesis(fit_Maximal, 
           c("cor_Sub__Correct_ConditionStudy:FeedbackFALSE__InitRet_Intercept:FeedbackFALSE = 0",
             "cor_Sub__Correct_ConditionStudy:FeedbackTRUE__InitRet_Intercept:FeedbackTRUE = 0",
             "cor_Sub__Correct_ConditionStudy:FeedbackFALSE__InitRet_Intercept:FeedbackFALSE = 
           cor_Sub__Correct_ConditionStudy:FeedbackTRUE__InitRet_Intercept:FeedbackTRUE"), 
           class = NULL)

# Differences in `Sub`-level correlations between `NumPracCorrect` and `Study` 
# (i.e., testing effect) as a function of feedback
hypothesis(fit_Maximal, 
           c("cor_Sub__Correct_ConditionStudy:FeedbackFALSE__NumPracCorrect_Intercept:FeedbackFALSE = 0",
             "cor_Sub__Correct_ConditionStudy:FeedbackTRUE__NumPracCorrect_Intercept:FeedbackTRUE = 0",
             "cor_Sub__Correct_ConditionStudy:FeedbackFALSE__NumPracCorrect_Intercept:FeedbackFALSE = 
           cor_Sub__Correct_ConditionStudy:FeedbackTRUE__NumPracCorrect_Intercept:FeedbackTRUE"), 
           class = NULL)

## Compare the maximal and minimal model
loo_Maximal <- loo(fit_Maximal, resp = "Correct")
loo_Intercept <- loo(fit_Intercept, resp = "Correct")
loo_compare(loo_Maximal, loo_Intercept)

###################################
## STEP 4: Calculate effect sizes
###################################

## Intercept model
d_int <-
  as_draws_df(fit_Intercept, variable = c("^b_", "sd_"), regex = TRUE) %>%
  select(contains("_Correct")) %>%
  mutate(across(contains("sd"), ~.^2)) %>%
  mutate(fb = b_Correct_ConditionStudy + `b_Correct_ConditionStudy:FeedbackTRUE`,
         d_fb = fb/sqrt(rowSums(.[5:10]))*-1,
         d_nofb = b_Correct_ConditionStudy/sqrt(rowSums(.[5:10]))*-1)

d_int_fb <- c(mean(d_int$d_fb),
              sd(d_int$d_fb),
  quantile(d_int$d_fb, c(.025, .975)))

d_int_nofb <- c(mean(d_int$d_nofb),
                sd(d_int$d_nofb),
  quantile(d_int$d_nofb, c(.025, .975)))

rm(d_int)

## Maximal model
d_max <-
  as_draws_df(fit_Maximal, variable = c("^b_", "sd_"), regex = TRUE) %>%
  select(contains("_Correct")) %>%
  mutate(across(contains("sd"), ~.^2)) %>%
  mutate(fb = b_Correct_ConditionStudy + `b_Correct_ConditionStudy:FeedbackTRUE`,
         d_fb = fb/sqrt(rowSums(.[5:20]))*-1,
         d_nofb = b_Correct_ConditionStudy/sqrt(rowSums(.[5:20]))*-1)

d_max_fb <- c(mean(d_max$d_fb),
              sd(d_max$d_fb),
              quantile(d_max$d_fb, c(.025, .975)))

d_max_nofb <- c(mean(d_max$d_nofb),
                sd(d_max$d_nofb),
                quantile(d_max$d_nofb, c(.025, .975)))

rm(d_max)
