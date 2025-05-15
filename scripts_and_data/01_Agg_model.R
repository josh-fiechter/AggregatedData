library(tidyverse)
library(here)
library(grafify)
library(brms)
library(posterior)
library(tidybayes)
ncores <- parallel::detectCores()
options(mc.cores = ncores)
set.seed(123)

#########################
## STEP 1: Read in data
#########################

dat <- read_csv(here("testing_effects.csv"))

######################################
## STEP 1: Filter and aggregate data
######################################

# Get rid of Swahili:English pairs from Fiechter and Benjamin (2018)--no variance in `Feedback`
dat %>% 
  filter(!Experiment %in% c("2A", "2B")) -> dat

# Get rid of studies that used a non-24-hr. retention interval
dat %>% 
  filter(!Experiment %in% c("1A", "1C")) -> dat

## Aggregate data
dat %>%
  group_by(Feedback, Experiment, Sub, Condition) %>%
  summarize(Correct = mean(Correct),
            InitRet = mean(InitRet),
            NumPracCorrect = mean(NumPracCorrect),
            NumPracTrial = mean(NumPracTrial)) %>%
  group_by(Feedback, Experiment, Sub) %>%
  summarize(test_effect = Correct[Condition == "Ret"]-Correct[Condition == "Study"],
            InitRet = InitRet[Condition == "Ret"],
            NumPracCorrect = NumPracCorrect[Condition == "Ret"],
            NumPracTrial = NumPracTrial[Condition == "Ret"]) %>%
  ungroup() %>%
  mutate(PropPracCorrect = NumPracCorrect/NumPracTrial,
         InitRet_z = scale(InitRet)[,1],
         PropPracCorrect_z = scale(PropPracCorrect)[,1]) %>%
  ungroup() -> dat_agg

###############################
## STEP 2: Model the data
###############################

bf_1 <- bf(test_effect ~ 0 + Intercept + Feedback*InitRet_z*PropPracCorrect_z,
          family = gaussian)

prior_1 <- prior(cauchy(0,1), class = b)

fit_Aggregate <- brm(bf_1,
              data = dat_agg,
              prior = prior_1,
              sample_prior = TRUE,
              iter = 5000,
              seed = 123,
              save_pars = save_pars(all = TRUE),
              threads = threading(floor(ncores/4)))

# Fit time: 1 s
# No warnings

rm(bf_1, prior_1)

###############################
## STEP 3: Hypothesis tests
###############################

# No feedback
# Feedback

bf_Aggregate <- hypothesis(fit_Aggregate, 
                      c("Intercept = 0", 
                        "Intercept + FeedbackTRUE = 0"))

1/bf_Aggregate$hypothesis$Evid.Ratio

###################################
## STEP 4: Calculate effect sizes
###################################

d_agg <-
  as_draws_df(fit_Aggregate) %>%
  mutate(fb = b_Intercept + b_FeedbackTRUE,
         d_fb = fb/sigma,
         d_nofb = b_Intercept/sigma)

d_agg_fb <- c(mean(d_agg$d_fb),
              sd(d_agg$d_fb),
  quantile(d_agg$d_fb, c(.025, .975)))
  
d_agg_nofb <- c(mean(d_agg$d_nofb),
                sd(d_agg$d_nofb),
  quantile(d_agg$d_nofb, c(.025, .975)))

rm(d_agg)