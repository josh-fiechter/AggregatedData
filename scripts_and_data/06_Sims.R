library(tidyverse)
library(brms)
library(grafify)
library(posterior)
library(tidybayes)
ncores = parallel::detectCores()
options(mc.cores = ncores)
set.seed(123)

source("03_Transformed_posterior_fxns.R")
source("05_Sims_fxns.R")

dat_12_1 <- sim_fits(100, 120, 12, 2)
dat_12_2 <- sim_fits(100, 120, 12, 1.66)
dat_12_3 <- sim_fits(100, 120, 12, 1.33)
dat_12_4 <- sim_fits(100, 120, 12, 1)
dat_12_5 <- sim_fits(100, 120, 12, .66)
dat_12_6 <- sim_fits(100, 120, 12, .33)

dat_36_1 <- sim_fits(100, 120, 36, 2)
dat_36_2 <- sim_fits(100, 120, 36, 1.66)
dat_36_3 <- sim_fits(100, 120, 36, 1.33)
dat_36_4 <- sim_fits(100, 120, 36, 1)
dat_36_5 <- sim_fits(100, 120, 36, .66)
dat_36_6 <- sim_fits(100, 120, 36, .33)

rm(fit_Aggregate, fit_Intercept, fit_Maximal,
   sim_data, sim_fits, 
   transformed_effects, transformed_groups)

dat_12_1 %>%
  bind_rows(dat_12_2,
            dat_12_3,
            dat_12_4,
            dat_12_5,
            dat_12_6,
            dat_36_1, 
            dat_36_2, 
            dat_36_3,
            dat_36_4,
            dat_36_5,
            dat_36_6) -> dat_sims

rm(dat_12_1, dat_12_2, dat_12_3, dat_12_4, dat_12_5, dat_12_6,
   dat_36_1, dat_36_2, dat_36_3, dat_36_4, dat_36_5, dat_36_6)

## Overview of findings
dat_sims %>% group_by(Feedback, n_items, model) %>%
  summarize(Aggregate_rep = mean(Aggregate_rep),
            Intercept_rep = mean(Intercept_rep),
            Maximal_rep = mean(Maximal_rep))

## Model the findings
bf1 <- bf(mvbind(Aggregate_rep, Intercept_rep, Maximal_rep) ~ Feedback*n_items*model +
            (Feedback*n_items*model|r1|group_sd),
          family = zero_one_inflated_beta)

prior1 <- prior(cauchy(0,1), class = b, resp = Aggregaterep) +
  prior(cauchy(0,1), class = b, resp = Interceptrep) +
  prior(cauchy(0,1), class = b, resp = Maximalrep)

dat_sims <- dat_sims[sample(nrow(dat_sims)),]
dat_sims$n_items <- factor(dat_sims$n_items)

fit <- brm(bf1 + set_rescor(FALSE),
           prior = prior1,
           sample_prior = TRUE,
           data = dat_sims,
           init = 0,
           control = list(adapt_delta = .99),
           seed = 123,
           threads = threading(floor(ncores/4)))

rm(bf1, prior1)

## 12 items--Replication of Aggregate model estimates
hypothesis(fit,
           c("Aggregaterep_modelIntercept = 0",
             "Aggregaterep_modelMaximal = 0",
             "Aggregaterep_modelIntercept = Aggregaterep_modelMaximal",
             "Aggregaterep_modelIntercept + Aggregaterep_Feedbackyes + Aggregaterep_Feedbackyes:modelIntercept = Aggregaterep_Feedbackyes",
             "Aggregaterep_modelMaximal + Aggregaterep_Feedbackyes + Aggregaterep_Feedbackyes:modelMaximal = Aggregaterep_Feedbackyes",
             "Aggregaterep_modelIntercept + Aggregaterep_Feedbackyes + Aggregaterep_Feedbackyes:modelIntercept = Aggregaterep_modelMaximal + Aggregaterep_Feedbackyes + Aggregaterep_Feedbackyes:modelMaximal"))

## 36 items--Replication of Aggregate model estimates
hypothesis(fit,
           c("Aggregaterep_modelIntercept + Aggregaterep_n_items36 + Aggregaterep_n_items36:modelIntercept = 0",
             "Aggregaterep_modelMaximal  + Aggregaterep_n_items36 + Aggregaterep_n_items36:modelMaximal = 0",
             "Aggregaterep_modelIntercept + Aggregaterep_n_items36 + Aggregaterep_n_items36:modelIntercept = Aggregaterep_modelMaximal  + Aggregaterep_n_items36 + Aggregaterep_n_items36:modelMaximal",
             "Aggregaterep_modelIntercept + Aggregaterep_Feedbackyes + Aggregaterep_n_items36 + Aggregaterep_Feedbackyes:n_items36 + Aggregaterep_Feedbackyes:modelIntercept + Aggregaterep_n_items36:modelIntercept + Aggregaterep_Feedbackyes:n_items36:modelIntercept = Aggregaterep_Feedbackyes + Aggregaterep_n_items36 + Aggregaterep_Feedbackyes:n_items36",
             "Aggregaterep_modelMaximal + Aggregaterep_Feedbackyes + Aggregaterep_n_items36 + Aggregaterep_Feedbackyes:n_items36 + Aggregaterep_Feedbackyes:modelMaximal + Aggregaterep_n_items36:modelMaximal + Aggregaterep_Feedbackyes:n_items36:modelMaximal = Aggregaterep_Feedbackyes + Aggregaterep_n_items36 + Aggregaterep_Feedbackyes:n_items36",
             "Aggregaterep_modelIntercept + Aggregaterep_Feedbackyes + Aggregaterep_n_items36 + Aggregaterep_Feedbackyes:n_items36 + Aggregaterep_Feedbackyes:modelIntercept + Aggregaterep_n_items36:modelIntercept + Aggregaterep_Feedbackyes:n_items36:modelIntercept = Aggregaterep_modelMaximal + Aggregaterep_Feedbackyes + Aggregaterep_n_items36 + Aggregaterep_Feedbackyes:n_items36 + Aggregaterep_Feedbackyes:modelMaximal + Aggregaterep_n_items36:modelMaximal + Aggregaterep_Feedbackyes:n_items36:modelMaximal"))

##############################################################################

## 12 items--Replication of Intercept model estimates
hypothesis(fit,
           c("Interceptrep_modelIntercept = 0",
             "Interceptrep_modelMaximal = 0",
             "Interceptrep_modelIntercept = Interceptrep_modelMaximal",
             "Interceptrep_modelIntercept + Interceptrep_Feedbackyes + Interceptrep_Feedbackyes:modelIntercept = Interceptrep_Feedbackyes",
             "Interceptrep_modelMaximal + Interceptrep_Feedbackyes + Interceptrep_Feedbackyes:modelMaximal = Interceptrep_Feedbackyes",
             "Interceptrep_modelIntercept + Interceptrep_Feedbackyes + Interceptrep_Feedbackyes:modelIntercept = Interceptrep_modelMaximal + Interceptrep_Feedbackyes + Interceptrep_Feedbackyes:modelMaximal"))

## 36 items--Replication of Intercept model estimates
hypothesis(fit,
           c("Interceptrep_modelIntercept + Interceptrep_n_items36 + Interceptrep_n_items36:modelIntercept = 0",
             "Interceptrep_modelMaximal  + Interceptrep_n_items36 + Interceptrep_n_items36:modelMaximal = 0",
             "Interceptrep_modelIntercept + Interceptrep_n_items36 + Interceptrep_n_items36:modelIntercept = Interceptrep_modelMaximal  + Interceptrep_n_items36 + Interceptrep_n_items36:modelMaximal",
             "Interceptrep_modelIntercept + Interceptrep_Feedbackyes + Interceptrep_n_items36 + Interceptrep_Feedbackyes:n_items36 + Interceptrep_Feedbackyes:modelIntercept + Interceptrep_n_items36:modelIntercept + Interceptrep_Feedbackyes:n_items36:modelIntercept = Interceptrep_Feedbackyes + Interceptrep_n_items36 + Interceptrep_Feedbackyes:n_items36",
             "Interceptrep_modelMaximal + Interceptrep_Feedbackyes + Interceptrep_n_items36 + Interceptrep_Feedbackyes:n_items36 + Interceptrep_Feedbackyes:modelMaximal + Interceptrep_n_items36:modelMaximal + Interceptrep_Feedbackyes:n_items36:modelMaximal = Interceptrep_Feedbackyes + Interceptrep_n_items36 + Interceptrep_Feedbackyes:n_items36",
             "Interceptrep_modelIntercept + Interceptrep_Feedbackyes + Interceptrep_n_items36 + Interceptrep_Feedbackyes:n_items36 + Interceptrep_Feedbackyes:modelIntercept + Interceptrep_n_items36:modelIntercept + Interceptrep_Feedbackyes:n_items36:modelIntercept = Interceptrep_modelMaximal + Interceptrep_Feedbackyes + Interceptrep_n_items36 + Interceptrep_Feedbackyes:n_items36 + Interceptrep_Feedbackyes:modelMaximal + Interceptrep_n_items36:modelMaximal + Interceptrep_Feedbackyes:n_items36:modelMaximal"))

##############################################################################

## 12 items--Replication of Maximal model estimates
hypothesis(fit,
           c("Maximalrep_modelIntercept = 0",
             "Maximalrep_modelMaximal = 0",
             "Maximalrep_modelIntercept = Maximalrep_modelMaximal",
             "Maximalrep_modelIntercept + Maximalrep_Feedbackyes + Maximalrep_Feedbackyes:modelIntercept = Maximalrep_Feedbackyes",
             "Maximalrep_modelMaximal + Maximalrep_Feedbackyes + Maximalrep_Feedbackyes:modelMaximal = Maximalrep_Feedbackyes",
             "Maximalrep_modelIntercept + Maximalrep_Feedbackyes + Maximalrep_Feedbackyes:modelIntercept = Maximalrep_modelMaximal + Maximalrep_Feedbackyes + Maximalrep_Feedbackyes:modelMaximal"))

## 36 items--Replication of Maximal model estimates
hypothesis(fit,
           c("Maximalrep_modelIntercept + Maximalrep_n_items36 + Maximalrep_n_items36:modelIntercept = 0",
             "Maximalrep_modelMaximal  + Maximalrep_n_items36 + Maximalrep_n_items36:modelMaximal = 0",
             "Maximalrep_modelIntercept + Maximalrep_n_items36 + Maximalrep_n_items36:modelIntercept = Maximalrep_modelMaximal  + Maximalrep_n_items36 + Maximalrep_n_items36:modelMaximal",
             "Maximalrep_modelIntercept + Maximalrep_Feedbackyes + Maximalrep_n_items36 + Maximalrep_Feedbackyes:n_items36 + Maximalrep_Feedbackyes:modelIntercept + Maximalrep_n_items36:modelIntercept + Maximalrep_Feedbackyes:n_items36:modelIntercept = Maximalrep_Feedbackyes + Maximalrep_n_items36 + Maximalrep_Feedbackyes:n_items36",
             "Maximalrep_modelMaximal + Maximalrep_Feedbackyes + Maximalrep_n_items36 + Maximalrep_Feedbackyes:n_items36 + Maximalrep_Feedbackyes:modelMaximal + Maximalrep_n_items36:modelMaximal + Maximalrep_Feedbackyes:n_items36:modelMaximal = Maximalrep_Feedbackyes + Maximalrep_n_items36 + Maximalrep_Feedbackyes:n_items36",
             "Maximalrep_modelIntercept + Maximalrep_Feedbackyes + Maximalrep_n_items36 + Maximalrep_Feedbackyes:n_items36 + Maximalrep_Feedbackyes:modelIntercept + Maximalrep_n_items36:modelIntercept + Maximalrep_Feedbackyes:n_items36:modelIntercept = Maximalrep_modelMaximal + Maximalrep_Feedbackyes + Maximalrep_n_items36 + Maximalrep_Feedbackyes:n_items36 + Maximalrep_Feedbackyes:modelMaximal + Maximalrep_n_items36:modelMaximal + Maximalrep_Feedbackyes:n_items36:modelMaximal"))

##############################################################################

## Visualize fits
dat_sims %>% 
  distinct(Feedback, n_items, model) %>%
  bind_cols(fitted(fit, resp = "Aggregaterep", newdata = ., re_formula = NA)) %>%
  unite("Aggregate_fits", Estimate:Q97.5) %>%
  bind_cols(fitted(fit, resp = "Interceptrep", newdata = ., re_formula = NA)) %>%
  unite("Intercept_fits", Estimate:Q97.5) %>%
  bind_cols(fitted(fit, resp = "Maximalrep", newdata = ., re_formula = NA)) %>%
  unite("Maximal_fits", Estimate:Q97.5) %>%
  pivot_longer(contains("fits"),
               names_to = "Replication model",
               values_to = "Estimates") %>%
  separate_wider_delim(Estimates,
                       delim = "_",
                       names = c("Estimate", "Error", "Q2.5", "Q97.5")) %>%
  mutate(across(Estimate:Q97.5, parse_number),
         Feedback = ifelse(Feedback == "no", 
                           "No Feedback",
                           "Feedback"),
         `Replication model` = str_remove_all(`Replication model`, "_fits")) -> dat_viz

ggplot(dat_viz, aes(x = model, y = Estimate)) +
  facet_wrap(~n_items + Feedback, nrow = 2) +
  geom_pointrange(aes(ymin = Q2.5, ymax = Q97.5,
                      color = `Replication model`,
                      fill = `Replication model`),
                  position = position_dodge(.5),
                  alpha = .5) +
  scale_color_grafify() +
  scale_fill_grafify() +
  ylab("P(Estimate capture)") +
  xlab("Reference model") +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black"),
        panel.grid.major.y = element_line(linetype = "dotted", color = "light gray"))

ggsave("reproducibility_sims.png",
       dpi = 320,
       height = 5,
       width = 6)
