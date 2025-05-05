
###############################
## Population-level estimates
###############################

dat_population <- transformed_groups(dat)

############################################
## Compare # of captured group-level means
############################################

## Subject-level estimates
dat_subs <- transformed_groups(dat, "Sub")

## Item-level estimates
dat_items <- transformed_groups(dat, "Cue")

## Trial-level estimates
dat_trials <- transformed_groups(dat, "TrialGr")

## Experiment-level estimates
dat_experiments <- transformed_groups(dat, "Experiment")

## Put group-level estimates together
dat_items %>%
  bind_rows(dat_subs,
            dat_trials,
            dat_experiments) -> dat_groups

rm(dat_items, dat_subs, dat_trials, dat_experiments)

########################################
## Compare estimates of testing effect
########################################

dat_population %>%
  rename("Estimate" = mean, "CI.Lower" = q2.5, "CI.Upper" = q97.5) %>%
  mutate(Feedback = ifelse(Feedback == TRUE, "Feedback", "No Feedback")) %>%
  mutate(N = case_when(
    model == "Aggregate" & Feedback == "Feedback" ~ 313,
    model == "Aggregate" & Feedback == "No Feedback" ~ 312,
    model != "Aggregate" & Feedback == "Feedback" ~ 3512,
    model != "Aggregate" & Feedback == "No Feedback" ~ 3520
  )) %>%
  arrange(Feedback, model) -> dat_population

ggplot(dat_population, aes(y = Estimate, 
                  ymin = CI.Lower, 
                  ymax = CI.Upper, 
                  x = model)) +
  geom_hline(yintercept = 0, 
             color = "red",
             linetype = "dotted") +
  geom_errorbar(aes(color = model),
                width = .3) +
  geom_point(aes(color = model,
                 fill = model,
                 size = N),
             alpha = .5) +
  geom_point(aes(color = model),
             shape = 4) +
  scale_color_grafify() +
  facet_wrap(~Feedback, 
             nrow = 1) +
  ylab("Testing effect") +
  xlab("Model") +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, 
                                    color = "black"),
        panel.grid.major.y = element_line(linetype = "dotted", 
                                          color = "light gray"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1))

ggsave("test_effects.png",
       dpi = 320,
       height = 3,
       width = 3.5)

##################################################
## Compare estimates of standardized effect size
##################################################

dat_d <- bind_rows(as_tibble(t(cbind(d_agg_fb, d_int_fb, d_max_fb))),
                   as_tibble(t(cbind(d_agg_nofb, d_int_nofb, d_max_nofb)))) %>%
  rename("Estimate" = V1, "SD" = V2, "CI.Lower" = `2.5%`, "CI.Upper" = `97.5%`) %>%
  mutate(Model = c(replicate(2, c("Aggregate", "Intercept", "Maximal"))),
         Feedback = c(replicate(3, "Feedback"), 
                      replicate(3, "No Feedback"))) %>%
  mutate(N = case_when(
    Model == "Aggregate" & Feedback == "Feedback" ~ 313,
    Model == "Aggregate" & Feedback == "No Feedback" ~ 312,
    Model != "Aggregate" & Feedback == "Feedback" ~ 3512,
    Model != "Aggregate" & Feedback == "No Feedback" ~ 3520
  )) %>%
  arrange(Feedback)

rm(d_agg_fb, d_agg_nofb, d_int_fb, d_int_nofb, d_max_fb, d_max_nofb)

##########################
## Compare Bayes factors
##########################

# These numbers would be best presented in a table owing to the different scales
dat_bf <- as_tibble(c(1/bf_Aggregate$hypothesis$Evid.Ratio,
                      1/bf_Intercept$hypothesis$Evid.Ratio,
                      1/bf_Maximal$hypothesis$Evid.Ratio)) %>%
  mutate(Model = c(replicate(2, "Aggregate"), 
                   replicate(2, "Intercept"), 
                   replicate(2, "Maximal")),
         Feedback = c(replicate(3, c("No Feedback", "Feedback")))) %>%
  arrange(Feedback)

dat_bf
