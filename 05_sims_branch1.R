n_subject <- 200
n_item <- 200
subject <- c(replicate(n_item, c(1:n_subject)))
sd_arr <- c(.3, .3, .3, .3, .3, .3, .3, .3)

tibble(subject) %>%
  group_by(subject) %>%
  mutate(item = seq_along(subject),
         feedback = sample(c("no", "yes"), 1)) %>%
  group_by(subject, item) %>%
  mutate(
    condition = sample(c("restudy","test"), 1)) %>%
  ungroup() %>%
  mutate(p_population = case_when(
    condition == "test" & feedback == "yes" ~ logit_scaled(.5),
    TRUE ~ logit_scaled(.25)
  )) %>%
  group_by(subject) %>%
  mutate(p_sub_Int = case_when(
    condition == "restudy" & feedback == "no" ~ rnorm(1, 0, sd_arr[1]),
    condition == "restudy" & feedback == "yes" ~ rnorm(1, 0, sd_arr[2]),
    TRUE ~ 0
  ),
  p_sub_Cond = case_when(
    condition == "test" & feedback == "no" ~ rnorm(1, 0, sd_arr[3]),
    TRUE ~ 0
  ),
  p_sub_Feedback = case_when(
    condition == "test" & feedback == "yes" ~ rnorm(1, 0, sd_arr[4]),
    TRUE ~ 0
  )) %>%
  group_by(item) %>%
  mutate(p_item_Int = case_when(
    condition == "restudy" & feedback == "no" ~ rnorm(1, 0, sd_arr[5]),
    TRUE ~ 0
  ),
  p_item_Cond = case_when(
    condition == "test" & feedback == "no" ~ rnorm(1, 0, sd_arr[6]),
    TRUE ~ 0
  ),
  p_item_Feedback = case_when(
    condition == "restudy" & feedback == "yes" ~ rnorm(1, 0, sd_arr[7]),
    TRUE ~ 0
  ),
  `p_item_Cond:Feedback` = case_when(
    condition == "test" & feedback == "yes" ~ rnorm(1, 0 , sd_arr[8]),
    TRUE ~ 0
  )) %>% 
  rowwise() %>%
  mutate(correct = rbinom(1, 1, inv_logit_scaled(sum(c(p_population, 
                                                       p_sub_Int,
                                                       p_sub_Cond, 
                                                       p_sub_Feedback, 
                                                       p_item_Int,
                                                       p_item_Cond,
                                                       p_item_Feedback,
                                                       `p_item_Cond:Feedback`))))) %>%
  arrange(subject, item) -> dat_

dat_ %>% ungroup() %>% summarize(across(p_sub_Int:`p_item_Cond:Feedback`, ~sd(.x)))

dat_ <- dat_[sample(nrow(dat_)),]

fit_max_ <- brm(bf_1,
               prior = prior_1,
               sample_prior = TRUE,
               data = dat_,
               threads = threading(floor(ncores/4)))

