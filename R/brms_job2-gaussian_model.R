# ----------------------------------------------------------------------------
# file: brms_job-gauss_gb_nz.R
# author: Jason Grafmiller
# date: 2022-11-15
# description:
# Fit linear model with raw ratings and z-scroed ratings in brms
# ----------------------------------------------------------------------------

library(here)
require(dplyr)
require(brms)

z. <- function(x){
  (x - mean(x))/(2*sd(x))
}

ratings_gb_nz_df <- here::here("data", "ratings_particle_verbs_2_vars.txt") |>
  read.delim() |> 
  droplevels() |> 
  mutate(z.Corpus_pred = z.(Corpus_pred))

# Use a wider prior for centered raw ratings
priors_gauss1 <- c(
  set_prior("normal(0, 20)", class = "b")
)

message("== fitting raw centered ratings ==\n")
brm_gauss_raw_gb_nz <- brm(
  formula = rating_c ~ (1|id) + (1|item) + (0 + z.Corpus_pred|id) + 
    Variety + z.Corpus_pred + Variety:z.Corpus_pred,
  data = ratings_gb_nz_df,
  prior = priors_gauss1,
  cores = 4,
  iter = 4000,
  save_model = here::here("stan", "brm-gauss1_raw.stan"),
  file = here::here("model_output", "brm-gauss1_raw.rds"),
  file_refit = "on_change",
  refresh = 0
)

if(interactive()){
  pp_check(brm_gauss_raw_gb_nz, ndraws = 100L)
}

message("== fitting scaled ratings ==\n")
# now scale the ratings by participant
ratings_gb_nz_df <- ratings_gb_nz_df |> 
  group_by(id) |> 
  mutate(
    rating_z = scale(rating_raw)[,1]
  ) |> 
  ungroup()

# Since ratings are now on the SD scale, we use a tighter prior
priors_gauss2 <- c(
  set_prior("normal(0, 3)", class = "b")
)

brm_gauss_scaled_gb_nz <- brm(
  formula = rating_z ~ (1|id) + (1|item) + (0 + z.Corpus_pred|id) + 
    Variety + z.Corpus_pred + Variety:z.Corpus_pred,
  data = ratings_gb_nz_df,
  prior = priors_gauss2,
  cores = 4,
  iter = 4000,
  save_model = here::here("stan", "brm-gauss1_scaled.stan"),
  file = here::here("model_output", "brm-gauss2_scaled.rds"),
  file_refit = "on_change",
  refresh = 0
)

if(interactive()){
  pp_check(brm_gauss_scaled_gb_nz, ndraws = 100L)
}



