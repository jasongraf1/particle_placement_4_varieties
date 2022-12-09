# ----------------------------------------------------------------------------
# file: brms_job-zoib_gb_nz.R
# author: Jason Grafmiller
# date: 2022-11-15
# description:
# Fit ZOIB model with brms
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

priors_gb <- c(
  set_prior("normal(0, 1.5)", class = "b"), # our mu parameter
  set_prior("normal(0, 2)", class = "b", dpar = "phi"),
  set_prior("normal(0, 1.5)", class = "b", dpar = "zoi"),
  set_prior("normal(0, 1.5)", class = "b", dpar = "coi")
)

model_zoib_gb_nz1 <- bf(
  rating_prop ~ (1|id) + (1|item) + (0 + z.Corpus_pred|id) +
    Variety + z.Corpus_pred + Variety:z.Corpus_pred,
  phi ~ (1|id) + (1|item) + (0 + z.Corpus_pred|id) +
    Variety + z.Corpus_pred + Variety:z.Corpus_pred,
  zoi ~ (1|id) + (1|item) + (0 + z.Corpus_pred|id) +
    Variety + z.Corpus_pred + Variety:z.Corpus_pred,
  coi ~ (1|id) + (1|item) + (0 + z.Corpus_pred|id) +
    Variety + z.Corpus_pred + Variety:z.Corpus_pred,
  family = zero_one_inflated_beta()
)

brm_zoib_gb_nz1 <- brm(
  formula = model_zoib_gb_nz1,
  data = ratings_gb_nz_df,
  prior = priors_gb,
  cores = parallel::detectCores(),
  iter = 4000,
  file = here::here("model_output", "brm-zoib_gb_nz1.rds"),
  file_refit = "on_change",
  refresh = 0
)

if(interactive()){
  pp_check(brm_zoib_gb_nz1, ndraws = 100L)
}

