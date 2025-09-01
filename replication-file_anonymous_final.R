
# DATA AND CODE MATERIALS FOR MANUSCRIPT: 
#   Comparing Residual- and Observation-Level Cross-Lagged Panel Models with Fixed Effects
#   Interpretive and Empirical Consequences


# Setup -------------------------------------------------------------------

# Packages
library(dplyr)
library(tibble)
library(lavaan)

# Clear working directory 
rm(list = ls())

# Set seed for replicability
set.seed(20240129)

# Calculate average time-varying coefficient 
mean_summary_func <- function(model) {
  
  phihatbar <- parTable(model) %>% 
    filter(substr(label, 1, 3) == "phi") %>% 
    select(est) %>% 
    lapply(., mean)
  
  gamhatbar <- parTable(model) %>% 
    filter(substr(label, 1, 3) == "gam") %>% 
    select(est) %>% 
    lapply(., mean)
  
  bethatbar <- parTable(model) %>% 
    filter(substr(label, 1, 3) == "bet") %>% 
    select(est) %>% 
    lapply(., mean)
  
  rhohatbar <- parTable(model) %>% 
    filter(substr(label, 1, 3) == "rho") %>% 
    select(est) %>% 
    lapply(., mean)
  
  return(list(phihatbar = phihatbar, 
              gamhatbar = gamhatbar, 
              bethatbar = bethatbar, 
              rhohatbar = rhohatbar))
  
}


# Generate data -----------------------------------------------------------

# --- Time-invariant DGP 

# For scenarios 1a-2b and 5a-6b

sim_func_ti <- function(n, 
                        sd_ux, sd_uy, 
                        sd_ex, sd_ey,
                        phi, gam, 
                        bet, rho,
                        obs = TRUE) {
  
  ux <- rnorm(n, 0, sd_ux)
  uy <- rnorm(n, 0, sd_uy)
  
  if(obs == TRUE) {
    
    df <- tibble(
      x0 = 1 * ux + rnorm(n, 0, sd_ex),
      y0 = 1 * uy + rnorm(n, 0, sd_ey),
      # x0 = rnorm(n, 0, sd_ex),
      # y0 = rnorm(n, 0, sd_ey),
      x1 = phi * x0 + gam * y0 + 1 * ux + rnorm(n, 0, sd_ex),
      y1 = bet * x0 + rho * y0 + 1 * uy + rnorm(n, 0, sd_ey),
      x2 = phi * x1 + gam * y1 + 1 * ux + rnorm(n, 0, sd_ex),
      y2 = bet * x1 + rho * y1 + 1 * uy + rnorm(n, 0, sd_ey),
      x3 = phi * x2 + gam * y2 + 1 * ux + rnorm(n, 0, sd_ex),
      y3 = bet * x2 + rho * y2 + 1 * uy + rnorm(n, 0, sd_ey),
      x4 = phi * x3 + gam * y3 + 1 * ux + rnorm(n, 0, sd_ex),
      y4 = bet * x3 + rho * y3 + 1 * uy + rnorm(n, 0, sd_ey),
      x5 = phi * x4 + gam * y4 + 1 * ux + rnorm(n, 0, sd_ex),
      y5 = bet * x4 + rho * y4 + 1 * uy + rnorm(n, 0, sd_ey)
    )
    
  } else {
    
    df <- tibble(
      x0 = 1 * ux + rnorm(n, 0, sd_ex),
      y0 = 1 * uy + rnorm(n, 0, sd_ey),
      # x0 = rnorm(n, 0, sd_ex),
      # y0 = rnorm(n, 0, sd_ey),
      x1 = phi * x0 + gam * y0 + (1 - phi) * ux - gam * uy + rnorm(n, 0, sd_ex),
      y1 = bet * x0 + rho * y0 - bet * ux + (1 - rho) * uy + rnorm(n, 0, sd_ey), 
      x2 = phi * x1 + gam * y1 + (1 - phi) * ux - gam * uy + rnorm(n, 0, sd_ex),
      y2 = bet * x1 + rho * y1 - bet * ux + (1 - rho) * uy + rnorm(n, 0, sd_ey), 
      x3 = phi * x2 + gam * y2 + (1 - phi) * ux - gam * uy + rnorm(n, 0, sd_ex),
      y3 = bet * x2 + rho * y2 - bet * ux + (1 - rho) * uy + rnorm(n, 0, sd_ey), 
      x4 = phi * x3 + gam * y3 + (1 - phi) * ux - gam * uy + rnorm(n, 0, sd_ex),
      y4 = bet * x3 + rho * y3 - bet * ux + (1 - rho) * uy + rnorm(n, 0, sd_ey), 
      x5 = phi * x4 + gam * y4 + (1 - phi) * ux - gam * uy + rnorm(n, 0, sd_ex),
      y5 = bet * x4 + rho * y4 - bet * ux + (1 - rho) * uy + rnorm(n, 0, sd_ey)
    )
    
  }
  
  return(df)
  
}

df_obs_ti <- sim_func_ti(n = 1000000, 
                         sd_ux = 1, sd_uy = 1, 
                         sd_ex = 1, sd_ey = 1, 
                         phi = 0.2, gam = 0.5,
                         bet = 0.5, rho = 0.2, 
                         obs = TRUE)
head(df_obs_ti)

df_res_ti <- sim_func_ti(n = 1000000, 
                         sd_ux = 1, sd_uy = 1, 
                         sd_ex = 1, sd_ey = 1, 
                         phi = 0.2, gam = 0.5,
                         bet = 0.5, rho = 0.2, 
                         obs = FALSE)
head(df_res_ti)


# --- Time-varying DGP 

# For scenarios 3a-4b and 7a-8b

sim_func_ftv <- function(n, 
                         sd_ux, sd_uy, 
                         sd_ex, sd_ey,
                         phi, gam, 
                         bet, rho,
                         obs = TRUE) {
  
  ux <- rnorm(n, 0, sd_ux)
  uy <- rnorm(n, 0, sd_uy)
  
  phi <- rnorm(n = 5, mean = phi, sd = 0.01)
  gam <- rnorm(n = 5, mean = gam, sd = 0.01)
  bet <- rnorm(n = 5, mean = bet, sd = 0.01)
  rho <- rnorm(n = 5, mean = rho, sd = 0.01)
  
  lamx0 <- rnorm(1, mean = 1, sd = 0.01)
  lamy0 <- rnorm(1, mean = 1, sd = 0.01)  
  
  lamx <- rnorm(n = 5, mean = 1, sd = 0.01)
  lamy <- rnorm(n = 5, mean = 1, sd = 0.01)
  
  if(obs == TRUE) {
    
    df <- tibble(
      x0 = lamx0 * ux + rnorm(n, 0, sd_ex),
      y0 = lamy0 * uy + rnorm(n, 0, sd_ey),
      # x0 = rnorm(n, 0, sd_ex),
      # y0 = rnorm(n, 0, sd_ey),
      x1 = phi[1] * x0 + gam[1] * y0 + lamx[1] * ux + rnorm(n, 0, sd_ex),
      y1 = bet[1] * x0 + rho[1] * y0 + lamy[1] * uy + rnorm(n, 0, sd_ey),
      x2 = phi[2] * x1 + gam[2] * y1 + lamx[2] * ux + rnorm(n, 0, sd_ex),
      y2 = bet[2] * x1 + rho[2] * y1 + lamy[2] * uy + rnorm(n, 0, sd_ey),
      x3 = phi[3] * x2 + gam[3] * y2 + lamx[3] * ux + rnorm(n, 0, sd_ex),
      y3 = bet[3] * x2 + rho[3] * y2 + lamy[3] * uy + rnorm(n, 0, sd_ey),
      x4 = phi[4] * x3 + gam[4] * y3 + lamx[4] * ux + rnorm(n, 0, sd_ex),
      y4 = bet[4] * x3 + rho[4] * y3 + lamy[4] * uy + rnorm(n, 0, sd_ey),
      x5 = phi[5] * x4 + gam[5] * y4 + lamx[5] * ux + rnorm(n, 0, sd_ex),
      y5 = bet[5] * x4 + rho[5] * y4 + lamy[5] * uy + rnorm(n, 0, sd_ey)
    )
    
  } else {
    
    df <- tibble(
      x0 = lamx0 * ux + rnorm(n, 0, sd_ex),
      y0 = lamy0 * uy + rnorm(n, 0, sd_ey),
      # x0 = rnorm(n, 0, sd_ex),
      # y0 = rnorm(n, 0, sd_ey),
      x1 = phi[1] * x0 + gam[1] * y0 + (lamx[1] - phi[1] * lamx0) * ux - (gam[1] * lamy0) * uy + rnorm(n, 0, sd_ex),
      y1 = bet[1] * x0 + rho[1] * y0 - (bet[1] * lamx0) * ux + (lamy[1] - rho[1] * lamy0) * uy + rnorm(n, 0, sd_ey), 
      x2 = phi[2] * x1 + gam[2] * y1 + (lamx[2] - phi[2] * lamx[1]) * ux - (gam[2] * lamy[1]) * uy + rnorm(n, 0, sd_ex),
      y2 = bet[2] * x1 + rho[2] * y1 - (bet[2] * lamx[1]) * ux + (lamy[2] - rho[2] * lamy[1]) * uy + rnorm(n, 0, sd_ey), 
      x3 = phi[3] * x2 + gam[3] * y2 + (lamx[3] - phi[3] * lamx[2]) * ux - (gam[3] * lamy[2]) * uy + rnorm(n, 0, sd_ex),
      y3 = bet[3] * x2 + rho[3] * y2 - (bet[3] * lamx[2]) * ux + (lamy[3] - rho[3] * lamy[2]) * uy + rnorm(n, 0, sd_ey), 
      x4 = phi[4] * x3 + gam[4] * y3 + (lamx[4] - phi[4] * lamx[3]) * ux - (gam[4] * lamy[3]) * uy + rnorm(n, 0, sd_ex),
      y4 = bet[4] * x3 + rho[4] * y3 - (bet[4] * lamx[3]) * ux + (lamy[4] - rho[4] * lamy[3]) * uy + rnorm(n, 0, sd_ey), 
      x5 = phi[5] * x4 + gam[5] * y4 + (lamx[5] - phi[5] * lamx[4]) * ux - (gam[5] * lamy[4]) * uy + rnorm(n, 0, sd_ex),
      y5 = bet[5] * x4 + rho[5] * y4 - (bet[5] * lamx[4]) * ux + (lamy[5] - rho[5] * lamy[4]) * uy + rnorm(n, 0, sd_ey)
    )
    
  }
  
  return(df)
  
}

df_obs_tv <- sim_func_ftv(n = 1000000, 
                           sd_ux = 1, sd_uy = 1, 
                           sd_ex = 1, sd_ey = 1, 
                           phi = 0.2, gam = 0.5,
                           bet = 0.5, rho = 0.2, 
                           obs = TRUE)
head(df_obs_tv)

df_obs_tv %>% 
  summarise(across(starts_with("x"), mean))
df_obs_tv %>% 
  summarise(across(starts_with("y"), mean))

df_res_tv <- sim_func_ftv(n = 1000000, 
                           sd_ux = 1, sd_uy = 1, 
                           sd_ex = 1, sd_ey = 1, 
                           phi = 0.2, gam = 0.5,
                           bet = 0.5, rho = 0.2, 
                           obs = FALSE)
head(df_res_tv)

df_res_tv %>% 
  summarise(across(starts_with("x"), mean))
df_res_tv %>% 
  summarise(across(starts_with("y"), mean))


# Simulation analysis  ----------------------------------------------------

{
  # --- Time invariant models 
  
  # For scenarios oddly-numbered scenarios
  
  mdp_ti <- "
  # Individual effects
  ux =~ 1*x2 + 1*x3 + 1*x4 + 1*x5
  uy =~ 1*y2 + 1*y3 + 1*y4 + 1*y5
  # Regressions (time invariant)
  x2 ~ phi*x1 + gam*y1
  y2 ~ bet*x1 + rho*y1
  x3 ~ phi*x2 + gam*y2
  y3 ~ bet*x2 + rho*y2
  x4 ~ phi*x3 + gam*y3
  y4 ~ bet*x3 + rho*y3
  x5 ~ phi*x4 + gam*y4
  y5 ~ bet*x4 + rho*y4
  # Correlations 
  ux ~~ uy + x1 + y1
  uy ~~ x1 + y1
  x1 ~~ y1
  # Residual variances constant over time 
  x2 ~~ ex*x2
  x3 ~~ ex*x3
  x4 ~~ ex*x4
  x5 ~~ ex*x5
  y2 ~~ ey*y2
  y3 ~~ ey*y3
  y4 ~~ ey*y4
  y5 ~~ ey*y5
"
  
  mri_ti <- "
  # Individual effects
  ux =~ 1*x1 + 1*x2 + 1*x3 + 1*x4 + 1*x5
  uy =~ 1*y1 + 1*y2 + 1*y3 + 1*y4 + 1*y5
  # Identify residuals 
  x1w =~ 1*x1
  x2w =~ 1*x2
  x3w =~ 1*x3
  x4w =~ 1*x4
  x5w =~ 1*x5
  y1w =~ 1*y1
  y2w =~ 1*y2
  y3w =~ 1*y3
  y4w =~ 1*y4
  y5w =~ 1*y5
  # Constrain original variances to zero 
  y1 ~~ 0*y1
  y2 ~~ 0*y2
  y3 ~~ 0*y3
  y4 ~~ 0*y4
  y5 ~~ 0*y5
  x1 ~~ 0*x1
  x2 ~~ 0*x2
  x3 ~~ 0*x3
  x4 ~~ 0*x4
  x5 ~~ 0*x5
  # Regressions (time invariant)
  x2w ~ phi*x1w + gam*y1w 
  y2w ~ bet*x1w + rho*y1w
  x3w ~ phi*x2w + gam*y2w 
  y3w ~ bet*x2w + rho*y2w
  x4w ~ phi*x3w + gam*y3w 
  y4w ~ bet*x3w + rho*y3w
  x5w ~ phi*x4w + gam*y4w 
  y5w ~ bet*x4w + rho*y4w
  # Correlations   
  ux ~~ uy + 0*x1w + 0*y1w 
  uy ~~ 0*x1w + 0*y1w
  x1w ~~ y1w
  # Residual variances constant over time
  x2w ~~ ex*x2w
  x3w ~~ ex*x3w
  x4w ~~ ex*x4w
  x5w ~~ ex*x5w
  y2w ~~ ey*y2w
  y3w ~~ ey*y3w
  y4w ~~ ey*y4w
  y5w ~~ ey*y5w
  # Need to constrain to zero (bug?)
  x5w ~~ 0*y5w
"
  
  # --- Time-varying models 
  
  # For evenly-numbered scenarios 
  
  mdp_tv <- "
  # Individual effects
  ux =~ x2 + x3 + x4 + x5
  uy =~ y2 + y3 + y4 + y5
  # Regressions (time invariant)
  x2 ~ phi2*x1 + gam2*y1
  y2 ~ bet2*x1 + rho2*y1
  x3 ~ phi3*x2 + gam3*y2
  y3 ~ bet3*x2 + rho3*y2
  x4 ~ phi4*x3 + gam4*y3
  y4 ~ bet4*x3 + rho4*y3
  x5 ~ phi5*x4 + gam5*y4
  y5 ~ bet5*x4 + rho5*y4
  # Correlations 
  ux ~~ uy + x1 + y1
  uy ~~ x1 + y1
  x1 ~~ y1
  # Residual variances constant over time 
  x2 ~~ ex*x2
  x3 ~~ ex*x3
  x4 ~~ ex*x4
  x5 ~~ ex*x5
  y2 ~~ ey*y2
  y3 ~~ ey*y3
  y4 ~~ ey*y4
  y5 ~~ ey*y5
"
  
  mri_tv <- "
  # Individual effects
  ux =~ x1 + x2 + x3 + x4 + x5
  uy =~ y1 + y2 + y3 + y4 + y5
  # Identify residuals 
  x1w =~ 1*x1
  x2w =~ 1*x2
  x3w =~ 1*x3
  x4w =~ 1*x4
  x5w =~ 1*x5
  y1w =~ 1*y1
  y2w =~ 1*y2
  y3w =~ 1*y3
  y4w =~ 1*y4
  y5w =~ 1*y5
  # Constrain original variances to zero 
  y1 ~~ 0*y1
  y2 ~~ 0*y2
  y3 ~~ 0*y3
  y4 ~~ 0*y4
  y5 ~~ 0*y5
  x1 ~~ 0*x1
  x2 ~~ 0*x2
  x3 ~~ 0*x3
  x4 ~~ 0*x4
  x5 ~~ 0*x5
  # Regressions (time invariant)
  x2w ~ phi2*x1w + gam2*y1w 
  y2w ~ bet2*x1w + rho2*y1w
  x3w ~ phi3*x2w + gam3*y2w 
  y3w ~ bet3*x2w + rho3*y2w
  x4w ~ phi4*x3w + gam4*y3w 
  y4w ~ bet4*x3w + rho4*y3w
  x5w ~ phi5*x4w + gam5*y4w 
  y5w ~ bet5*x4w + rho5*y4w
  # Correlations   
  ux ~~ uy + 0*x1w + 0*y1w 
  uy ~~ 0*x1w + 0*y1w
  x1w ~~ y1w
  # Residual variances constant over time
  x2w ~~ ex*x2w
  x3w ~~ ex*x3w
  x4w ~~ ex*x4w
  x5w ~~ ex*x5w
  y2w ~~ ey*y2w
  y3w ~~ ey*y3w
  y4w ~~ ey*y4w
  y5w ~~ ey*y5w
  # Need to constrain to zero (bug?)
  x5w ~~ 0*y5w
"
  
}


# Simulation  -------------------------------------------------------------

# --- Run models 

# 1a
mdp_ti_obs_ti.fit <- sem(mdp_ti, data = df_obs_ti, estimator = "ML")
summary(mdp_ti_obs_ti.fit)

# 1b 
mri_ti_obs_ti.fit <- sem(mri_ti, data = df_obs_ti, estimator = "ML")
summary(mri_ti_obs_ti.fit)

# 2a 
mdp_tv_obs_ti.fit <- sem(mdp_tv, data = df_obs_ti, estimator = "ML")
summary(mdp_tv_obs_ti.fit)

# 2b 
mri_tv_obs_ti.fit <- sem(mri_tv, data = df_obs_ti, estimator = "ML")
summary(mri_tv_obs_ti.fit)



# 3a
mdp_ti_obs_tv.fit <- sem(mdp_ti, data = df_obs_tv, estimator = "ML")
summary(mdp_ti_obs_tv.fit)

# 3b 
mri_ti_obs_tv.fit <- sem(mri_ti, data = df_obs_tv, estimator = "ML")
summary(mri_ti_obs_tv.fit)

# 4a 
mdp_tv_obs_tv.fit <- sem(mdp_tv, data = df_obs_tv, estimator = "ML")
summary(mdp_tv_obs_tv.fit)

# 4b 
mri_tv_obs_tv.fit <- sem(mri_tv, data = df_obs_tv, estimator = "ML")
summary(mri_tv_obs_tv.fit)



# 5a
mdp_ti_res_ti.fit <- sem(mdp_ti, data = df_res_ti, estimator = "ML")
summary(mdp_ti_res_ti.fit)

# 5b 
mri_ti_res_ti.fit <- sem(mri_ti, data = df_res_ti, estimator = "ML")
summary(mri_ti_res_ti.fit)

anova(mri_ti_res_ti.fit, mdp_ti_res_ti.fit)

# 6a 
mdp_tv_res_ti.fit <- sem(mdp_tv, data = df_res_ti, estimator = "ML")
summary(mdp_tv_res_ti.fit)

# 6b 
mri_tv_res_ti.fit <- sem(mri_tv, data = df_res_ti, estimator = "ML")
summary(mri_tv_res_ti.fit)



# 7a
mdp_ti_res_tv.fit <- sem(mdp_ti, data = df_res_tv, estimator = "ML")
summary(mdp_ti_res_tv.fit)

# 7b 
mri_ti_res_tv.fit <- sem(mri_ti, data = df_res_tv, estimator = "ML")
summary(mri_ti_res_tv.fit)

anova(mri_ti_res_tv.fit, mdp_ti_res_tv.fit)

# 8a 
mdp_tv_res_tv.fit <- sem(mdp_tv, data = df_res_tv, estimator = "ML")
summary(mdp_tv_res_tv.fit)

# 8b 
mri_tv_res_tv.fit <- sem(mri_tv, data = df_res_tv, estimator = "ML")
summary(mri_tv_res_tv.fit)

anova(mri_tv_res_tv.fit, mdp_tv_res_tv.fit)




# --- Overview 

# | Nr. | DGP | DGP-Type | Dataframe | Model   | Model-Type | Model name  | 
# |-----|-----|----------|-----------|---------|------------|-------------|
# | 1a  | OBS | TI       | df_obs_ti | DPM     | TI         | mdp_ti      |
# | 1b  | OBS | TI       | df_obs_ti | RI-CLPM | TI         | mri_ti      |
# | 2a  | OBS | TI       | df_obs_ti | DPM     | TV         | mdp_tv      |
# | 2b  | OBS | TI       | df_obs_ti | RI-CLPM | TV         | mri_tv      |
# |-----|-----|----------|-----------|---------|------------|-------------|
# | 3a  | OBS | TV       | df_obs_tv | DPM     | TI         | mdp_ti      |
# | 3b  | OBS | TV       | df_obs_tv | RI-CLPM | TI         | mri_ti      |
# | 4a  | OBS | TV       | df_obs_tv | DPM     | TV         | mdp_tv      |
# | 4b  | OBS | TV       | df_obs_tv | RI-CLPM | TV         | mri_tv      |
# |-----|-----|----------|-----------|---------|------------|-------------|
# | 5a  | RES | TI       | df_res_ti | DPM     | TI         | mdp_ti      |
# | 5b  | RES | TI       | df_res_ti | RI-CLPM | TI         | mri_ti      |
# | 6a  | RES | TI       | df_res_ti | TV      | TV         | mdp_tv      |
# | 6b  | RES | TI       | df_res_ti | TV      | TV         | mri_tv      |
# |-----|-----|----------|-----------|---------|------------|-------------|
# | 7a  | RES | TV       | df_res_tv | DPM     | TI         | mdp_ti      |
# | 7b  | REs | TV       | df_res_tv | RI-CLPM | TI         | mri_ti      |
# | 8a  | RES | TV       | df_res_tv | DPM     | TV         | mdp_tv      | 
# | 8b  | RES | TV       | df_res_tv | RI-CLPM | TV         | mri_tv      |


# --- Pull mean estimates over time for Table 2 

models <- paste0(rep(paste0(rep(c("mdp", "mri"), 2), "_", rep(c("ti", "tv"), each = 2)), 4),
                 "_",
                 c(rep("obs_ti", 4), rep("obs_tv", 4), rep("res_ti", 4), rep("res_tv", 4)),
                 ".fit")

params <- c("phi", "gam", "bet", "rho")

# Function: takes parameter prefix and a lavaan object
new_func <- function(param, model) {
  
  parameterestimates(model) %>% 
    filter(substr(label, 1, 3) == param) %>% 
    summarise(mean_est = mean(est, na.rm = TRUE)) %>% 
    pull(mean_est)
  }

# Loop over combinations
pm_df <- expand.grid(param = params, model = models, stringsAsFactors = FALSE)

results <- mapply(function(param, modelname) {
  
  model_obj <- get(modelname)  
  
  new_func(param, model_obj)
  }, 
  pm_df$param, pm_df$model)

# Turn into matrix
res_matrix <- matrix(results, nrow = length(params), ncol = length(models), byrow = FALSE)
rownames(res_matrix) <- params
colnames(res_matrix) <- models

round(res_matrix, 3)

# --- Get chisq

new_func2 <- function(model) {
  
  fm <- fitmeasures(model)
  
  c(chisq = round(fm["chisq"], 3),
    df = fm["df"],
    pvalue = round(fm["pvalue"], 3))
  }

results <- lapply(models, function(modelname) {
  
  model_obj <- get(modelname)
  
  new_func2(model_obj)
  }
)

fit_indices <- do.call(rbind, results)
rownames(fit_indices) <- models
fit_indices <- as.data.frame(fit_indices)

fit_indices


# Longer DGP for 4b -------------------------------------------------------

# For scenarios 3a-4b and 7a-8b

sim_func_ftv <- function(n, 
                         sd_ux, sd_uy, 
                         sd_ex, sd_ey,
                         phi, gam, 
                         bet, rho,
                         obs = TRUE) {
  
  ux <- rnorm(n, 0, sd_ux)
  uy <- rnorm(n, 0, sd_uy)
  
  phi <- rnorm(n = 10, mean = phi, sd = 0.01)
  gam <- rnorm(n = 10, mean = gam, sd = 0.01)
  bet <- rnorm(n = 10, mean = bet, sd = 0.01)
  rho <- rnorm(n = 10, mean = rho, sd = 0.01)
  
  lamx0 <- rnorm(1, mean = 1, sd = 0.01)
  lamy0 <- rnorm(1, mean = 1, sd = 0.01)  
  
  lamx <- rnorm(n = 10, mean = 1, sd = 0.01)
  lamy <- rnorm(n = 10, mean = 1, sd = 0.01)
  
  if(obs == TRUE) {
    
    df <- tibble(
      x0 = lamx0 * ux + rnorm(n, 0, sd_ex),
      y0 = lamy0 * uy + rnorm(n, 0, sd_ey),
      # x0 = rnorm(n, 0, sd_ex),
      # y0 = rnorm(n, 0, sd_ey),
      x1 = phi[1] * x0 + gam[1] * y0 + lamx[1] * ux + rnorm(n, 0, sd_ex),
      y1 = bet[1] * x0 + rho[1] * y0 + lamy[1] * uy + rnorm(n, 0, sd_ey),
      x2 = phi[2] * x1 + gam[2] * y1 + lamx[2] * ux + rnorm(n, 0, sd_ex),
      y2 = bet[2] * x1 + rho[2] * y1 + lamy[2] * uy + rnorm(n, 0, sd_ey),
      x3 = phi[3] * x2 + gam[3] * y2 + lamx[3] * ux + rnorm(n, 0, sd_ex),
      y3 = bet[3] * x2 + rho[3] * y2 + lamy[3] * uy + rnorm(n, 0, sd_ey),
      x4 = phi[4] * x3 + gam[4] * y3 + lamx[4] * ux + rnorm(n, 0, sd_ex),
      y4 = bet[4] * x3 + rho[4] * y3 + lamy[4] * uy + rnorm(n, 0, sd_ey),
      x5 = phi[5] * x4 + gam[5] * y4 + lamx[5] * ux + rnorm(n, 0, sd_ex),
      y5 = bet[5] * x4 + rho[5] * y4 + lamy[5] * uy + rnorm(n, 0, sd_ey),
      x6 = phi[6] * x5 + gam[6] * y5 + lamx[6] * ux + rnorm(n, 0, sd_ex),
      y6 = bet[6] * x5 + rho[6] * y5 + lamy[6] * uy + rnorm(n, 0, sd_ey),
      x7 = phi[7] * x6 + gam[7] * y6 + lamx[7] * ux + rnorm(n, 0, sd_ex),
      y7 = bet[7] * x6 + rho[7] * y6 + lamy[7] * uy + rnorm(n, 0, sd_ey),
      x8 = phi[8] * x7 + gam[8] * y7 + lamx[8] * ux + rnorm(n, 0, sd_ex),
      y8 = bet[8] * x7 + rho[8] * y7 + lamy[8] * uy + rnorm(n, 0, sd_ey),
      x9 = phi[9] * x8 + gam[9] * y8 + lamx[9] * ux + rnorm(n, 0, sd_ex),
      y9 = bet[9] * x8 + rho[9] * y8 + lamy[9] * uy + rnorm(n, 0, sd_ey)
    )
    
  } else {
    
    df <- tibble(
      x0 = lamx0 * ux + rnorm(n, 0, sd_ex),
      y0 = lamy0 * uy + rnorm(n, 0, sd_ey),
      # x0 = rnorm(n, 0, sd_ex),
      # y0 = rnorm(n, 0, sd_ey),
      x1 = phi[1] * x0 + gam[1] * y0 + (lamx[1] - phi[1] * lamx0) * ux - (gam[1] * lamy0) * uy + rnorm(n, 0, sd_ex),
      y1 = bet[1] * x0 + rho[1] * y0 - (bet[1] * lamx0) * ux + (lamy[1] - rho[1] * lamy0) * uy + rnorm(n, 0, sd_ey), 
      x2 = phi[2] * x1 + gam[2] * y1 + (lamx[2] - phi[2] * lamx[1]) * ux - (gam[2] * lamy[1]) * uy + rnorm(n, 0, sd_ex),
      y2 = bet[2] * x1 + rho[2] * y1 - (bet[2] * lamx[1]) * ux + (lamy[2] - rho[2] * lamy[1]) * uy + rnorm(n, 0, sd_ey), 
      x3 = phi[3] * x2 + gam[3] * y2 + (lamx[3] - phi[3] * lamx[2]) * ux - (gam[3] * lamy[2]) * uy + rnorm(n, 0, sd_ex),
      y3 = bet[3] * x2 + rho[3] * y2 - (bet[3] * lamx[2]) * ux + (lamy[3] - rho[3] * lamy[2]) * uy + rnorm(n, 0, sd_ey), 
      x4 = phi[4] * x3 + gam[4] * y3 + (lamx[4] - phi[4] * lamx[3]) * ux - (gam[4] * lamy[3]) * uy + rnorm(n, 0, sd_ex),
      y4 = bet[4] * x3 + rho[4] * y3 - (bet[4] * lamx[3]) * ux + (lamy[4] - rho[4] * lamy[3]) * uy + rnorm(n, 0, sd_ey), 
      x5 = phi[5] * x4 + gam[5] * y4 + (lamx[5] - phi[5] * lamx[4]) * ux - (gam[5] * lamy[4]) * uy + rnorm(n, 0, sd_ex),
      y5 = bet[5] * x4 + rho[5] * y4 - (bet[5] * lamx[4]) * ux + (lamy[5] - rho[5] * lamy[4]) * uy + rnorm(n, 0, sd_ey), 
      x6 = phi[6] * x5 + gam[6] * y5 + (lamx[6] - phi[6] * lamx[5]) * ux - (gam[6] * lamy[5]) * uy + rnorm(n, 0, sd_ex),
      y6 = bet[6] * x5 + rho[6] * y5 - (bet[6] * lamx[5]) * ux + (lamy[6] - rho[6] * lamy[5]) * uy + rnorm(n, 0, sd_ey),
      x7 = phi[7] * x6 + gam[7] * y6 + (lamx[7] - phi[7] * lamx[6]) * ux - (gam[7] * lamy[6]) * uy + rnorm(n, 0, sd_ex),
      y7 = bet[7] * x6 + rho[7] * y6 - (bet[7] * lamx[6]) * ux + (lamy[7] - rho[7] * lamy[6]) * uy + rnorm(n, 0, sd_ey),
      x8 = phi[8] * x7 + gam[8] * y7 + (lamx[8] - phi[8] * lamx[7]) * ux - (gam[8] * lamy[7]) * uy + rnorm(n, 0, sd_ex),
      y8 = bet[8] * x7 + rho[8] * y7 - (bet[8] * lamx[7]) * ux + (lamy[8] - rho[8] * lamy[7]) * uy + rnorm(n, 0, sd_ey),
      x9 = phi[9] * x8 + gam[9] * y8 + (lamx[9] - phi[9] * lamx[8]) * ux - (gam[9] * lamy[8]) * uy + rnorm(n, 0, sd_ex),
      y9 = bet[9] * x8 + rho[9] * y8 - (bet[9] * lamx[8]) * ux + (lamy[9] - rho[9] * lamy[8]) * uy + rnorm(n, 0, sd_ey)
    )
    
  }
  
  return(df)
  
}

df_obs_tv_late <- sim_func_ftv(n = 1000000, 
                          sd_ux = 1, sd_uy = 1, 
                          sd_ex = 1, sd_ey = 1, 
                          phi = 0.2, gam = 0.5,
                          bet = 0.5, rho = 0.2, 
                          obs = TRUE)
head(df_obs_tv_late)

df_obs_tv %>% 
  summarise(across(starts_with("x"), mean))
df_obs_tv %>% 
  summarise(across(starts_with("y"), mean))

mri_tv <- "
  # Individual effects
  ux =~ x5 + x6 + x7 + x8 + x9
  uy =~ y5 + y6 + y7 + y8 + y9
  # Identify residuals 
  x5w =~ 1*x5
  x6w =~ 1*x6
  x7w =~ 1*x7
  x8w =~ 1*x8
  x9w =~ 1*x9
  y5w =~ 1*y5
  y6w =~ 1*y6
  y7w =~ 1*y7
  y8w =~ 1*y8
  y9w =~ 1*y9
  # Constrain original variances to zero 
  y5 ~~ 0*y5
  y6 ~~ 0*y6
  y7 ~~ 0*y7
  y8 ~~ 0*y8
  y9 ~~ 0*y9
  x5 ~~ 0*x5
  x6 ~~ 0*x6
  x7 ~~ 0*x7
  x8 ~~ 0*x8
  x9 ~~ 0*x9
  # Regressions (time invariant)
  x6w ~ phi6*x5w + gam6*y5w 
  y6w ~ bet6*x5w + rho6*y5w
  x7w ~ phi7*x6w + gam7*y6w 
  y7w ~ bet7*x6w + rho7*y6w
  x8w ~ phi8*x7w + gam8*y7w 
  y8w ~ bet8*x7w + rho8*y7w
  x9w ~ phi9*x8w + gam9*y8w 
  y9w ~ bet9*x8w + rho9*y8w
  # Correlations   
  ux ~~ uy + 0*x5w + 0*y5w 
  uy ~~ 0*x5w + 0*y5w
  x5w ~~ y5w
  # Residual variances constant over time
  x6w ~~ ex*x6w
  x7w ~~ ex*x7w
  x8w ~~ ex*x8w
  x9w ~~ ex*x9w
  y6w ~~ ey*y6w
  y7w ~~ ey*y7w
  y8w ~~ ey*y8w
  y9w ~~ ey*y9w
  # Need to constrain to zero (bug?)
  x9w ~~ 0*y9w
"

# 4bstar 
mri_tv_obs_tv_star.fit <- sem(mri_tv, data = df_obs_tv_late, estimator = "ML")
summary(mri_tv_obs_tv_star.fit)

# --- Average estimates over time

# Phi 
parameterestimates(mri_tv_obs_tv_star.fit) %>% 
  filter(substr(label, 1, 1) == "p") %>% 
  summarise(mean_est = mean(est))
# = .207

# Gamma 
parameterestimates(mri_tv_obs_tv_star.fit) %>% 
  filter(substr(label, 1, 1) == "g") %>% 
  summarise(mean_est = mean(est))
# = .502

# Beta
parameterestimates(mri_tv_obs_tv_star.fit) %>% 
  filter(substr(label, 1, 1) == "b") %>% 
  summarise(mean_est = mean(est))
# = .505

# Rho 
parameterestimates(mri_tv_obs_tv_star.fit) %>% 
  filter(substr(label, 1, 1) == "r") %>% 
  summarise(mean_est = mean(est))
# = .202


# Make Figure 8 -----------------------------------------------------------

# Packages
library(ggplot2)
library(reshape2)

# Make long format dataframe of results
dfl <- tibble(melt(res_matrix)); dfl

# Some recoding for labelling purposes
dfl <- dfl %>% 
  mutate(dgp = substr(Var2, 8, 10), 
         dgp_tv = substr(Var2, 12, 13), 
         mod = substr(Var2, 1, 3),
         mod_tv = substr(Var2, 5, 6), 
         dgp = ifelse(dgp == "obs", "OBS", "RES"),
         dgp_tv = ifelse(dgp_tv == "ti", "TI", "TV"),
         mod = ifelse(mod == "mdp", "DPM", "RI-CLPM"),
         mod_tv = ifelse(mod_tv == "ti", "TI", "TV"), 
         Var1 = ifelse(Var1 == "phi", "Phi", 
                       ifelse(Var1 == "gam", "Gamma", 
                              ifelse(Var1 == "bet", "Beta", 
                                     ifelse(Var1 == "rho", "Rho", NA)))), 
         dgp_full = paste0(dgp, " ", dgp_tv),
         collab = paste0("Model: ", mod_tv)); dfl

# Express in deviations from known true (mean) parameter values
dfl$dev <- ifelse(dfl$Var1 %in% c("Phi", "Rho"), dfl$value - .2, dfl$value - .5); dfl

# Plot
dfl %>% 
  ggplot(aes(x = factor(dgp_full, levels = rev(c("OBS TI", "RES TI", "OBS TV", "RES TV"))), y = dev, color = mod, shape = mod)) + 
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") + 
  scale_y_continuous(name = "Deviation from true parameter value") + 
  scale_x_discrete(name = "DGP") + 
  scale_color_discrete(name = "Model") + 
  scale_shape_discrete(name = "Model") + 
  coord_flip() + 
  facet_grid(rows = vars(Var1), cols = vars(collab))


# Make Figure 3 -----------------------------------------------------------

# Packages
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(grid)
library(ggpubr)

# Clear working directory
rm(list = ls())

# Set a seed
set.seed(20240124)

# Just five observations for sake of visualization
n <- 5L

# Set the autoregressive paramter
rho <- 0.6

# Initial realization is fixed
eta <- rnorm(n, mean = 0, sd = 1.0)
nu0 <- -2:2

dfw <- tibble(id = 1:n,
              eta = eta, 
              y0 = nu0, 
              y1 = eta + rho * y0 + rnorm(n, mean = 0, sd = 0.1),
              y2 = eta + rho * y1 + rnorm(n, mean = 0, sd = 0.1),
              y3 = eta + rho * y2 + rnorm(n, mean = 0, sd = 0.1),
              y4 = eta + rho * y3 + rnorm(n, mean = 0, sd = 0.1),
              y5 = eta + rho * y4 + rnorm(n, mean = 0, sd = 0.1),
              y6 = eta + rho * y5 + rnorm(n, mean = 0, sd = 0.1),
              y7 = eta + rho * y6 + rnorm(n, mean = 0, sd = 0.1),
              y8 = eta + rho * y7 + rnorm(n, mean = 0, sd = 0.1),
              y9 = eta + rho * y8 + rnorm(n, mean = 0, sd = 0.1))

dfl <- tidyr::gather(dfw, key = time, value = y, y0:y9, factor_key = FALSE)
dfl$time <- as.numeric(substr(dfl$time, 2, 2))

p1 <- dfl %>% 
  ggplot(aes(x = time, y = y, color = as.factor(id))) + 
  geom_point() + 
  geom_line() + 
  geom_line(aes(y = (1 - rho)^-1 * eta), linetype = "dashed") + 
  scale_x_continuous(name = "Time", limits = c(0, 9), breaks = 0:9) + 
  # scale_y_continuous(limits = c(-2.5, 2.5), breaks = seq(-2.5, 2.5, 1)) + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) + 
  scale_color_discrete(name = "ID") 
p1

# Initial realization is completely random 
eta <- rnorm(n, mean = 0, sd = 1.0)
nu0 <- rnorm(n, mean = 0, sd = 0.1)

dfw <- tibble(id = 1:n, 
              eta = eta, 
              y0 = nu0, 
              y1 = eta + rho * y0 + rnorm(n, mean = 0, sd = 0.1),
              y2 = eta + rho * y1 + rnorm(n, mean = 0, sd = 0.1),
              y3 = eta + rho * y2 + rnorm(n, mean = 0, sd = 0.1),
              y4 = eta + rho * y3 + rnorm(n, mean = 0, sd = 0.1),
              y5 = eta + rho * y4 + rnorm(n, mean = 0, sd = 0.1),
              y6 = eta + rho * y5 + rnorm(n, mean = 0, sd = 0.1),
              y7 = eta + rho * y6 + rnorm(n, mean = 0, sd = 0.1),
              y8 = eta + rho * y7 + rnorm(n, mean = 0, sd = 0.1),
              y9 = eta + rho * y8 + rnorm(n, mean = 0, sd = 0.1))

dfl <- tidyr::gather(dfw, key = time, value = y, y0:y9, factor_key = FALSE)
dfl$time <- as.numeric(substr(dfl$time, 2, 2))

p2 <- dfl %>% 
  ggplot(aes(x = time, y = y, color = as.factor(id))) + 
  geom_point() + 
  geom_line() + 
  geom_line(aes(y = (1 - rho)^-1 * eta), linetype = "dashed") + 
  scale_x_continuous(name = "Time", limits = c(0, 9), breaks = 0:9) + 
  # scale_y_continuous(limits = c(-2.5, 2.5), breaks = seq(-2.5, 2.5, 1)) + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) + 
  scale_color_discrete(name = "ID")
p2

# Initial realization is random but correlated with individual effect 
eta <- rnorm(n, mean = 0, sd = 1.0)
nu0 <- 0.5 + eta + rnorm(n, mean = 0, sd = 0.1)

dfw <- tibble(id = 1:n, 
              eta = eta, 
              y0 = nu0, 
              y1 = eta + rho * y0 + rnorm(n, mean = 0, sd = 0.1),
              y2 = eta + rho * y1 + rnorm(n, mean = 0, sd = 0.1),
              y3 = eta + rho * y2 + rnorm(n, mean = 0, sd = 0.1),
              y4 = eta + rho * y3 + rnorm(n, mean = 0, sd = 0.1),
              y5 = eta + rho * y4 + rnorm(n, mean = 0, sd = 0.1),
              y6 = eta + rho * y5 + rnorm(n, mean = 0, sd = 0.1),
              y7 = eta + rho * y6 + rnorm(n, mean = 0, sd = 0.1),
              y8 = eta + rho * y7 + rnorm(n, mean = 0, sd = 0.1),
              y9 = eta + rho * y8 + rnorm(n, mean = 0, sd = 0.1))

dfl <- tidyr::gather(dfw, key = time, value = y, y0:y9, factor_key = FALSE)
dfl$time <- as.numeric(substr(dfl$time, 2, 2))

p3 <- dfl %>% 
  ggplot(aes(x = time, y = y, color = as.factor(id))) + 
  geom_point() + 
  geom_line() + 
  geom_line(aes(y = (1 - rho)^-1 * eta), linetype = "dashed") + 
  scale_x_continuous(name = "Time", limits = c(0, 9), breaks = 0:9) + 
  # scale_y_continuous(limits = c(-2.5, 2.5), breaks = seq(-2.5, 2.5, 1)) + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) + 
  scale_color_discrete(name = "ID")
p3

# Put plots together 
commonp <- ggarrange(p1 + rremove("ylab") + rremove("xlab"), 
                     p2 + rremove("ylab") + rremove("xlab"), 
                     p3 + rremove("ylab") + rremove("xlab"),
                     labels = NULL, 
                     ncol = 3, nrow = 1,
                     common.legend = TRUE, 
                     legend = "right",
                     align = "hv")
annotate_figure(commonp, 
                left = textGrob("", rot = 0),
                bottom = textGrob("Time"))
