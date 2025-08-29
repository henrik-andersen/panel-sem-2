
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


# Time-invariant DGP ------------------------------------------------------

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


# Models: time invariant ---

mri <- "
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

mdp <- "
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

# Data: OBS-TI, Models: TI
mdp_obs_ti_ti.fit <- sem(mdp, data = df_obs_ti, estimator = "ML")
mri_obs_ti_ti.fit <- sem(mri, data = df_obs_ti, estimator = "ML")

summary(mdp_obs_ti_ti.fit)
summary(mri_obs_ti_ti.fit)

# Data: RES-TI, Models: TI
mdp_res_ti_ti.fit <- sem(mdp, data = df_res_ti, estimator = "ML")
mri_res_ti_ti.fit <- sem(mri, data = df_res_ti, estimator = "ML")

summary(mdp_res_ti_ti.fit)
summary(mri_res_ti_ti.fit)


# Models: partially time invariant ---

mri <- "
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

mdp <- "
  # Individual effects
  ux =~ 1*x2 + 1*x3 + 1*x4 + 1*x5
  uy =~ 1*y2 + 1*y3 + 1*y4 + 1*y5
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

# Data: OBS-TI, Models: PTV
mdp_obs_ti_ptv.fit <- sem(mdp, data = df_obs_ti, estimator = "ML")
mri_obs_ti_ptv.fit <- sem(mri, data = df_obs_ti, estimator = "ML")

summary(mdp_obs_ti_ptv.fit)
summary(mri_obs_ti_ptv.fit)

# Data: RES-TI, Models: PTI
mdp_res_ti_ptv.fit <- sem(mdp, data = df_res_ti, estimator = "ML")
mri_res_ti_ptv.fit <- sem(mri, data = df_res_ti, estimator = "ML")

summary(mdp_res_ti_ptv.fit)
summary(mri_res_ti_ptv.fit)


# Models: Fully time invariant ---

mri <- "
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

mdp <- "
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

# Data: OBS-TI, Models: FTV
mdp_obs_ti_ftv.fit <- sem(mdp, data = df_obs_ti, estimator = "ML")
mri_obs_ti_ftv.fit <- sem(mri, data = df_obs_ti, estimator = "ML")

summary(mdp_obs_ti_ftv.fit); mean_summary_func(mdp_obs_ti_ftv.fit)
summary(mri_obs_ti_ftv.fit); mean_summary_func(mri_obs_ti_ftv.fit)


# Data: RES-TI, Models: FTV
mdp_res_ti_ftv.fit <- sem(mdp, data = df_res_ti, estimator = "ML")
mri_res_ti_ftv.fit <- sem(mri, data = df_res_ti, estimator = "ML")

summary(mdp_res_ti_ftv.fit); mean_summary_func(mdp_res_ti_ftv.fit)
summary(mri_res_ti_ftv.fit); mean_summary_func(mri_res_ti_ftv.fit)



# Partially time invariant DGP --------------------------------------------

sim_func_ptv <- function(n, 
                        sd_ux, sd_uy, 
                        sd_ex, sd_ey,
                        phi, gam, 
                        bet, rho,
                        obs = TRUE) {
  
  ux <- rnorm(n, 0, sd_ux)
  uy <- rnorm(n, 0, sd_uy)
  
  phi <- rnorm(n = 6, mean = phi, sd = 0.01)
  gam <- rnorm(n = 6, mean = gam, sd = 0.01)
  bet <- rnorm(n = 6, mean = bet, sd = 0.01)
  rho <- rnorm(n = 6, mean = rho, sd = 0.01)
  
  if(obs == TRUE) {
    
    df <- tibble(
      # x0 = 1 * ux + rnorm(n, 0, sd_ex),
      # y0 = 1 * uy + rnorm(n, 0, sd_ey),
      x0 = rnorm(n, 0, sd_ex),
      y0 = rnorm(n, 0, sd_ey),
      x1 = phi[2] * x0 + gam[2] * y0 + 1 * ux + rnorm(n, 0, sd_ex),
      y1 = bet[2] * x0 + rho[2] * y0 + 1 * uy + rnorm(n, 0, sd_ey),
      x2 = phi[3] * x1 + gam[3] * y1 + 1 * ux + rnorm(n, 0, sd_ex),
      y2 = bet[3] * x1 + rho[3] * y1 + 1 * uy + rnorm(n, 0, sd_ey),
      x3 = phi[4] * x2 + gam[4] * y2 + 1 * ux + rnorm(n, 0, sd_ex),
      y3 = bet[4] * x2 + rho[4] * y2 + 1 * uy + rnorm(n, 0, sd_ey),
      x4 = phi[5] * x3 + gam[5] * y3 + 1 * ux + rnorm(n, 0, sd_ex),
      y4 = bet[5] * x3 + rho[5] * y3 + 1 * uy + rnorm(n, 0, sd_ey),
      x5 = phi[6] * x4 + gam[6] * y4 + 1 * ux + rnorm(n, 0, sd_ex),
      y5 = bet[6] * x4 + rho[6] * y4 + 1 * uy + rnorm(n, 0, sd_ey)
    )
    
  } else {
    
    df <- tibble(
      # x0 = 1 * ux + rnorm(n, 0, sd_ex),
      # y0 = 1 * uy + rnorm(n, 0, sd_ey),
      x0 = rnorm(n, 0, sd_ex),
      y0 = rnorm(n, 0, sd_ey),
      x1 = phi[2] * x0 + gam[2] * y0 + (1 - phi[2]) * ux - gam[2] * uy + rnorm(n, 0, sd_ex),
      y1 = bet[2] * x0 + rho[2] * y0 - bet[2] * ux + (1 - rho[2]) * uy + rnorm(n, 0, sd_ey), 
      x2 = phi[3] * x1 + gam[3] * y1 + (1 - phi[3]) * ux - gam[3] * uy + rnorm(n, 0, sd_ex),
      y2 = bet[3] * x1 + rho[3] * y1 - bet[3] * ux + (1 - rho[3]) * uy + rnorm(n, 0, sd_ey), 
      x3 = phi[4] * x2 + gam[4] * y2 + (1 - phi[4]) * ux - gam[4] * uy + rnorm(n, 0, sd_ex),
      y3 = bet[4] * x2 + rho[4] * y2 - bet[4] * ux + (1 - rho[4]) * uy + rnorm(n, 0, sd_ey), 
      x4 = phi[5] * x3 + gam[5] * y3 + (1 - phi[5]) * ux - gam[5] * uy + rnorm(n, 0, sd_ex),
      y4 = bet[5] * x3 + rho[5] * y3 - bet[5] * ux + (1 - rho[5]) * uy + rnorm(n, 0, sd_ey), 
      x5 = phi[6] * x4 + gam[6] * y4 + (1 - phi[6]) * ux - gam[6] * uy + rnorm(n, 0, sd_ex),
      y5 = bet[6] * x4 + rho[6] * y4 - bet[6] * ux + (1 - rho[6]) * uy + rnorm(n, 0, sd_ey)
    )
    
  }
  
  return(df)
  
}

df_obs_ptv <- sim_func_ptv(n = 1000000, 
                           sd_ux = 1, sd_uy = 1, 
                           sd_ex = 1, sd_ey = 1, 
                           phi = 0.2, gam = 0.5,
                           bet = 0.5, rho = 0.2, 
                           obs = TRUE)
head(df_obs_ti)

df_res_ptv <- sim_func_ptv(n = 1000000, 
                           sd_ux = 1, sd_uy = 1, 
                           sd_ex = 1, sd_ey = 1, 
                           phi = 0.2, gam = 0.5,
                           bet = 0.5, rho = 0.2, 
                           obs = FALSE)
head(df_res_ti)


# Models: time invariant ---

mri <- "
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

mdp <- "
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

# Data: OBS-PTV, Models: TI
mdp_obs_ptv_ti.fit <- sem(mdp, data = df_obs_ptv, estimator = "ML")
mri_obs_ptv_ti.fit <- sem(mri, data = df_obs_ptv, estimator = "ML")

summary(mdp_obs_ptv_ti.fit)
summary(mri_obs_ptv_ti.fit)

# Data: RES-PTV, Models: TI
mdp_res_ptv_ti.fit <- sem(mdp, data = df_res_ptv, estimator = "ML")
mri_res_ptv_ti.fit <- sem(mri, data = df_res_ptv, estimator = "ML")

summary(mdp_res_ptv_ti.fit)
summary(mri_res_ptv_ti.fit)


# Models: partially time invariant ---

mri <- "
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

mdp <- "
  # Individual effects
  ux =~ 1*x2 + 1*x3 + 1*x4 + 1*x5
  uy =~ 1*y2 + 1*y3 + 1*y4 + 1*y5
  # Regressions (time invariant)
  x2 ~ phi2*x1 + gam2*y1
  y2 ~ bet2*x1 + rho2*y1
  x3 ~ phi3*x2 + gam3*y2
  y3 ~ bet3*x2 + rho3*y2
  x4 ~ phi4*x3 + gam4*y3
  y4 ~ bet4*x3 + rho4*y3
  x5 ~ phi4*x4 + gam5*y4
  y5 ~ bet4*x4 + rho5*y4
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

# Data: OBS-PTV, Models: PTV
mdp_obs_ptv_ptv.fit <- sem(mdp, data = df_obs_ptv, estimator = "ML")
mri_obs_ptv_ptv.fit <- sem(mri, data = df_obs_ptv, estimator = "ML")

summary(mdp_obs_ptv_ptv.fit)
summary(mri_obs_ptv_ptv.fit)

# Data: RES-PTV, Models: PTV
mdp_res_ptv_ptv.fit <- sem(mdp, data = df_res_ptv, estimator = "ML")
mri_res_ptv_ptv.fit <- sem(mri, data = df_res_ptv, estimator = "ML")

summary(mdp_res_ptv_ptv.fit)
summary(mri_res_ptv_ptv.fit)


# Models: Fully time invariant ---

mri <- "
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

mdp <- "
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
  x5 ~ phi4*x4 + gam5*y4
  y5 ~ bet4*x4 + rho5*y4
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

# Data: OBS-PTV, Models: FTV
mdp_obs_ptv_ftv.fit <- sem(mdp, data = df_obs_ti, estimator = "ML")
mri_obs_ptv_ftv.fit <- sem(mri, data = df_obs_ti, estimator = "ML")

summary(mdp_obs_ptv_ftv.fit)
summary(mri_obs_ptv_ftv.fit)

# Data: RES-PTV, Models: FTV
mdp_res_ptv_ftv.fit <- sem(mdp, data = df_res_ti, estimator = "ML")
mri_res_ptv_ftv.fit <- sem(mri, data = df_res_ti, estimator = "ML")

summary(mdp_res_ptv_ftv.fit)
summary(mri_res_ptv_ftv.fit)




# Fully time invariant DGP ------------------------------------------------

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

df_obs_ftv <- sim_func_ftv(n = 1000000, 
                           sd_ux = 1, sd_uy = 1, 
                           sd_ex = 1, sd_ey = 1, 
                           phi = 0.2, gam = 0.5,
                           bet = 0.5, rho = 0.2, 
                           obs = TRUE)
head(df_obs_ftv)

df_obs_ftv %>% 
  summarise(across(starts_with("x"), mean))
df_obs_ftv %>% 
  summarise(across(starts_with("y"), mean))

df_res_ftv <- sim_func_ftv(n = 1000000, 
                           sd_ux = 1, sd_uy = 1, 
                           sd_ex = 1, sd_ey = 1, 
                           phi = 0.2, gam = 0.5,
                           bet = 0.5, rho = 0.2, 
                           obs = FALSE)
head(df_res_ftv)

df_res_ftv %>% 
  summarise(across(starts_with("x"), mean))
df_res_ftv %>% 
  summarise(across(starts_with("y"), mean))

# Models: time invariant ---

mri <- "
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

mdp <- "
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

# Data: OBS-FTV, Models: TI
mdp_obs_ftv_ti.fit <- sem(mdp, data = df_obs_ftv, estimator = "ML")
mri_obs_ftv_ti.fit <- sem(mri, data = df_obs_ftv, estimator = "ML")

summary(mdp_obs_ftv_ti.fit); mean_summary_func(mdp_obs_ftv_ti.fit)
summary(mri_obs_ftv_ti.fit); mean_summary_func(mri_obs_ftv_ti.fit)

# Data: RES-FTV, Models: TI
mdp_res_ftv_ti.fit <- sem(mdp, data = df_res_ftv, estimator = "ML")
mri_res_ftv_ti.fit <- sem(mri, data = df_res_ftv, estimator = "ML")

summary(mdp_res_ftv_ti.fit); mean_summary_func(mdp_res_ftv_ti.fit) 
summary(mri_res_ftv_ti.fit); mean_summary_func(mri_res_ftv_ti.fit) 


# Models: partially time invariant ---

mri <- "
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

mdp <- "
  # Individual effects
  ux =~ 1*x2 + 1*x3 + 1*x4 + 1*x5
  uy =~ 1*y2 + 1*y3 + 1*y4 + 1*y5
  # Regressions (time invariant)
  x2 ~ phi2*x1 + gam2*y1
  y2 ~ bet2*x1 + rho2*y1
  x3 ~ phi3*x2 + gam3*y2
  y3 ~ bet3*x2 + rho3*y2
  x4 ~ phi4*x3 + gam4*y3
  y4 ~ bet4*x3 + rho4*y3
  x5 ~ phi4*x4 + gam5*y4
  y5 ~ bet4*x4 + rho5*y4
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

# Data: OBS-FTV, Models: PTV
mdp_obs_ftv_ptv.fit <- sem(mdp, data = df_obs_ftv, estimator = "ML")
mri_obs_ftv_ptv.fit <- sem(mri, data = df_obs_ftv, estimator = "ML")

summary(mdp_obs_ftv_ptv.fit)
summary(mri_obs_ftv_ptv.fit)

# Data: RES-FTV, Models: PTV
mdp_res_ftv_ptv.fit <- sem(mdp, data = df_res_ftv, estimator = "ML")
mri_res_ftv_ptv.fit <- sem(mri, data = df_res_ftv, estimator = "ML")

summary(mdp_res_ftv_ptv.fit)
summary(mri_res_ftv_ptv.fit)


# Models: Fully time invariant ---

mri <- "
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

mdp <- "
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

# Data: OBS-FTV, Models: FTV
mdp_obs_ftv_ftv.fit <- sem(mdp, data = df_obs_ftv, estimator = "ML")
mri_obs_ftv_ftv.fit <- sem(mri, data = df_obs_ftv, estimator = "ML")

summary(mdp_obs_ftv_ftv.fit); mean_summary_func(mdp_obs_ftv_ftv.fit)
summary(mri_obs_ftv_ftv.fit); mean_summary_func(mri_obs_ftv_ftv.fit)

# Data: RES-FTV, Models: FTV
mdp_res_ftv_ftv.fit <- sem(mdp, data = df_res_ftv, estimator = "ML")
mri_res_ftv_ftv.fit <- sem(mri, data = df_res_ftv, estimator = "ML")

summary(mdp_res_ftv_ftv.fit); mean_summary_func(mdp_res_ftv_ftv.fit)
summary(mri_res_ftv_ftv.fit); mean_summary_func(mri_res_ftv_ftv.fit)

anova(mdp_obs_ftv_ftv.fit, mri_obs_ftv_ftv.fit)
anova(mdp_res_ftv_ftv.fit, mri_res_ftv_ftv.fit)

mri_obs <- "
  # Individual effects
  ux =~ lx1*x1 + lx2*x2 + lx3*x3 + lx4*x4 + lx5*x5
  uy =~ ly1*y1 + ly2*y2 + ly3*y3 + ly4*y4 + ly5*y5
  # Regressions (time invariant)
  x2 ~ phi2*x1 + gam2*y1 + cxux2*ux + cxuy2*uy
  y2 ~ bet2*x1 + rho2*y1 + cyuy2*uy + cyux2*ux
  x3 ~ phi3*x2 + gam3*y2 + cxux3*ux + cxuy3*uy
  y3 ~ bet3*x2 + rho3*y2 + cyuy3*uy + cyux3*ux
  x4 ~ phi4*x3 + gam4*y3 + cxux4*ux + cxuy4*uy
  y4 ~ bet4*x3 + rho4*y3 + cyuy4*uy + cyux4*ux
  x5 ~ phi5*x4 + gam5*y4 + cxux5*ux + cxuy5*uy
  y5 ~ bet5*x4 + rho5*y4 + cyuy5*uy + cyux5*ux
  # Correlations   
  ux ~~ uy + 0*x1 + 0*y1
  uy ~~ 0*x1 + 0*y1
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
  # Need to constrain to zero (bug?)
  # x5w ~~ 0*y5w
  # Constraints
  cxux2 == lx2 - phi2 * lx1
  cxux3 == lx3 - phi3 * lx2
  cxux4 == lx4 - phi4 * lx3
  cxux5 == lx5 - phi5 * lx4
  cxuy2 == -gam2 * ly1
  cxuy3 == -gam3 * ly2
  cxuy4 == -gam4 * ly3
  cxuy5 == -gam5 * ly4
  cyuy2 == lx2 - phi2 * lx1
  cyuy3 == lx3 - phi3 * lx2
  cyuy4 == lx4 - phi4 * lx3
  cyuy5 == lx5 - phi5 * lx4
  cyux2 == -bet2 * lx1
  cyux3 == -bet3 * lx2
  cyux4 == -bet4 * lx3
  cyux5 == -bet5 * lx4
"
mri_obs_obs_ftv_ftv.fit <- sem(mri_obs, data = df_obs_ftv, estimator = "ML")
summary(mri_obs_obs_ftv_ftv.fit)

mri_obs_res_ftv_ftv.fit <- sem(mri_obs, data = df_res_ftv, estimator = "ML")


mco <- "
  # Individual effects
  ux =~ lx2*x2 + lx3*x3 + lx4*x4 + lx5*x5
  uy =~ ly2*y2 + ly3*y3 + ly4*y4 + ly5*y5
  # Regressions (time invariant)
  x2 ~ phi2*x1 + gam2*y1
  y2 ~ bet2*x1 + rho2*y1
  x3 ~ phi3*x2 + gam3*y2
  y3 ~ bet3*x2 + rho3*y2
  x4 ~ phi4*x3 + gam4*y3
  y4 ~ bet4*x3 + rho4*y3
  x5 ~ phi4*x4 + gam5*y4
  y5 ~ bet4*x4 + rho5*y4
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


