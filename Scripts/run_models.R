
# Load in library:
library(deSolve)
library(tidyverse)

# ---- Load models ----
source("R/tree_monomolecular.R")
source("R/millennial_model.R")
source("R/millennial_model_detritivory.R")

# ---- Load config utilities ----
source("R/load_config.R")
source("R/make_tree_forcing.R")
source("R/derive_millennial_parms.R")
source("R/init_millennial_state.R")
source("R/plot_ode_out.R")

# ---- Load and merge configs ----
parms  <- load_config("tree_monomolecular")
parms  <- modifyList(parms, yaml::read_yaml("config/millennial.yml"))
parms  <- derive_millennial_parms(parms)

parms$tree_forcing <- make_tree_forcing_equilibrium(parms)

millennial_eqm = rootSolve::stode(
  y     = init_millennial_state(),
  func  = millennial_model,
  parms = parms
)

detritivore_eqm = rootSolve::stode(
  y     = init_millennial_state(T),
  func  = millennial_model_detritivory,
  parms = parms
)

# ---- Build forcing ONCE ----
parms$tree_forcing <- make_tree_forcing(parms)


# ---- Run Millennial model ----
times <- seq(0, 365 * 15, by = 5)

out_millennial <- ode(
  y = millennial_eqm$y,
  times = times, 
  func = millennial_model, 
  parms = parms)


millennial_eqm_det$y

# ---- Run Detritivore model ----
times <- seq(0, 365 * 15, by = 5)
#fist attempt
out_detritivore_eqm <- ode(
  y = detritivore_eqm$y,
  times = times, 
  func = millennial_model_detritivory, 
  parms = parms)

plot_ode_output(out_detritivore_eqm, variable_cols = names(millennial_eqm_det$y))

timeEND <- as.data.frame(out_detritivore_eqm) %>%
  slice_tail(n = 1) %>%
  select(-time) %>%
  unlist()

timeCheck_stable <- tibble(data.frame(out_detritivore_eqm)) %>%
  filter(
    time %in% c(
      max(times),
      max(times)-365,
      max(times)-365*2
    )
  )

#second attempt
times <- seq(0, 365 * 500, by = 5)
out_new <- ode(
  y = timeEND,
  times = times,
  func = millennial_model_detritivory,
  parms = parms
)

plot_ode_output(out_new, variable_cols = names(millennial_eqm_det$y))

timeCheck_stable2 <- tibble(data.frame(out_new)) %>%
  filter(
    time %in% c(
      max(times),
      max(times)-365,
      max(times)-365*2
    )
  )

timeEND2 <- as.data.frame(out_new) %>%
  slice_tail(n = 1) %>%
  select(-time) %>%
  unlist()

#third attemps
#2000 years simulation
times <- seq(0, 365 * 2000, by = 5)
out_new2 <- ode(
  y = timeEND2,
  times = times,
  func = millennial_model_detritivory,
  parms = parms
)

plot_ode_output(out_new2, variable_cols = names(millennial_eqm_det$y))

timeCheck_stable3 <- tibble(data.frame(out_new2)) %>%
  filter(
    time %in% c(
      max(times),
      max(times)-365,
      max(times)-365*2
    )
  )
timeEND3 <- as.data.frame(out_new2) %>%
  slice_tail(n = 1) %>%
  select(-time) %>%
  unlist()

#this is at equilibrium, use as the starting point
#in case it gets lost
timeEND3 <- c(
  Litter      =  94.5075870,
  CWD         = 356.5237003,
  Organic     = 878.2591231,
  DOM         =   4.5505607,
  MIC         =   0.1310657,
  P           = 108.8946685,
  L           =   5.1868545,
  A           = 311.5078745,
  M           = 929.1426153,
  B           =   7.8094434,
  Detritivore =   0.1967198,
  Predator    =   0.1389734
)


# model is stable- do scenarios 

#first make sure that trees and herbs are at equilibrium.
parms2 = parms
parms2$B0 = parms2$TBTmax
parms2$tree_forcing <- make_tree_forcing(parms2)

#reset the forcing function
t(sapply(1:365*3, parms2$tree_forcing)) %>% View()


#create the different scenarios
fullModel <- timeEND3
noPred    <- replace(timeEND3, "Predator", 0) #setting predators = 0 
noPredDet <- replace(timeEND3, c("Predator", "Detritivore"), 0) #setting predators + detritivores = 0 

#simulate the scenarios for just 10 years 
times <- seq(0, 365 * 10, by = 1)
out_new_scenarios <- ode(
  y = timeEND2,
  times = times,
  func = millennial_model_detritivory,
  parms = parms
)

plot_ode_output(out_new2, variable_cols = names(millennial_eqm_det$y))

#but check for 100 years for reviewer 2
times2 <- seq(0, 365 * 100, by = 5)
out_new_scenarios_long <- ode(
  y = timeEND2,
  times = times2,
  func = millennial_model_detritivory,
  parms = parms
)


plot_ode_output(out_new_scenarios_long, variable_cols = names(millennial_eqm_det$y))


# Run the three scenarios - 10 years
out_full <- ode(y = timeEND3,  times = times, func = millennial_model_detritivory, parms = parms2)
out_noPred    <- ode(y = noPred,    times = times, func = millennial_model_detritivory, parms = parms2)
out_noPredDet <- ode(y = noPredDet, times = times, func = millennial_model_detritivory, parms = parms2)


# Run the three scenarios - 100 years

out_full2 <- ode(y = timeEND3,  times = times2, func = millennial_model_detritivory, parms = parms2)
out_noPred2    <- ode(y = noPred,    times = times2, func = millennial_model_detritivory, parms = parms2)
out_noPredDet2 <- ode(y = noPredDet, times = times2, func = millennial_model_detritivory, parms = parms2)


vars <- names(millennial_eqm_det$y)

# Helper to reshape an ode output the same way the function does internally
to_long <- function(out, scenario) {
  df <- as.data.frame(out)
  names(df)[1] <- "time"
  df |>
    select(time, all_of(vars)) |>
    pivot_longer(-time, names_to = "state", values_to = "value") |>
    mutate(scenario = scenario)
}

# Base plot from function (full food web)
p <- plot_ode_output(out_full, variable_cols = vars)

# other two scenarios
p <- p+
  geom_line(data = to_long(out_full,       "Full food web"), aes(time, value, color = scenario)) +
  geom_line(data = to_long(out_noPred,     "No predator"),   aes(time, value, color = scenario)) +
  geom_line(data = to_long(out_noPredDet,  "No pred + det"), aes(time, value, color = scenario)) +
  scale_color_manual(values = c("Full food web" = "black",
                                "No predator"   = "blue",
                                "No pred + det" = "lightgreen")) +
  labs(color = "Scenario")
p
ggsave("outputs/simulations.png",  plot = p,  width = 9.5, height = 6, dpi = 300)

#plot the 100 year simulation scenarios, but 
# Keep only one day per year (default: day 1)
keep_yearly <- function(out, day_of_year = 1) {
  out[out[, "time"] %% 365 == day_of_year, ]
}

p2 <- plot_ode_output(keep_yearly(out_full2), variable_cols = vars)

p2 <- p2 +
  geom_line(data = to_long(keep_yearly(out_full2),      "Full food web"), aes(time, value, color = scenario)) +
  geom_line(data = to_long(keep_yearly(out_noPred2),    "No predator"),   aes(time, value, color = scenario)) +
  geom_line(data = to_long(keep_yearly(out_noPredDet2), "No pred + det"), aes(time, value, color = scenario)) +
  scale_color_manual(values = c("Full food web" = "black",
                                "No predator"   = "blue",
                                "No pred + det" = "lightgreen")) +
  labs(color = "Scenario")
p2

#ggsave("outputs/simulations_100yrs.png",  plot = p2,  width = 9.5, height = 6, dpi = 300)


# Combine the three scenario outputs into one long data frame
to_long <- function(out, scenario) {
  as.data.frame(out) |>
    pivot_longer(-time, names_to = "Pool", values_to = "Value") |>
    mutate(Scenario = scenario)
}

scenario_long <- bind_rows(
  to_long(out_full,      "Full model"),
  to_long(out_noPred,    "No predators"),
  to_long(out_noPredDet, "No predators/detritivores")
)

scenario_colors <- c("No predators"              = "blue",
                     "No predators/detritivores" = "lightgreen")

# ---- Fig 1: endpoint (max time) difference from Full model ----
full_endpoint <- scenario_long %>%
  filter(Scenario == "Full model", time == max(time)) %>%
  select(Pool, full_val = Value)

diff_end <- scenario_long %>%
  filter(Scenario != "Full model", time == max(time)) %>%
  left_join(full_endpoint, by = "Pool") %>%
  mutate(Diff = Value - full_val)

fig_end <- ggplot(diff_end, aes(x = Scenario, y = Diff, fill = Scenario)) +
  geom_col(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  facet_wrap(~ Pool, scales = "free_y") +
  scale_fill_manual(values = scenario_colors) +
  labs(y = expression(Delta ~ "endpoint pool size"),
       x = NULL, fill = NULL) +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

fig_end
ggsave("outputs/fig_end.png",  plot = fig_end,  width = 8, height = 6, dpi = 300)

# ---- Fig 2: whole-trajectory mean difference from Full model ----
full_mean <- scenario_long %>%
  filter(Scenario == "Full model") %>%
  group_by(Pool) %>%
  summarise(full_mean = mean(Value), .groups = "drop")

diff_mean <- scenario_long %>%
  filter(Scenario != "Full model") %>%
  group_by(Scenario, Pool) %>%
  summarise(mean_value = mean(Value), .groups = "drop") %>%
  left_join(full_mean, by = "Pool") %>%
  mutate(Diff = mean_value - full_mean)

fig_mean <- ggplot(diff_mean, aes(x = Scenario, y = Diff, fill = Scenario)) +
  geom_col(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  facet_wrap(~ Pool, scales = "free_y") +
  scale_fill_manual(values = scenario_colors) +
  labs(y = expression(Delta ~ "mean pool size"),
       x = NULL, fill = NULL) +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

fig_mean
ggsave("outputs/fig_mean.png", plot = fig_mean, width = 8, height = 6, dpi = 300)


#effect sizes

# Per-scenario, per-pool summary (endpoint + time-mean)
summary_df <- scenario_long %>%
  group_by(Scenario, Pool) %>%
  summarise(
    endpoint  = Value[time == max(time)],
    time_mean = mean(Value),
    .groups = "drop"
  )

# Baseline values (Full model)
baseline <- summary_df %>%
  filter(Scenario == "Full model") %>%
  select(Pool, baseline_end = endpoint, baseline_mean = time_mean)

# Effect sizes for the removal scenarios
effect_sizes <- summary_df %>%
  filter(Scenario != "Full model") %>%
  left_join(baseline, by = "Pool") %>%
  mutate(
    # Endpoint effects
    delta_end       = endpoint - baseline_end,
    pct_end         = 100 * delta_end / baseline_end,
    LRR_end         = log(endpoint / baseline_end),
    # Time-averaged effects
    delta_mean      = time_mean - baseline_mean,
    pct_mean        = 100 * delta_mean / baseline_mean,
    LRR_mean        = log(time_mean / baseline_mean)
  ) %>%
  select(Scenario, Pool,
         baseline_end, endpoint, delta_end, pct_end, LRR_end,
         baseline_mean, time_mean, delta_mean, pct_mean, LRR_mean)

effect_sizes
write.csv(effect_sizes, "outputs/effect_sizes.csv", row.names = FALSE)

effect_sizes <- effect_sizes %>%
  filter(!(Scenario == "No predators"              & Pool == "Predator"),
         !(Scenario == "No predators/detritivores" & Pool %in% c("Predator", "Detritivore")))

library(scales)

percent_change <- ggplot(effect_sizes, aes(x = Scenario, y = Pool, fill = pct_mean)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%+.1f%%", pct_mean)), size = 3) +
  scale_fill_gradient2(
    low = "#2C7BB6", mid = "white", high = "#D7191C",
    midpoint = 0,
    limits = c(-50, 50),     # cap the color scale here
    oob = squish,            # anything beyond limits gets the extreme color
    name = "% change\n(time-mean)"
  ) +
  labs(x = NULL, y = NULL,
       title = "Effect of soil animals on soil C pools") +
  theme_minimal() +
  theme(panel.grid = element_blank())

ggsave("outputs/percent_change.png", plot = percent_change, width = 4.5, height = 7, dpi = 300)

#Adding in standard deviations model scenarios----

#----Run Up----


# ---- user-set options --------------------------------------------------
NTOT         <- 2                 # number of parameter draws (test value)
yts          <- 5                 # years per integration chunk
tmax2        <- 100               # scenario simulation length (days)
sdlog_val    <- 0.1536            # lognormal SD on animal parameters
rel_tol      <- 1e-3              # relative convergence tolerance (0.1%/yr)
max_attempts <- 50                # safety cap (up to ~250 yrs per run)

# Animal parameters to perturb
animal_params <- c("d_detritivores")
# When ready, expand to: c("d_detritivores", "a_detritivores", "p_detritivores",
#                          "c_detritivores", "d_predator", "a_predator",
#                          "p_predator",     "c_predator")

eq_times <- seq(1, 365 * yts, by = 73)

# ---- singlerun ---------------------------------------------------------
singlerun <- function(idx) {
  cat(sprintf("--- idx %d ---\n", idx))
  
  ID  <- round(runif(1) * 1e8, 0)
  ID2 <- round(runif(1) * 1e8, 0)
  
  paramscur <- parms2
  
  # Draw parameters ONCE per run (fixed across the equilibration chunks)
  for (nm in animal_params) {
    paramscur[[nm]] <- rlnorm(1,
                              meanlog = log(parms2[[nm]]),
                              sdlog   = sdlog_val)
  }
  
  # --- iteratively extend integration until equilibrium ---
  ystable <- timeEND3
  ERRRRR  <- TRUE
  attempt <- 0
  
  while (ERRRRR && attempt < max_attempts) {
    attempt <- attempt + 1
    
    stablerun <- ode(y = ystable, times = eq_times,
                     func = millennial_model_detritivory, parms = paramscur)
    
    if (nrow(stablerun) == length(eq_times)) {
      n        <- nrow(stablerun)
      ystable  <- stablerun[n, -1]   # continue from last state (same params)
      
      last_t   <- stablerun[n, "time"]
      prev_idx <- which.min(abs(stablerun[, "time"] - (last_t - 365)))
      
      if (last_t - stablerun[prev_idx, "time"] >= 300) {
        rel_delta <- max(abs(stablerun[n, -1] - stablerun[prev_idx, -1]) /
                           pmax(abs(stablerun[n, -1]), 1e-6))
        ERRRRR <- rel_delta > rel_tol
        cat(sprintf("  attempt %d  rel_delta = %.4g  %s\n",
                    attempt, rel_delta,
                    if (ERRRRR) "(still drifting, extending)" else "(equilibrated)"))
      }
    } else {
      cat(sprintf("  attempt %d  ode() returned wrong number of rows\n", attempt))
    }
  }
  
  if (ERRRRR) {
    warning(sprintf("idx %d never equilibrated after %d attempts",
                    idx, max_attempts))
    return(NULL)
  }
  ystable2 <- ystable
  
  # --- three scenarios (initial conditions) ---
  y_full      <- ystable
  y_noPred    <- ystable; y_noPred["Predator"] <- 0
  y_noPredDet <- ystable; y_noPredDet[c("Predator", "Detritivore")] <- 0
  
  out_full      <- ode(y = y_full,      times = 1:tmax2,
                       func = millennial_model_detritivory, parms = paramscur)
  out_noPred    <- ode(y = y_noPred,    times = 1:tmax2,
                       func = millennial_model_detritivory, parms = paramscur)
  out_noPredDet <- ode(y = y_noPredDet, times = 1:tmax2,
                       func = millennial_model_detritivory, parms = paramscur)
  
  # --- stack and annotate ---
  out  <- rbind(out_full, out_noPred, out_noPredDet)
  out1 <- as.data.frame(out)
  out1$Scenario <- rep(c("Full", "NoPred", "NoPredDet"),
                       times = c(nrow(out_full),
                                 nrow(out_noPred),
                                 nrow(out_noPredDet)))
  out1$Run  <- ID
  out1$Run2 <- ID2
  attr(out1, "PARS")    <- unlist(paramscur[animal_params])
  attr(out1, "YSTABLE") <- ystable2
  
  return(out1)
}

# ---- run + save --------------------------------------------------------
repseq <- seq_len(NTOT)
out1   <- lapply(repseq, singlerun)

# Drop any NULLs from runs that failed to equilibrate
out1   <- out1[!vapply(out1, is.null, logical(1))]
outf   <- do.call("rbind", out1)

fname  <- paste0("outputs/modelcluster_",
                 round(runif(1), 7) * 1e7,
                 round(runif(1), 7) * 1e7,
                 ".csv")
write.csv(outf, fname, row.names = FALSE)

#----END----