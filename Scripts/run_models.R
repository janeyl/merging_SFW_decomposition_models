
# Load in library:
library(deSolve)
library(tidyverse)

# ---- Load models ----
source("R/tree_monomolecular.R")
source("R/millennial_model.R")

# ---- Load config utilities ----
source("R/load_config.R")
source("R/make_tree_forcing.R")
source("R/derive_millennial_parms.R")
source("R/init_millennial_state.R")

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

# ---- Build forcing ONCE ----
parms$tree_forcing <- make_tree_forcing(parms)


# ---- Run model ----
times <- seq(0, 365 * 10, by = 1)

out_millennial <- ode(
  y = millennial_eqm$y,
  times = times, 
  func = millennial_model, 
  parms = parms)



# ---- Load code ----
source("R/tree_monomolecular.R")
source("R/make_tree_forcing.R")
source("R/century_model.R")
source("R/load_config.R")
source("R/init_century_state.R")

# ---- Load parameters ----
parms <- load_config("tree_monomolecular")
parms <- modifyList(parms, yaml::read_yaml("config/century.yml"))

#---- Calculate equilibrium ----
parms$tree_forcing <- make_tree_forcing_equilibrium(parms)

century_eqm = rootSolve::stode(
  y     = init_century_state(),
  func  = century_model,
  parms = parms
)


# ---- Build forcing once ----
parms$tree_forcing <- make_tree_forcing(parms)

# ---- Run ----
out_century <- ode(
  y     = century_eqm$y,
  times = times,
  func  = century_model,
  parms = parms
)


# ---- Load code ----
source("R/tree_monomolecular.R")
source("R/make_tree_forcing.R")
source("R/MIMICS_model.R")
source("R/derive_MIMICS_parms.R")
source("R/load_config.R")
source("R/init_MIMICS_state.R")

parms <- load_config("tree_monomolecular")
parms <- modifyList(parms, yaml::read_yaml("config/MIMICS.yml"))
parms <- derive_MIMICS_parms(parms)

#---- Calculate equilibrium ----
parms$tree_forcing <- make_tree_forcing_equilibrium(parms)

MIMICS_eqm = rootSolve::stode(
  y     = init_MIMICS_state(),
  func  = MIMICS_model,
  parms = parms
)


# ---- Build forcing once ----
parms$tree_forcing <- make_tree_forcing(parms)

# ---- Run ----
out_MIMICS <- ode(
  y     = MIMICS_eqm$y,
  times = times,
  func  = MIMICS_model,
  parms = parms
)


tibble(data.frame(out_millennial))

tibble(data.frame(out_century))
