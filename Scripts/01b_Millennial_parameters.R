# Load in the model parameters:

# Parameter vector
parms <- list(
  force_eqm = 0,
  B0 = 0,
  # Drivers and allocation
  k_monomol = 0.014,      # Saturation control for the monomolecular model, from the reference cited in the model for now.
  TBTmax = 10000,      # Maximum biomass estimated from the portions. Could get higher, but need to explore.
  a_leaf = 0.75*0.1,
  a_wood = 0.75*0.9,
  a_root = 0.25,
  a_root_herb = 0.5,
  # Plant turnover (CBM-inspired)
  k_litterfall_ann = 0.35,   # yr^-1 (foliar litter; HW may approach ~1 yr^-1)
  litter_peak_doy = 288, # Fall
  litter_width_d = 30, # width of peak
  k_mort_leaf  = 0.02/365,   # yr^-1 (non-litterfall leaf mortality)
  k_mort_wood  = 0.025/365,   # yr^-1 (wood mortality to CWD)
  k_mort_root  = 0.08/365,   # yr^-1 (aggregate root mortality)
  k_exudate    = 0.01/365,   # yr^-1 (placeholder; tune to data)
  
  # Fragmentation (physical transfers; donor-controlled)
  k_frag_litter    = 0.8/365,   # Litter -> Organic (placeholder)
  k_frag_CWD    = 0.05/365,      # CWD    -> Organic (placeholder)
  k_frag_organic = 0.08/365,    # Organic-> POM (analogous to AG slow -> BG slow transfer)
  
  # Microbial decomposition in the organic horizon:
  K_ol = 10,
  alpha_ol = 4e10,
  K_ob = 10,
  alpha_ob = 4e10,
  
  k_MICd = 0.036,  # m^2 gC^-1 d^-1; density-dependent mortality
  
  k_l_o = 0.15, # d^-1, higher than in the mineral soil.
  
  # Mortality partitioning fractions
  root_to_organic = 0.7,
  faeces_to_organic= 0.5,
  carcass_to_organic= 0.5,
  herbivory_leaf= 0,
  herbivory_wood= 0,
  herbivory_root= 0,
  detritivory_litter= 0,
  detritivory_CWD = 0,
  detritivory_organic= 0,
  leaf_harvest = 0,
  wood_harvest = 0,
  root_harvest = 0,
  faeces = 0,
  carcass = 0, 
  extra_mineral_input_test = 0,
  
  # Forcings:
  MAT = 6, # From the occupancy model sites
  T_amp = 26, # Set to get summer and winter temperatures correct
  MAtheta = 0.30, # m^3 m^-3; moderately moist soil
  theta_amp = 0.15,
  
  ## --- sorption capacity/soil physical state ---
  depth        = 0.20,    # m; sampling/active layer depth; common 0–30 cm stock. Paper used site-specific depths. [Ref Ambramoff]
  BD           = 1300,    # kg soil m^-3; bulk density from my data and NFI.
  pct_claysilt = 70,      # %; moderate clay+silt content; from my data.
  p_c          = 0.86,    # unitless; coefficient to convert clay+silt to Qmax (Eq. 11). [Ref Ambramoff]
  rho_p        = 3250,    # kg m^-3; particle density; used to estimate porosity from BD; value chosen to get phi_por = 0.6 as in the original model.
  
  ## --- chemistry ---
  pH           = 6,     # typical pH from my data
  
  ## --- moisture scalars ---
  psi_matric   = -15,     # kPa; matric potential ψ. Paper shows λ·φ in Eq. 15 but Table A1 distinguishes total porosity φ and matric potential φ (likely ψ). We use ψ = -15 kPa (field cap-ish). [Ref Ambramoff]
  lambda_mat   = 2.1e-4,  # kPa^-1; dependence of rate on matric potential (λ). Table A1. [Ref Ambramoff]
  k_a_min      = 0.20,    # minimum relative rate in saturated soil (oxygen limitation). Eq. 15; Table A1. [Ref Ambramoff]
  
  ## --- kinetics switches ---
  kinetics     = 1,    # 1 "MM" (default V2), 2 "ECA", or 3 "LIN". See Section 2.3. [Ref Ambramoff]
  
  ## --- Arrhenius constants (temperature sensitivity) ---
  Rgas         = 8.31446, # universal gas constant J K^-1 mol^-1. Table A1. [Ref Ambramoff]
  alpha_pl     = 1.8e12,  # pre-exponential for V_pl (Eq. 3; Table A1). Units consistent with model rate. [Ref Ambramoff]
  Ea_pl        = 63909,   # J mol^-1; activation energy for depolymerization. Table A1. [Ref Ambramoff]
  alpha_lb     = 2.3e12,  # pre-exponential for V_lb (Eq. 14; Table A1). [Ref Ambramoff]
  Ea_lb        = 57865,   # J mol^-1; activation energy for microbial uptake. Table A1. [Ref Ambramoff]
  
  ## --- half-saturation constants ---
  K_pl         = 6443,   # g C m^-2; half-saturation for POM depolymerization (Eq. 2). Table A1 default. [Ref Ambramoff]
  K_lb         = 774.6,     # g C m^-2; half-saturation for microbial uptake (Eq. 13). Table A1 default. [Ref Ambramoff]
  
  ## --- adsorption/desorption ---
  p1           = 0.12,   # coefficient from Mayes et al. (2012) to compute binding affinity term (Eq. 10). [Ref Ambramoff]
  p2           = 0.216,   # second coefficient (Eq. 10). [Ref Ambramoff]
  K_ld         = 1e-3,    # d^-1; we set a small desorption rate coefficient. Table A1 shows mg C L^-1 d^-1,
  # but for a single-layer model we treat it as an effective rate constant.
  # NOTE: Paper lists "F_ld = K_ld * M / Qmax" (Eq. 12), units ambiguous. We adopt
  # d^-1 scaling to keep magnitudes reasonable; adjust as needed with data. [Ref Ambramoff]
  
  ## --- aggregation rates ---
  k_pa         = 0.018,   # d^-1; aggregate formation from POM (Eq. 5). Table A1; Segoli et al., 2013. [Ref Ambramoff]
  k_b          = 0.02,   # d^-1; aggregate breakdown (Eq. 6). Table A1; Segoli et al., 2013. [Ref Ambramoff]
  k_ma         = 0.0048,   # d^-1; aggregate formation from MAOM (Eq. 18). Table A1. [Ref Ambramoff]
  
  ## --- microbial processes ---
  k_bd         = 0.0045,  # m^2 gC^-1 d^-1; density-dependent mortality coefficient (Eq. 16). Table A1; Georgiou et al., 2017. [Ref Ambramoff]
  CUE_ref      = 0.19,    # reference carbon use efficiency (Eqs. 21–22). Table A1. [Ref Ambramoff]
  CUE_T        = 0.012,   # °C^-1; temperature dependence of CUE (linear term). Table A1. [Ref Ambramoff]
  T_ref        = 15,      # °C; reference temperature for CUE. Table A1. [Ref Ambramoff]
  
  ## --- leaching ---
  k_l          = 0.0015,  # d^-1; LMWC leaching loss (Eq. 8). Table A1. [Ref Ambramoff]
  
  ## --- input partitioning ---
  # p_i          = 0.66,    # fraction of C input to POM (Eq. 1). Table A1 (CLM analog). [Ref Ambramoff] # REMOVED BECAUSE I CALCULATE THIS IN THE COUPLED MODEL.
  p_a          = 0.33,    # fraction of aggregate breakdown allocated to POM (Eq. 1; 19). Table A1. [Ref Ambramoff]
  p_b          = 0.50     # necromass partitioning to MAOM (vs L) (Eqs. 7, 19). Table A1. [Ref Ambramoff]
)


parms$phi_por      = 1 - parms$BD/parms$rho_p # total porosity φ (m^3 m^-3); computed from BD & particle density. (Eq. 4) [Ref Ambramoff]
