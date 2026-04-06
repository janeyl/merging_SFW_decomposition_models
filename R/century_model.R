
# ------------------------------------------------------------
# Updated Century
# ------------------------------------------------------------
#Author: Robert Buchkowski after Rose Abramoff, after Parton et al. (1987)
#Date: Jan 25, 2026; Sep 11, 2021

#This function contains the system of equations for the Century model, which was developed in Parton et al. (1987). 
#The equation numbers correspond to those in Abramoff et al. (2021) Appendix B, with parameters defined in Table A2.

century_model <- function(time, state, parms){
  
  with(c(state, parms), {
    # ----------------------------
    # ---- Get vegetation + climate forcing ----
    # ----------------------------
    forcing <- tree_forcing(time)
    
    # Extract required drivers
    TotalBiomassTree <- forcing["B_tree"]
    litterfall       <- forcing["litterfall_f"]
    leaf_mortality   <- forcing["leaf_mort"]
    wood_mortality   <- forcing["wood_mort"]
    root_mortality   <- forcing["root_mort"]
    exudates         <- forcing["exudates_f"]
    # ----------------------------
    # Forcings at current time (Fi, T, theta)
    # ----------------------------
    T_t              <- forcing["Temp"] # °C
    theta_t          <- forcing["theta"] # m^3 m^-3
    
    # ----------------------------------#
    # Abiotic scalars for decomposition:
    # ----------------------------------#
    
    #Equation B1
    t_scalar <- (t2 + (t3 / pi) * atan(pi * t4 * (T_t - t1))) /
      (t2 + (t3 / pi) * atan(pi * t4 *(30.0 - t1)))
    
    #Equation B2
    w_scalar <- 1.0 / (1.0 + w1 * exp(-w2 * theta_t/field_cap))
    
    #Equation B3
    f_TEX = c1 - c2*pct_claysilt*0.01
    
    #Equation B4
    f_StrLitter = StrLitter * k_strlitter * t_scalar * w_scalar * exp(-3*LigFrac)
    
    #Equation B5
    f_MetLitter = MetLitter * k_metlitter * t_scalar * w_scalar  
    
    #Equation B6
    f_ACTIVE <- ACTIVE * k_active * t_scalar * w_scalar * f_TEX
    
    #Equation B7 
    f_SLOW <- SLOW * k_slow * t_scalar * w_scalar
    
    #Equation B8
    f_PASSIVE <- PASSIVE * k_passive * t_scalar * w_scalar
    
    #Equation B9
    dStrLitter = wood_mortality + root_to_organic*root_mortality - f_StrLitter
    
    #Equation B10
    dMetLitter = litterfall + leaf_mortality + exudates - f_MetLitter
    
    #Equation B11
    dACTIVE <- (1 - root_to_organic)*root_mortality + (1-LigFrac) * strlitter_to_active * f_StrLitter + metlitter_to_active * f_MetLitter  + f_SLOW * slow_to_active + f_PASSIVE * passive_to_active - f_ACTIVE
    
    #Equation B12
    dSLOW <-  LigFrac * strlitter_to_slow * f_StrLitter + f_ACTIVE * (1-f_TEX-active_to_passive) - f_SLOW
    
    #Equation B13
    dPASSIVE <- f_ACTIVE * active_to_passive + f_SLOW * slow_to_passive - f_PASSIVE

    # ---------------------------
    # Return list for deSolve
    # ---------------------------
    list(
      c(dStrLitter, dMetLitter, dACTIVE, dSLOW, dPASSIVE)
    )
  })
}


