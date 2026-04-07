# ------------------------------------------------------------
# derive_mimics_parms.R
# Compute all derived MIMICS parameters
# ------------------------------------------------------------

derive_MIMICS_parms <- function(parms) {
  
  # ---- Protection fractions ----
  fPHYS <- c(
    parms$fPHYS_coeff[1] * exp(parms$fPHYS_exp[1] * parms$fCLAY),
    parms$fPHYS_coeff[2] * exp(parms$fPHYS_exp[2] * parms$fCLAY)
  )
  
  # ---- Desorption ----
  desorb <- parms$desorb_coeff * exp(-parms$desorb_exp * parms$fCLAY)
  
  # ---- Texture modifiers ----
  pSCALAR <- parms$a_texture * exp(-parms$k_texture * sqrt(parms$fCLAY))
  
  MOD2 <- parms$MOD2_base
  MOD2[c(3, 6)] <- c(4 * pSCALAR, 6 * pSCALAR)

  # ---- Return expanded parameter list ----
  modifyList(
    parms,
    list(
      fPHYS = fPHYS,
      desorb = desorb,
      pSCALAR = pSCALAR,
      MOD2 = MOD2
    )
  )
}