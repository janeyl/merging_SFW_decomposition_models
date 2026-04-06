derive_millennial_parms <- function(parms) {
  
  parms$phi_por <- 1 - parms$BD / parms$rho_p
  
  return(parms)
}