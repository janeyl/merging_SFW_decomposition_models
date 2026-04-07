# ------------------------------------------------------------
# Updated MIMICS model ODE
# ------------------------------------------------------------

# Original license for MIMICS code:
# The MIT License (MIT)
# 
# Copyright (c) 2015 will wieder
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#   The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, ITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# Modified by Robert Buchkowski

MIMICS_model <- function(time, state, parms) {
  
  with(as.list(c(state, parms)), {
    # -------------------------------------------------
    # Tree + climate forcing (ONE interface, always)
    # -------------------------------------------------
    tf <- tree_forcing(time)
    
    TotalBiomassTree <- tf["B_tree"]
    litterfall       <- tf["litterfall_f"]
    leaf_mortality   <- tf["leaf_mort"]
    wood_mortality   <- tf["wood_mort"]
    root_mortality   <- tf["root_mort"]
    exudates         <- tf["exudates_f"]
    
    T_t     <- tf["Temp"]    # °C
    theta_t <- tf["theta"]  # m3 m-3
    
    # -------------------------------------------------
    # Optional external biological fluxes (unchanged)
    # -------------------------------------------------
    detritivory_litter  <- if (exists("detritivory_litter")) detritivory_litter else 0
    detritivory_CWD     <- if (exists("detritivory_CWD")) detritivory_CWD else 0
    detritivory_organic <- if (exists("detritivory_organic")) detritivory_organic else 0
    
    faeces  <- if (exists("faeces")) faeces else 0
    carcass <- if (exists("carcass")) carcass else 0
    
    extra_mineral_input_test <-
      if (exists("extra_mineral_input_test")) extra_mineral_input_test else 0
    
    # ------------------------------#
    # MIMICS Model code:
    # ------------------------------#
 
    I <- c(
      litterfall + leaf_mortality + exudates,
      wood_mortality + root_mortality
    )
    
    fMET = I[1]/I[2]
    
    # ---- Temperature‑dependent parameters ----
    Vmax <- exp(T_t * parms$Vslope + parms$Vint) * parms$aV
    Km   <- exp(parms$Kslope * T_t + parms$Kint) * parms$aK
    
    tao <- c(
      parms$tao_lit_coeff[1] * exp(parms$tao_lit_exp[1] * fMET),
      parms$tao_lit_coeff[2] * exp(parms$tao_lit_exp[2] * fMET)
    ) * Tao_MOD1
    
    
    fCHEM <- c(
      parms$fCHEM_coeff[1] * exp(-parms$fCHEM_exp[1] * fMET),
      parms$fCHEM_coeff[2] * exp(-parms$fCHEM_exp[2] * fMET)
    )
    
    fAVAI <- 1 - (fPHYS + fCHEM)
    
    # ---- Final kinetic matrices ----
    VMAX <- Vmax * parms$MOD1
    KM   <- Km / MOD2
    
    # ------------------------------#
    # MIMICS ODEs:
    # ------------------------------#
    
    LITmin = c(NA, NA)
    MICtrn = c(NA, NA, NA)
    SOMmin = c(NA, NA)
    
    #Flows to and from MIC_1
    LITmin[1] = MIC_1 * VMAX[1] * LIT_1 / (KM[1] + LIT_1)   #MIC_1 decomp of MET lit
    LITmin[2] = MIC_1 * VMAX[2] * LIT_2 / (KM[2] + LIT_2)   #MIC_1 decomp of STRUC lit
    MICtrn[1] = MIC_1 * tao[1]  * fPHYS[1]                  #MIC_1 turnover to PHYSICAL SOM 
    MICtrn[2] = MIC_1 * tao[1]  * fCHEM[1]                  #MIC_1 turnover to CHEMICAL SOM  
    MICtrn[3] = MIC_1 * tao[1]  * fAVAI[1]                  #MIC_1 turnover to AVAILABLE SOM  
    SOMmin[1] = MIC_1 * VMAX[3] * SOM_3 / (KM[3] + SOM_3)   #decomp of SOMa by MIC_1
    
    #Flows to and from MIC_2
    LITmin[3] = MIC_2 * VMAX[4] * LIT_1 / (KM[4] + LIT_1)   #decomp of MET litter
    LITmin[4] = MIC_2 * VMAX[5] * LIT_2 / (KM[5] + LIT_2)   #decomp of SRUCTURAL litter
    MICtrn[4] = MIC_2 * tao[2]  * fPHYS[2]                  #MIC_2 turnover to PHYSICAL  SOM 
    MICtrn[5] = MIC_2 * tao[2]  * fCHEM[2]                  #MIC_2 turnover to CHEMICAL  SOM  
    MICtrn[6] = MIC_2 * tao[2]  * fAVAI[2]                  #MIC_2 turnover to AVAILABLE SOM  
    SOMmin[2] = MIC_2 * VMAX[6] * SOM_3 / (KM[6] + SOM_3)   #decomp of SOMa by MIC_2
    
    DEsorb    = SOM_1 * desorb  #* (MIC_1 + MIC_2)		#desorbtion of PHYS to AVAIL (function of fCLAY)
    OXIDAT    = ((MIC_2 * VMAX[5] * SOM_2 / (KO[2]*KM[5] + SOM_2)) +
                   (MIC_1 * VMAX[2] * SOM_2 / (KO[1]*KM[2] + SOM_2)))  #oxidation of C to A
    #can make fluxes from CHEM a function of microbial biomass size?
    
    dLIT_1 = I[1]*(1-FI[1]) - LITmin[1] - LITmin[3]
    dMIC_1 = CUE[1]*(LITmin[1]+ SOMmin[1]) + CUE[2]*(LITmin[2]) - sum(MICtrn[1:3])
    dSOM_1 = I[1]*FI[1] + MICtrn[1] + MICtrn[4]- DEsorb 
    
    dLIT_2 = I[2] * (1-FI[2]) - LITmin[2] - LITmin[4]
    dMIC_2 = CUE[3]*(LITmin[3]+ SOMmin[2]) + CUE[4]*(LITmin[4]) - sum(MICtrn[4:6])  
    dSOM_2 = I[2]*FI[2] + MICtrn[2] + MICtrn[5] - OXIDAT
    
    dSOM_3  = MICtrn[3] + MICtrn[6] + DEsorb + OXIDAT - SOMmin[1] - SOMmin[2]
    
    list(c(dLIT_1, dLIT_2, dMIC_1, dMIC_2, dSOM_1, dSOM_2, dSOM_3))
  })
}
