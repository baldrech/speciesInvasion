# mizerEncounter with temperature dependence on search volume


mizerEncounterTemp <- function(params, n, n_pp, n_other, t, ...) {
  
  # idx_sp are the index values of params@w_full such that
  # params@w_full[idx_sp] = params@w
  idx_sp <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
  
  # temperature component
  scalar <- tempFun(params@w, temperature = params@other_params$other$temperature[t+1], # t starts at 0
                    t_d = params@species_params$t_d,
                    t_ref = params@species_params$t_ref,
                    Ea = params@species_params$Ea,
                    c_a = params@species_params$c_a,
                    Ed = params@species_params$Ed,
                    c_d = params@species_params$c_d)
  
  scalar <- t(drop(scalar))
  
  # If the the user has set a custom pred_kernel we can not use fft.
  # In this case we use the code from mizer version 0.3
  if (!is.null(comment(params@pred_kernel))) {
    # n_eff_prey is the total prey abundance by size exposed to each
    # predator (prey not broken into species - here we are just working out
    # how much a predator eats - not which species are being eaten - that is
    # in the mortality calculation
    # \sum_j \theta_{ij} N_j(w_p) w_p dw_p
    n_eff_prey <- sweep(params@interaction %*% n, 2, 
                        params@w * params@dw, "*", check.margin = FALSE) 
    # pred_kernel is predator species x predator size x prey size
    # So multiply 3rd dimension of pred_kernel by the prey biomass density
    # Then sum over 3rd dimension to get consumption rate of each predator by 
    # predator size
    # This line is a bottle neck
    phi_prey_species <- rowSums(sweep(
      params@pred_kernel[, , idx_sp, drop = FALSE],
      c(1, 3), n_eff_prey, "*", check.margin = FALSE), dims = 2)
    # Eating the background
    # This line is a bottle neck
    phi_prey_background <- params@species_params$interaction_resource *
      rowSums(sweep(
        params@pred_kernel, 3, params@dw_full * params@w_full * n_pp,
        "*", check.margin = FALSE), dims = 2)
    encounter <- scalar * params@search_vol * (phi_prey_species + phi_prey_background)
  } else {
    prey <- outer(params@species_params$interaction_resource, n_pp)
    prey[, idx_sp] <- prey[, idx_sp] + params@interaction %*% n
    # The vector prey equals everything inside integral (3.4) except the feeding
    # kernel phi_i(w_p/w).
    prey <- sweep(prey, 2, params@w_full * params@dw_full, "*")
    # Eq (3.4) is then a convolution integral in terms of prey[w_p] and phi[w_p/w].
    # We approximate the integral by the trapezoidal method. Using the
    # convolution theorem we can evaluate the resulting sum via fast fourier
    # transform.
    # mvfft() does a Fourier transform of each column of its argument, but
    # we need the Fourier transforms of each row, so we need to apply mvfft()
    # to the transposed matrices and then transpose again at the end.
    avail_energy <- Re(base::t(mvfft(base::t(params@ft_pred_kernel_e) * 
                                       mvfft(base::t(prey)),
                                     inverse = TRUE))) / length(params@w_full)
    # Only keep the bit for fish sizes
    avail_energy <- avail_energy[, idx_sp, drop = FALSE]
    # Due to numerical errors we might get negative or very small entries that
    # should be 0
    avail_energy[avail_energy < 1e-18] <- 0
    
    encounter <- scalar * params@search_vol * avail_energy
  }
  
  # Add contributions from other components
  for (i in seq_along(params@other_encounter)) {
    encounter <- encounter + 
      do.call(params@other_encounter[[i]], 
              list(params = params,
                   n = n, n_pp = n_pp, n_other = n_other,
                   component = names(params@other_encounter)[[i]], ...))
  }
  return(encounter)
}


# mizerFeedingLevel with tempetature dependence on intake_max

mizerFeedingLevelTemp <- function(params, n, n_pp, n_other, t, encounter, ...) {
  # temperature component
  scalar <- tempFun(params@w, temperature = params@other_params$other$temperature[t+1],
                    t_d = params@species_params$t_d,
                    t_ref = params@species_params$t_ref,
                    Ea = params@species_params$Ea,
                    c_a = params@species_params$c_a,
                    Ed = params@species_params$Ed,
                    c_d = params@species_params$c_d)
  
  scalar <- t(drop(scalar))
  
  return(encounter / (encounter + params@intake_max * scalar))
}


# mizerPredRate with temperature dependence on search_volume

mizerPredRateTemp <- function(params, n, n_pp, n_other, t, feeding_level, ...) {
  no_sp <- dim(params@interaction)[1]
  no_w <- length(params@w)
  no_w_full <- length(params@w_full)
  
  # temperature component
  scalar <- tempFun(params@w, temperature = params@other_params$other$temperature[t+1],
                    t_d = params@species_params$t_d,
                    t_ref = params@species_params$t_ref,
                    Ea = params@species_params$Ea,
                    c_a = params@species_params$c_a,
                    Ed = params@species_params$Ed,
                    c_d = params@species_params$c_d)
  
  scalar <- t(drop(scalar))
  
  # If the the user has set a custom pred_kernel we can not use fft.
  # In this case we use the code from mizer version 0.3
  if (!is.null(comment(params@pred_kernel))) {
    n_total_in_size_bins <- sweep(n, 2, params@dw, '*', check.margin = FALSE)
    # The next line is a bottle neck
    pred_rate <- sweep(params@pred_kernel, c(1,2),
                       (1 - feeding_level) * params@search_vol * scalar *
                         n_total_in_size_bins,
                       "*", check.margin = FALSE)
    # integrate over all predator sizes
    pred_rate <- colSums(aperm(pred_rate, c(2, 1, 3)), dims = 1)
    return(pred_rate)
  }
  
  # Get indices of w_full that give w
  idx_sp <- (no_w_full - no_w + 1):no_w_full
  # We express the result as a a convolution  involving
  # two objects: Q[i,] and ft_pred_kernel_p[i,].
  # Here Q[i,] is all the integrand of (3.12) except the feeding kernel
  # and theta
  Q <- matrix(0, nrow = no_sp, ncol = no_w_full)
  # We fill the end of each row of Q with the proper values
  Q[, idx_sp] <- sweep( (1 - feeding_level) * params@search_vol * scalar * n, 2,
                        params@dw, "*")
  
  # We do our spectral integration in parallel over the different species
  pred_rate <- Re(t(mvfft(t(params@ft_pred_kernel_p) *
                            mvfft(t(Q)), inverse = TRUE))) / no_w_full
  # Due to numerical errors we might get negative or very small entries that
  # should be 0
  pred_rate[pred_rate < 1e-18] <- 0
  
  return(pred_rate * params@ft_mask)
}


tempFun <- function(w, temperature, t_d, t_ref, Ea, c_a = 0, Ed = 0, c_d = 0) # 
{
  # function needs to manage multiple formats but for now just doing one case
  # assuming that we get one value per species and arguments are vector

  
  scalarArray <- array(data = NA, dim = c(length(temperature), length(w), length(Ea)), 
                       dimnames = list("temperature" = temperature, "w" = w, "species" = seq(1:length(Ea))))
  
  
  # tempFun returns a matrix with w (size) as columns and temperature as rows
  k = 8.617332e-5 # Boltzmann constant 
  # equation
  # exp((-Ea/k)*((1/temperature) - (1/t_ref))) *
  # (1 + exp((-Ed/k)*((1/temperature) - (1/t_d))))^-1 *
  # (1 + exp((-Ed/k)*((1/Tref) - (1/t_d))))
  
  # converting to Kelvin from Celcius
  temperature <- temperature + 273 
  t_ref <- t_ref + 273
  t_d <- t_d + 273
  
  for(iSpecies in 1:length(Ea))
  {
    scalarArray[,,iSpecies] <- t(sapply(w,FUN = function(x){x^(c_a[iSpecies]*(temperature - t_ref[iSpecies]))}) * 
                                   exp(-Ea[iSpecies]/k*(1/temperature - 1/t_ref[iSpecies]))) *
      t(1/(sapply(w,FUN = function(x){x^(c_d[iSpecies]*(temperature - t_d[iSpecies]))}) *
             (exp(-Ed[iSpecies]/k*(1/temperature - 1/t_d[iSpecies])) + 1)))*
      sapply(w,FUN = function(x){x^(c_d[iSpecies]*(t_ref[iSpecies] - t_d[iSpecies]))}) *
      (exp(-Ed[iSpecies]/k*(1/t_ref[iSpecies] - 1/t_d[iSpecies])) + 1) 
  }
  return(scalarArray)
}



#other examples of temperature functions

### temperature functions ###

# expFun calculate the temperature scalar by size depending on temperature, activation energy (var1) and mass corrected temperature scaling (var2) using an exponential method
# Ea is activation energy of the rate we want to look at between intake/mortality/metabolism/maturation
# c_a is mass-correction of the temperature scalar in the rising part of the scalar. When 0, rates scale with temperatures equally for all size bins
# T_ref is the reference temperature (at which the temperature scalar = 1)
# c_d is mass-correction of the temperature scalar in the deactivation-part of the scalar. When 0, rates scale with temperatures equally for all size bins
# t_d is the temperature where deactivation starts
# object is mizer object with all necessary parameters (so we might get var1 to 3 in object directly)
# temperature is integer

# add parameters: Ea, Ed, c_a, t_ref, c_d, t_d * metabolism, maturation, mortality, intake

# a bit is missing in there, check papers folder in IMAS to get the formula from Anna

# tempFun <- function(temperature, t_ref, Ea, Ed = 0, c_a, c_d = 0, tmax, w) # default are 0 for now as deactivation is buggy
# {
#   # tempFun returns a matrix with w (size) as columns and temperature as rows
#   
#   # t_d is the beginning of the deactivation curve
#   if (Ed == 0) t_d <- temperature else t_d <- Ed * tmax / (Ed - tmax * 8.617332e-5 * log(Ed/Ea -1)) # so it does not crash as Ed > Ea to work
#   
#   # equation
#   # (w^(c_a*(temperature-t_ref)))  *exp((-Ea/8.617332e-5)*((1/temperature) - (1/t_ref)))
#   # *(1/(w^(c_d*(temperature-t_d)))*exp((-Ed/8.617332e-5)*((1/temperature) - (1/t_d))))
#   
#   temperatureScalar <- t(sapply(w,FUN = function(x){x^(c_a*(temperature-t_ref))}) *exp((-Ea/8.617332e-5)*((1/temperature) - (1/t_ref))) * (1/(sapply(w,FUN = function(x){x^(c_d*(temperature-t_d))}))*exp((-Ed/8.617332e-5)*((1/temperature) - (1/t_d)))))
#   
#   return(temperatureScalar)
# }



# that one is just activation part

# tempFun <- function(w, temperature, t_d = 25, t_ref, Ea, c_a, Ed = 0, c_d = 0) # default are 0 for now as deactivation is buggy
# {
#   # tempFun returns a matrix with w (size) as columns and temperature as rows
#   k = 8.617332e-5 # Boltzmann constant 
#   # equation
#   # (w^(c_a*(temperature-t_ref)))  *exp((-Ea/8.617332e-5)*((1/temperature) - (1/t_ref)))
#   # *(1/(w^(c_d*(temperature-t_d)))*exp((-Ed/8.617332e-5)*((1/temperature) - (1/t_d))))
#   
#   # converting to Kelvin from Celcius
#   temperature <- temperature + 273 
#   t_ref <- t_ref + 273
#   t_d <- t_d + 273
#   
#   
#   
#   temperatureScalar <- t(sapply(w,FUN = function(x){x^(c_a*(temperature-(t_ref)))}) *
#                            exp((-Ea/k)*((1/temperature) - (1/(t_ref)))))
#   
#   
#   return(temperatureScalar)
# }


# the one from my PhD


# tempFun <- function(w, temperature, t_d = 25, t_ref, Ea, c_a = 0, Ed = 0, c_d = 0) # 
# {
#   # tempFun returns a matrix with w (size) as columns and temperature as rows
#   k = 8.617332e-5 # Boltzmann constant 
#   # equation
#   # exp((-Ea/k)*((1/temperature) - (1/t_ref))) *
#   # (1 + exp((-Ed/k)*((1/temperature) - (1/t_d))))^-1 *
#   # (1 + exp((-Ed/k)*((1/Tref) - (1/t_d))))
#   
#   # converting to Kelvin from Celcius
#   temperature <- temperature + 273 
#   t_ref <- t_ref + 273
#   t_d <- t_d + 273
#   
#   temperatureScalar <- t(sapply(w,FUN = function(x){x^(c_a*(temperature - t_ref))}) * exp(-Ea/k*(1/temperature - 1/t_ref)   )) *
#     t(1/(sapply(w,FUN = function(x){x^(c_d*(temperature - t_d  ))}) *(exp(-Ed/k*(1/temperature - 1/t_d)) + 1)))*
#     sapply(w,FUN = function(x){x^(c_d*(t_ref       - t_d  ))}) *(exp(-Ed/k*(1/t_ref       - 1/t_d)) + 1) 
#   return(temperatureScalar)
# }
