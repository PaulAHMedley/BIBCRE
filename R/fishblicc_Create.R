##  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>
# ><> #  # ><> #  # ><> #  # ><> #  # ><> #
# ><> #  # ><> ####  HCR & POPULATION FUNCTIONS
# ><> #  # ><> #  # ><> #  # ><> #  # ><> #
##  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>


#' Creates a Beverton & Holt recruitment function
#'
#' A function is produced that given a SSB and recruitment process error
#' calculates the corresponding recruitment. The model allows for
#' autocorrelated errors.
#'
#' Expected recruitment R = R_0 * S / (R_b + S)
#'
#' The function produced takes a SSB vector of length equal to the number of
#' simulations and a lRd vector of 2 log recruitment deviates, the second being
#' the current log-deviate. This allows autocorrelation to be included.
#'
#' @param lR0  Log maximum recruitment
#' @param Rb   Second BH parameter defining an inflection point
#' @param lRs  Log-recruitment normal standard deviation
#' @param Rac  Log-recruitment autocorrelation coefficient
#' @return A function taking 2 parameters, SSB and lRd - a vector of 2 log
#' recruitment deviates
#' @export
#'
BH_recruitment <- function(lR0, Rb, lRs, Rac = 0) {
  force(lR0)
  lRs1 <- lRs * sqrt(1 - Rac * Rac)
  lRs2 <- lRs * Rac
  rm(Rac)
  if (lRs2 > 0) {
    return(function(SSB, lRd) {
      rec <- exp(lR0 + log(SSB) - log(Rb + SSB) +
                   lRd[, 2L] * lRs1 + lRd[, 1L] * lRs2)
      return(rec)
    })
  } else {
    return(function(SSB, lRd) {
      rec <- exp(lR0 + log(SSB) - log(Rb + SSB) +
                   lRd[, 2L] * lRs1)
      return(rec)
    })
  }
  
}

#' Creates a fixed recruitment function
#'
#' A function is produced that given recruitment process error
#' calculates the corresponding recruitment. The model allows for
#' autocorrelated errors. The SSB is ignored.
#' 
#' @inheritParams BH_recruitment
#' @export
#' 
fixed_recruitment <- function(lR0, lRs, Rac = 0) {
  force(lR0)
  lRs1 <- lRs * sqrt(1 - Rac * Rac)
  lRs2 <- lRs * Rac
  rm(Rac)
  if (lRs2 > 0) {
    return(function(SSB, lRd) {
      rec <- exp(lR0 + lRd[, 2L] * lRs1 + lRd[, 1L] * lRs2)
      return(rec)
    })
  } else {
    return(function(SSB, lRd) {
      rec <- exp(lR0 + lRd[, 2L] * lRs1)
      return(rec)
    })
  }
}  

#' Creates a function to estimate the mean length from length frequency data
#' matrices
#'
#' The function is based on the fishblicc data list: blicc_ld. The resulting
#' function takes a list of matrices with the same dimensions, one for each
#' gear, with the columns aligned to the length bins and rows the number of MCMC
#' simulations.
#'
#' @param  blicc_ld A list of settings used in fitted the length-based catch
#'   curve fishblicc model
#' @param  gears2use Either an integer vector of gears or "All"
#' @return A function taking 1 parameter, a matrix list of length frequency data
#' @export
#' 
meanLength_func <- function(blicc_ld, gears2use = "All") {
  xLMP <- with(blicc_ld, rep(LMP, each = NG))
  if (gears2use == "All")
    gears2use <- 1L:blicc_ld$NG
  return(function(lenfq) {
    lenfq <- lenfq[, , gears2use, drop = FALSE]
    dm <- dim(lenfq)
    ml <- vapply(seq_len(dm[2]), function(i) {
      sum(as.vector(lenfq[, i, ]) * xLMP) / sum(lenfq[, i, ])
    }, numeric(1))
    return(ml)
  })
}


#' Creates a function to estimate the SPR from length frequency data matrices
#'
#'
#' @inheritParams meanLength_func
#' @export
#'
SPR_func <- function(blicc_ld, gears2use = "All") {
  if (gears2use == "All")
    gears2use <- 1L:blicc_ld$NG
  blicc_ld$NG <- length(gears2use)
  blicc_ld$NQ <- blicc_ld$NG
  return(
    # Catch curve
    function(lenfq) {
      dm <- dim(lenfq)
      LN <- dm[1]
      NG_1 <- length(gears2use)
      SPR <- vapply(seq_len(dm[2]), function(i) {
        #blicc_ld$fq <- Map(\(x) x[, i, , drop = TRUE], lenfq)
        blicc_ld$fq <- split(lenfq[, i, gears2use], rep(NG_1, each = LN))
        res <- tryCatch(
          fishblicc::blicc_mpd(blicc_ld),
          error = function(e)
            NULL
        )
        if (is.null(res))
          return(NA_real_)
        val <- with(res, mpd[par == "SPR[1]"])
        if (length(val) == 0)
          return(NA_real_)  # handle no match/error
        return(as.numeric(val[1]))              # enforce numeric scalar
      }, numeric(1))
      return(SPR)
    }
  )
}

#' Creates a function to estimate the SPR from length frequency data matrices
#' 
#' This is currently too slow to use
#'
#' @inheritParams meanLength_func
#' @export
#'
SPR_func_para <- function(blicc_ld, gears2use = "All") {
  if (gears2use == "All")
    gears2use <- 1L:blicc_ld$NG
  future::plan(sequential)
  gc(verbose = FALSE)
  #blicc_ld$NG <- length(gears2use)
  #blicc_ld$NQ <- blicc_ld$NG
  return(function(lenfq) {
    tmp <- rowSums(lenfq)
    retain <- rep(TRUE, length(tmp))
    i <- 1L
    while (tmp[i] == 0) {
      retain[i] <- FALSE
      i <- i + 1L
    }
    if (i > 1L)
      retain[i - 1L] <- TRUE
    i <- length(tmp)
    while (tmp[i] == 0) {
      retain[i] <- FALSE
      i <- i - 1L
    }
    if (i < length(tmp))
      retain[i + 1L] <- TRUE
    
    blicc_ld$NB <- sum(retain)
    blicc_ld$LLB <- blicc_ld$LLB[retain]
    blicc_ld$LMP <- blicc_ld$LMP[retain]
    blicc_ld$wt_L <- blicc_ld$wt_L[retain]
    blicc_ld$ma_L <- blicc_ld$ma_L[retain]
    blicc_ld$M_L <- blicc_ld$M_L[retain]
    lenfq <- lenfq[retain, , , drop = FALSE]
    dm <- dim(lenfq)
    LN <- dm[1]
    NG_1 <- 1:dm[3]
    SPR <- future_sapply(seq_len(dm[2]), function(i) {
      blicc_ld$fq <- split(lenfq[, i, ], rep(NG_1, each = LN))
      res <- NULL
      tryCatch({
        res <- fishblicc::blicc_mpd(blicc_ld, refresh = 0)
      }, error = function(e)
        NULL)
      if (is.null(res))
        return(NA_real_)
      if (!tibble::is_tibble(res))
        return(NA_real_)
      val <- with(res, mpd[par == "SPR[1]"])
      if (length(val) == 0)
        return(NA_real_)
      as.numeric(val[1])
    }, future.seed = TRUE)
    return(SPR)
  })
}


#' Creates a function to estimate the SPR from length frequency data matrices
#'
#' Parallel version. Only estimates fishing mortality. Selectivity and other parameters are fixed to
#' the prior modes. Still too slow to use.
#'
#' @inheritParams meanLength_func
#' @export
#'
SPR_onlyF_para <- function(blicc_ld, gears2use = "All") {
  if (gears2use == "All")
    gears2use <- 1L:blicc_ld$NG
  future::plan(sequential)
  gc(verbose = FALSE)
  stmod <- rstan::stan_model(file = here("R", "BLICC_onlyF.stan"),
                             model_name = "blicc_onlyF")
  return(function(lenfq) {
    tmp <- rowSums(lenfq)
    retain <- rep(TRUE, length(tmp))
    i <- 1L
    while (tmp[i] == 0) {
      retain[i] <- FALSE
      i <- i + 1L
    }
    if (i > 1L)
      retain[i - 1L] <- TRUE
    i <- length(tmp)
    while (tmp[i] == 0) {
      retain[i] <- FALSE
      i <- i - 1L
    }
    if (i < length(tmp))
      retain[i + 1L] <- TRUE
    
    blicc_ld$NB <- sum(retain)
    blicc_ld$LLB <- blicc_ld$LLB[retain]
    blicc_ld$LMP <- blicc_ld$LMP[retain]
    blicc_ld$wt_L <- blicc_ld$wt_L[retain]
    blicc_ld$ma_L <- blicc_ld$ma_L[retain]
    blicc_ld$M_L <- blicc_ld$M_L[retain]
    lenfq <- lenfq[retain, , , drop = FALSE]
    dm <- dim(lenfq)
    LN <- dm[1]
    NG_1 <- 1:dm[3]
    SPR <- future_sapply(seq_len(dm[2]), function(i) {
      blicc_ld$fq <- split(lenfq[, i, ], rep(NG_1, each = LN))
      fit <- rstan::optimizing(
        stmod,
        data = blicc_ld,
        init = 0,
        hessian = TRUE,
        as_vector = FALSE,
        verbose = FALSE,
        iter = 10000,
        refresh = 0,
        tol_obj = 1e-12,
        tol_rel_obj = 1e3,
        tol_grad = 1e-8,
        tol_rel_grad = 1e5,
        tol_param = 1e-8
      )
      if (is.null(fit))
        return(NA_real_)
      val <- fit$par$SPR
      if (length(val) == 0)
        return(NA_real_)
      return(as.numeric(val[1]))
    }, future.seed = TRUE)
    return(SPR)
  })
}


#' Create the function to update the HCR index
#'
#' The HCR indicator function that is created calculates an indicator on
#' stock status that can be compared to reference points. The indicator
#' is either a function using length frequency. An optional moving average is
#' applied to the index.
#'
#' Either a mean length or length-based catch curve indicator is used.
#'
#' @param  ma  Moving average parameter for index time series
#' @param  index_func An index function that takes a single 3d length
#'   data matrices `lenfq` with dimensions length, sims, gears
#' @param  default_factor A default factor to allow recovery with no valid index
#'   (defaults to 1.05 = 5% increase)
#' @return A function taking the previous index and a 3d length frequency matrix that can be
#'   used to calculate the HCR index
#' @export
#'
fb_create_HCR_index_func <- function(ma, index_func, default_factor = 1.05) {
  return(function(indx_1, lenfq) {
    # lenfq is a 3d matrix of length frequency data: length, simulation, gear
    indx_0 <- index_func(lenfq)
    invalid <- is.na(indx_0) | (indx_0 <= 0)
    indx_0[invalid] <- indx_1[invalid] * default_factor
    return(ma * indx_0 + (1 - ma) * indx_1)
  })
}


#' Create a control function to be applied based on the index value above
#'
#' This function returns a function to calculate a control based on the index,
#' The only control currently supported
#' is an effort control. However, the function supports controls on multiple
#' gears.
#'
#' The trControl must be a matrix with a column for each gear
#'
#' @param refControl The F matrix with rows the number of simulations and
#'   columns the number of gears.
#' @param trIndex A single sorted vector of index inflection points from low to
#'   high. Same rows as `trControl`
#' @param trControl A matrix (or vector) of sorted vectors of control
#'   inflection points defining linear HCR. If there is no separate column for
#'   each gear, the vector is repeated.
#' @param change_limit An annual limit on the control changes as a proportion
#' @param control_type Must be "Effort". For future use.
#' @return A function to calculate the new controls as a matrix from the index
#'   and previous year's control. The matrix is number of simulations * the
#'   number of gears.
#' @export
#'
fb_create_control_func <- function(refControl,
                                   trIndex,
                                   trControl,
                                   change_limit = NA) {
  #if (control_type != "Effort") stop("Only effort control is currently supported.")
  trIndex <- c(0, trIndex, Inf)
  refControl <- as.array(refControl)
  trControl <- as.array(trControl)
  if (length(dim(trControl)) == 1L & length(dim(refControl) == 2L)) {
    trControl <- array(rep(trControl, dim(refControl)[2]), dim = c(length(trControl), dim(refControl)[2]))
  }
  if (dim(refControl)[2] != dim(trControl)[2])
    stop("trControl has the wrong dimensions (must be a vector or columns == ngears).")
  trControl <- rbind(trControl[1, ], trControl, trControl[nrow(trControl), ])
  if (is.na(change_limit))
    return(function(indx, prev_con) {
      # Calculate the control for each index value
      ii <- findInterval(indx, trIndex) #length nsim
      # Get nsim*ngear matrix of relative change
      rctrl <- trControl[ii, ] +
        (trControl[ii + 1L, ] - trControl[ii, ]) *
        (indx - trIndex[ii]) / (trIndex[ii + 1L] - trIndex[ii])
      ctrl <- rctrl * refControl
      return(ctrl)
    })
  else
    return(function(indx, prev_con) {
      # Calculate the control for each index value
      prop_change <- array(0, dim = dim(prev_con))
      ii <- findInterval(indx, trIndex)
      
      rctrl <- trControl[ii, ] +
        (trControl[ii + 1L, ] - trControl[ii, ]) *
        (indx - trIndex[ii]) / (trIndex[ii + 1L] - trIndex[ii])
      ctrl <- rctrl * refControl
      tochange <- prev_con > 0
      prop_change[tochange] <- (ctrl[tochange] - prev_con[tochange]) /
        prev_con[tochange]
      
      var_limit <- (abs(prop_change) > change_limit)
      if (any(var_limit)) {
        prop_change <- 1 + sign(prop_change[var_limit]) * change_limit
        ctrl[var_limit] <- prev_con[var_limit] * prop_change
      }
      return(ctrl)
    })
}


#' Extracts reference points from an HCR simulation
#'
#' The mean of all simulations is used where appropriate
#'
#' @param HCR_MSE HCR_MSE function created by [create_fishblicc_MSE]
#' @return List of MSY and MSY-proxy based reference points
#'   mortality for relative changes to that year
#' @export
#'
fb_get_refpt <- function(HCR_MSE) {
  ref_pt <- get("ref_pt", envir = environment(HCR_MSE))
  nsim <- get("nsim", envir = environment(HCR_MSE))
  ref_pt <- lapply(ref_pt, function(rp) {
    if (length(rp) == nsim)
      return(mean(rp))
    else
      return(rp)
  })
  return(ref_pt)
}


#' Estimates the MSY reference points from a simple age structured model at 
#' equilibrium
#'
#' @param Ngtg Number of growth groups
#' @param Nage Number of age groups
#' @param LN Number of length bins
#' @param yr_steps Number of steps within a year
#' @param ref_F Reference F for relative changes
#' @param selectivity Selectivity vector
#' @param mM Natural mortality
#' @param sim_idx Indices for each simulation for population vector
#' @param len_idx Indices for each length bin for population vector
#' @param wtmat Mature weight for each length bin
#' @param wt Weight for each length bin
#' @param R0 R0 for BH SR
#' @param Rb Rb for BH SR
#' @param R00 R00 for BH SR
#' @param SampleSize Effective sample size for length frequency data
#' @param IndexFunc Function to calculate HCR index
#' @param incrementAge logical indicating when the age is incremented with 
#'   yr_steps
#' @return List of reference points: SSB0, SSBMSY, index at MSY, FMSY and MSY.
#' @noRd
#'
fb_calc_MSY_refpt <- function(Ngtg,
                              Nage,
                              LN,
                              yr_steps,
                              ref_F,
                              selectivity,
                              mM,
                              sim_idx,
                              len_idx,
                              wtmat,
                              wt,
                              R0,
                              Rb,
                              R00,
                              SampleSize,
                              IndexFunc,
                              incrementAge
                              ) {
  nsim <- nrow(ref_F)
  NG <- ncol(ref_F)
  dt <- 1 / yr_steps
  msy_pop1 <- double(Ngtg * Nage)
  msy_age_idx <- rep(seq_len(Nage), each = Ngtg)         # ages
  #msy_gtg_idx <- rep(seq_len(Ngtg), Nage)             # growth groups
  Pi <- seq_len(Nage * Ngtg)
  msy_rec_idx <- Pi[msy_age_idx == 1L]
  msy_Age_from_idx <- Pi[msy_age_idx != Nage]
  msy_Age_to_idx <- Pi[msy_age_idx != 1L]
  msy_plusgp_idx <- Pi[msy_age_idx == Nage]
  rm(Pi)
  
  msy_wt <- msy_surv_idx <- list()
  SSB0 <- SSBMSY <- double(nsim)
  FMSY <- array(0, dim = dim(ref_F))
  MSY  <- IMSY <- double(nsim)
  LF <- array(0, dim = c(LN, nsim, NG))
  
  for (simi in seq_len(nsim)) {
    msy_wtmat <- wtmat[sim_idx == simi]
    for (mi in seq_len(yr_steps)) {
      msy_surv_idx[[mi]] <- len_idx[[mi]][sim_idx == simi]
      msy_wt[[mi]] <- wt[[mi]][sim_idx == simi]
    }
    msy_mM <- mM[(simi - 1L) * LN + (1:LN)]
    msy_FaL <- double(LN) # F at length
    for (gi in seq_len(NG)) {
      msy_FaL <- msy_FaL + selectivity[[gi]][, simi] * rep(ref_F[simi, gi] * dt, each = LN)
    }
    msy_R0 <- R0[simi]
    msy_Rb <- Rb[simi]
    
    res <- stats::optim(0.1, function(Fmul) {
      msy_ZaL <- Fmul * msy_FaL + msy_mM
      surv <- exp(-msy_ZaL)
      # Construct equilibrium population size / age structure
      msy_pop1 <- fb_build_equil(
        Ngtg,
        Nage,
        yr_steps,
        surv,
        msy_surv_idx,
        msy_rec_idx,
        msy_Age_to_idx,
        msy_Age_from_idx,
        msy_plusgp_idx,
        incrementAge
      )
      pB <- sum(msy_pop1 * msy_wtmat)
      recruits <- msy_R0 - msy_Rb / pB
      if (recruits <= 0)
        return(-recruits)
      msy_pop1 <- msy_pop1 * recruits
      recruits <- recruits / (sum(incrementAge)*Ngtg)
      # catch mortality proportion
      p <- (msy_ZaL - msy_mM) / msy_ZaL
      catch <- 0
      for (mi in seq_len(yr_steps)) {
        msy_pop2 <- msy_pop1 * surv[msy_surv_idx[[mi]]]
        mort <- msy_pop1 - msy_pop2
        Ca <- mort * p[msy_surv_idx[[mi]]] # catch numbers
        catch <- catch + sum(Ca * msy_wt[[mi]])
        if (incrementAge[mi]) {
          msy_pop1[msy_Age_to_idx] <- msy_pop2[msy_Age_from_idx]
          msy_pop1[msy_rec_idx] <- recruits 
          msy_pop1[msy_plusgp_idx] <- msy_pop2[msy_plusgp_idx] + msy_pop1[msy_plusgp_idx]
        } else {
          msy_pop1 <- msy_pop2
        }
      } #mi
      return(-catch)
    }, method = "Brent", lower = 0, upper = 10)
    
    surv <- exp(-msy_mM)
    msy_pop1 <- fb_build_equil(
      Ngtg,
      Nage,
      yr_steps,
      surv,
      msy_surv_idx,
      msy_rec_idx,
      msy_Age_to_idx,
      msy_Age_from_idx,
      msy_plusgp_idx,
      incrementAge
    )
    SSB0[simi] <- sum(msy_pop1 * msy_wtmat) * R00[simi]  # unexploited SSB
    
    msy_ZaL <- res$par * msy_FaL + msy_mM
    surv <- exp(-msy_ZaL)
    # build equilibrium population size / age structure
    msy_pop1 <- fb_build_equil(
      Ngtg,
      Nage,
      yr_steps,
      surv,
      msy_surv_idx,
      msy_rec_idx,
      msy_Age_to_idx,
      msy_Age_from_idx,
      msy_plusgp_idx,
      incrementAge
    )
    
    msy_len_idx <- msy_surv_idx # same thing for one simulation
    pB <- sum(msy_pop1 * msy_wtmat)
    recruits <- msy_R0 - msy_Rb / pB
    msy_pop1 <- msy_pop1 * recruits
    recruits <- recruits / (sum(incrementAge)*Ngtg)
    # Get catch
    freq <- array(0, dim = c(LN, NG))
    for (mi in seq_len(yr_steps)) {
      idx <- sort(unique(msy_len_idx[[mi]]))
      msy_pop2 <- msy_pop1 * surv[msy_surv_idx[[mi]]]
      mort <- msy_pop1 - msy_pop2
      for (gi in seq_len(NG)) {
        pg <- (selectivity[[gi]][, simi] * ref_F[simi, gi] * dt) / msy_ZaL
        Ca <- mort * pg[msy_surv_idx[[mi]]] # catch numbers
        freq[idx, gi] <- freq[idx, gi] +
          tapply(Ca, INDEX = list(msy_len_idx[[mi]]), FUN = "sum")
      }
      if (incrementAge[mi]) {
        msy_pop1[msy_Age_to_idx] <- msy_pop2[msy_Age_from_idx]
        msy_pop1[msy_rec_idx] <- recruits 
        msy_pop1[msy_plusgp_idx] <- msy_pop2[msy_plusgp_idx] + msy_pop1[msy_plusgp_idx]
      } else {
        msy_pop1 <- msy_pop2
      }
    } #mi
    
    SSBMSY[simi] <- sum(msy_pop1 * msy_wtmat)
    FMSY[simi, ] <- ref_F[simi, ] * res$par
    MSY[simi] <- -res$value
    LF[, simi, ] <- round(freq * SampleSize / sum(freq))
  } #simi
  IMSY <- IndexFunc(LF)
  return(list(
    SSB0 = SSB0,
    SSBMSY = SSBMSY,
    IMSY = IMSY,
    FMSY = FMSY,
    MSY = MSY
  ))
}



#' Calculates reference points for a model based on equilibrium
#'
#'
#' @inheritParams fb_calc_MSY_refpt
#' @param recalc_F Whether to recalculate the fshing mortality reference 
#'   points or not
#' @param SSBtar  Target SSB which depends on the type of reference point being 
#'   used
#' @noRd
fb_find_refpts <- function(Ngtg,
                        Nage,
                        LN,
                        yr_steps,
                        ref_F,
                        selectivity,
                        mM,
                        sim_idx,
                        len_idx,
                        wtmat,
                        wt,
                        R0,
                        Rb,
                        SampleSize,
                        IndexFunc,
                        incrementAge,
                        recalc_F,
                        SSBtar) {
  nsim <- nrow(ref_F)
  NG <- ncol(ref_F)
  dt <- 1 / yr_steps
  tar_pop1 <- double(Ngtg * Nage)
  tar_age_idx <- rep(seq_len(Nage), each = Ngtg)       # ages
  Pi <- seq_len(Nage * Ngtg)
  tar_rec_idx <- Pi[tar_age_idx == 1L]
  tar_Age_from_idx <- Pi[tar_age_idx != Nage]
  tar_Age_to_idx <- Pi[tar_age_idx != 1L]
  tar_plusgp_idx <- Pi[tar_age_idx == Nage]
  rm(Pi)
  if (!recalc_F)
    solveF <- list(root = 1)
  tar_wt <- tar_surv_idx <- list()
  F_tar <- array(0, dim = dim(ref_F))
  Ca_tar  <- I_tar <- double(nsim)
  LF <- array(0, dim = c(LN, nsim, NG))
  invalid <- logical(nsim)
  for (simi in seq_len(nsim)) {
    tar_wtmat <- wtmat[sim_idx == simi]
    for (mi in seq_len(yr_steps)) {
      tar_surv_idx[[mi]] <- len_idx[[mi]][sim_idx == simi]
      tar_wt[[mi]] <- wt[[mi]][sim_idx == simi]
    }
    tar_mM <- mM[(simi - 1L) * LN + (1:LN)]
    tar_FaL <- double(LN) # F at length
    for (gi in seq_len(NG)) {
      tar_FaL <- tar_FaL + selectivity[[gi]][, simi] * rep(ref_F[simi, gi] * dt, each = LN)
    }
    tar_R0 <- R0[simi]
    tar_Rb <- Rb[simi]
    if (recalc_F) {

      solveF <- tryCatch(stats::uniroot(function(Fmul) {
        tar_ZaL <- Fmul * tar_FaL + tar_mM
        surv <- exp(-tar_ZaL)
        # Construct equilibrium population size / age structure
        tar_pop1 <- fb_build_equil(
          Ngtg,
          Nage,
          yr_steps,
          surv,
          tar_surv_idx,
          tar_rec_idx,
          tar_Age_to_idx,
          tar_Age_from_idx,
          tar_plusgp_idx,
          incrementAge
        )
        
        pSSB <- sum(tar_pop1 * tar_wtmat)
        Rec <- tar_R0 - tar_Rb / pSSB
        SSB <- pSSB * Rec
        return(SSB - SSBtar[simi])
      }, lower = 0, upper = 200),
      error = function(e) list(root = -1))
    }
    
    if (solveF$root < 0) {
      invalid[simi] <- TRUE
      solveF$root <- 1
    }
    
    tar_ZaL <- solveF$root * tar_FaL + tar_mM
    surv <- exp(-tar_ZaL)
    # build equilibrium population size / age structure
    tar_pop1 <- fb_build_equil(
      Ngtg,
      Nage,
      yr_steps,
      surv,
      tar_surv_idx,
      tar_rec_idx,
      tar_Age_to_idx,
      tar_Age_from_idx,
      tar_plusgp_idx,
      incrementAge
    )
    
    tar_len_idx <- tar_surv_idx # same thing for one simulation
    pSSB <- sum(tar_pop1 * tar_wtmat)
    Rec <- tar_R0 - tar_Rb / pSSB
    tar_pop1 <- tar_pop1 * Rec
    # Get catch
    freq <- array(0, dim = c(LN, NG))
    for (mi in seq_len(yr_steps)) {
      idx <- sort(unique(tar_len_idx[[mi]]))
      tar_pop2 <- tar_pop1 * surv[tar_surv_idx[[mi]]]
      mort <- tar_pop1 - tar_pop2
      for (gi in seq_len(NG)) {
        pg <- (selectivity[[gi]][, simi] * ref_F[simi, gi] * dt) / tar_ZaL
        Ca <- mort * pg[tar_surv_idx[[mi]]] # catch numbers
        freq[idx, gi] <- freq[idx, gi] +
          tapply(Ca, INDEX = list(tar_len_idx[[mi]]), FUN = "sum")
        Ca_tar[simi] <- Ca_tar[simi] + sum(Ca * tar_wt[[mi]])
      }
      tar_pop1 <- tar_pop2
    } #mi
    F_tar[simi, ] <- ref_F[simi, ] * solveF$root
    LF[, simi, ] <- round(freq * SampleSize / sum(freq))
  } #simi

  I_tar <- IndexFunc(LF)
  return(list(
    I_tar = I_tar,
    F_tar = F_tar,
    Ca_tar = Ca_tar,
    invalid = invalid
  ))
}

#' Set the reference point type in the HCR function
#' 
#' Reference point types can be "SSB0" (relative to the unexploited state), 
#' "MSY" (relative to SSB at MSY) or "Current" (relative to the values estimated
#' in the stock assessment). MSY is generally not recommended in data poor 
#' situations as it will not be well estimated and can lead to risky advice. 
#' SSB0 with the target set at precautionary values is probably a better proxy
#' for MSY. The "Current" setting may be useful to guide HCR that improve on the
#' current status where there is little confidence in what that status is.
#' 
#' @param HCR_MSE HCR_MSE function created by [create_fishblicc_MSE]
#' @param rp_type Reference point type ("SSB0", "MSY", "Current")
#' @param TRP Target reference point relative to rp_type
#' @param LRP Limit reference point relative to rp_type
#' @export
#' 
fb_set_ref_pt_type <- function(HCR_MSE, rp_type, TRP, LRP = NULL) {
  if ((length(rp_type) != 1L) | !is.character(rp_type) |
      !(rp_type %in% c("SSB0", "MSY", "Current"))) {
    stop("rp_type must be either 'SSB0', 'MSY' or 'Current'\n")
  }
  if (!((length(TRP) == 1L) & is.numeric(TRP) & (TRP > 0)))
    stop(
      "TRP must be a single positive number of the proportion of the type
         (Current, MSY or SSB0) to target.\n"
    )
  if ((rp_type == "SSB0") & (TRP >= 1.0))
    stop("If rp_type is SSB0, TRP must be [0.0 < TRP < 1.0] .\n")
  
  new_ref_pt <- tryCatch(
    get("ref_pt", pos = environment(HCR_MSE)),
    error = function(e) {
      stop("HCR_MSE must be a function created by a 'create_fishblicc_MSE'.\n")
    }
  )
  
  if (is.null(LRP)) LRP <- new_ref_pt$LRP else new_ref_pt$LRP <- LRP
  
  Ngtg <- as.integer(get("Ngtg", pos = environment(HCR_MSE)))
  Nage <- get("Nage", pos = environment(HCR_MSE))
  LN <- get("LN", pos = environment(HCR_MSE))
  yr_steps <- get("yr_steps", pos = environment(HCR_MSE))
  selectivity <- get("selectivity", pos = environment(HCR_MSE))
  mM <- get("mM", pos = environment(HCR_MSE))
  sim_idx <- get("sim_idx", pos = environment(HCR_MSE))
  len_idx <- get("len_idx", pos = environment(HCR_MSE))
  wtmat <- get("wtmat", pos = environment(HCR_MSE))
  wt <- get("wt", pos = environment(HCR_MSE))
  R0 <- get("R0", pos = environment(HCR_MSE))
  Rb  <- get("Rb", pos = environment(HCR_MSE))
  SampleSize <- get("SampleSize", pos = environment(HCR_MSE))
  IndexFunc <- get("index_func", pos = environment(HCR_MSE))
  incrementAge <- get("incrementAge", pos = environment(HCR_MSE))
  
  ref_F <- switch(
    rp_type,
    Current = get("mF", pos = environment(HCR_MSE))[[1]],
    MSY = new_ref_pt$FMSY,
    SSB0 = get("mF", pos = environment(HCR_MSE))[[1]]
  )
  
  SSBtar <- switch(
    rp_type,
    Current = get("SSB", pos = environment(HCR_MSE))[[1]] * TRP,
    MSY = new_ref_pt$SSBMSY * TRP,
    SSB0 = new_ref_pt$SSB0 * TRP
  )

  rp <- fb_find_refpts(
    Ngtg,
    Nage,
    LN,
    yr_steps,
    ref_F,
    selectivity,
    mM,
    sim_idx,
    len_idx,
    wtmat,
    wt,
    R0,
    Rb,
    SampleSize = SampleSize,
    IndexFunc,
    incrementAge,
    recalc_F = TRUE,
    SSBtar
  )
  
  new_ref_pt$rp_type <- rp_type
  new_ref_pt$TRP <- TRP
  new_ref_pt$B_tar <- SSBtar
  new_ref_pt$B_lim <- SSBtar * LRP
  new_ref_pt$I_tar <- rp$I_tar
  new_ref_pt$F_tar <- rp$F_tar
  new_ref_pt$Ca_tar <- rp$Ca_tar
  assign("ref_pt", new_ref_pt, pos = environment(HCR_MSE))
  return(HCR_MSE)
}


#' Sets the recruitment deviates type, whether shared between simulations or 
#' independent
#' 
#' @inheritParams fb_set_ref_pt_type
#' @param sim_recruit Whether simulated recruitment deviates are independent 
#'   for each simulation or shared across simulations
#' @export
#'
set_sim_recruit <- function(HCR_MSE, sim_recruit) {
  if (!(sim_recruit %in% c("shared", "independent")))
    stop("Error: sim_recruit can only be 'shared' or 'independent'.")
  assign("sim_recruit", sim_recruit, pos = environment(HCR_MSE))
  return(HCR_MSE)
}



#' Builds population structure assuming equilibrium
#' 
#' @inheritParams fb_calc_MSY_refpt
#' @param eq_surv  Survival probability for each length
#' @param eq_surv_idx Survival allocation to the population vector
#' @param rec_idx Recruitment index for the population vector
#' @param age_to_idx Cohort transition (growth) to the population vector
#' @param age_from_idx Cohort transition from the population vector
#' @param plusgp_idx Age plus group index in the population vector
#' @noRd
#'
fb_build_equil <- function(Ngtg,
                           Nage,
                           yr_steps,
                           eq_surv,
                           eq_surv_idx,
                           rec_idx,
                           age_to_idx,
                           age_from_idx,
                           plusgp_idx,
                           incrementAge) {
  recruits <- rep(1 / (sum(incrementAge) * Ngtg), Ngtg)
  eq_pop1 <- double(length(eq_surv_idx[[1]]))
  eq_pop1[rec_idx] <- recruits # add recruits as split gtg age 1 for each sim
  mi <- 1L
  for (ai in seq_len(Nage + 1L)) {
    for (mi in seq_len(yr_steps)) {
      eq_pop1 <- eq_pop1 * eq_surv[eq_surv_idx[[mi]]]
      if (incrementAge[mi]) {
        eq_pop2 <- eq_pop1
        eq_pop1[age_to_idx] <- eq_pop2[age_from_idx]
        eq_pop1[rec_idx] <- recruits 
        eq_pop1[plusgp_idx] <- eq_pop2[plusgp_idx] + eq_pop1[plusgp_idx]
      }
    }
  }
  plusgp_n1 <- 0
  toler <- 0.0001
  repeat {
    plusgp_n2 <- plusgp_n1
    plusgp_n1 <- sum(eq_pop1[plusgp_idx])
    for (mi in seq_len(yr_steps)) {
      eq_pop1 <- eq_pop1 * eq_surv[eq_surv_idx[[mi]]]
      if (incrementAge[mi]) {
        eq_pop2 <- eq_pop1
        eq_pop1[age_to_idx] <- eq_pop2[age_from_idx]
        eq_pop1[rec_idx] <- recruits # add recruits as split gtg age 1 for each sim
        eq_pop1[plusgp_idx] <- eq_pop2[plusgp_idx] + eq_pop1[plusgp_idx]
      }
    }
    if (abs(plusgp_n2 - plusgp_n1) < toler)
      break
  }
  return(eq_pop1)
}


#' Create a HCR MSE function based on a fishblicc fit.
#'
#' The HCR function returns simulation results from HCR parameters using a
#' fishblicc model fit results.
#'
#' The the same state space lognormal random number can be applied across
#' simulations, so HCR performance is more comparible.
#' 
#' Uses fishblicc results from the `blicc_ref_pts` function for growth and 
#' mortality parameters, but rate parameter K, and recruitment parameters must 
#' be provided to convert fishblicc from a static model to a time series
#' model. The MCMC draws are used to provide a measure of uncertainty from
#' the posterior fit. 
#'
#' @details 
#' 
#' The model is a length and age structured simulation of a single species. 
#' Growth variability is simulated using 50 growth-type-groups over a 
#' Gamma-distributed Linf parameter based on the fishblicc fit. The supplied 
#' von Bertalanffy growth rate (K) can also be randomly drawn. Most parameters
#' apart from K are taken from the fishblicc fit. 
#' 
#' The stock recruitment relationship is the Beverton and Holt function, with 
#' steepness set by the user. A default precautionary steepness = 0.75 is 
#' applied This is raised if the population dynamics do not support this 
#' steepness level and a warning is given. This can occur when the estimated 
#' SPR is very low. The unexploited recruitment is estimated as scaling factor 
#' consistent with the catch (catch_weight).
#' 
#' Reference points are estimated. The default reference point is 40% 
#' unexploited SSB, which is precautionary and consistent with 40% SPR. However, 
#' SPR does not take account of a stock recruitment relationship which is 
#' included in this model. MSY reference points are also provided for reference. 
#' These can be used for the target, but in general this is not recommended in
#' data-poor situations as estimated values depend entirely on the assumed 
#' stock-recruitment and are generally not precautionary. You could treat MSY as
#' a limit reference point and set the target as two times MSY, which could work 
#' well in some circumstances (see [fb_set_ref_pt_type]). 
#' 
#' The HCR function that is produced takes six parameters. 
#'   - trIndex: a vector of HCR index inflection points
#'   - trControl: a single vector or list of vectors, one for each gear
#'   - change_limit: an annual proportional change limit on the control, if any  
#'   - ma: a moving average parameter for the index (default=0.5)
#'   - control_type: ignored (default="Effort")
#'   - ctrl_pF: Proportion of the fishing mortality affected by the control. So 
#'     (1-ctrl_pF)*F is fixed and unaffected. This could represent uncontrolled
#'     subsistence fishing for example (default=1). 
#' 
#' @param fishblicc_fit A list returned from package fishblicc `blicc_ref_pts`
#'   function (required)
#' @param vbK The von Bertalanfy growth rate parameter K (/year). Can be a 
#'   single value, or a vector the same length as nsim to allow for random 
#'   variation. (required)
#' @param catch_weight The total catch weight for the year fishblicc was fit 
#'   (last data year). Used for scaling values so catch values align in 
#'   magnitude. (default=1000)
#' @param control_type Type of control applied by HCR. Only "Effort" is 
#'   currently supported for this model.
#' @param index_type HCR index used (default="MeanLength")
#' @param nsim Number of simulations to run (default=1000)
#' @param proj_length  Projection length in years (default=50)
#' @param StartYear Year projections start, used for plot labels only (default=1).
#' @param AnnualRecruitment = TRUE,
#' @param yr_steps Time steps within each year. The higher number is more 
#'   accurate but slower simulations. Can be left to default so each step is 
#'   equivalent to max(2, vbK=0.2) (default=0, i.e. estimate). 
#' @param AssessmentPeriod Period within the year when the HCR is calculated for 
#'   the following year's control (default=yr_steps).
#' @param SampleSize The length frequency sample size. If not provided, the 
#'   fishblicc sample size is used  (default=0, i.e. estimate). 
#' @param SRDelay annual delay between SSB calculated subsequent recruitment 
#'   (default=1).
#' @param h_steepness Beverton and Holt stock recruitment relationship steepness 
#'   parameter (default=0.75).
#' @param autocorr Autocorrelation (0-1.0) for recruitment deviates (default=0).
#' @param recruit_cv Lognormal sigma (coefficient of variation) for recruitment 
#'   deviates (default=0.3)
#' @param sim_recruit recruitment deviates can be 'shared' across the simulated
#'   parameters set or generated 'independent'-ly for each parameter set.
#' @param ref_control A reference control (fishing mortality) used to multiply
#'   the HCR control. This is usually the current fishing mortality, so HCR 
#'   changes are relative to the current fishery (default = NULL i.e. current 
#'   fishing mortality estimated in fishblicc).
#' @param rp_type = Either "SSB0", "MSY" or current (default = "SSB0"). 
#'   See details.
#' @param rseed Random seed if needed for repeated simulations (default = Sys.time()). 
#'   This random seed is used for each HCR tested using the returned 
#'   function.
#' @return A function that will apply the HCR using the fishblicc model fit and
#'   specified HCR parameters: index-control inflection point vectors, the
#'   control change limit and moving average parameter. See details.
#' @export
#'
create_fishblicc_MSE <- function(fishblicc_fit,
                                 vbK,
                                 catch_weight = 1000,
                                 control_type = "Effort",
                                 index_type = "MeanLength",
                                 nsim = 1000,
                                 proj_length = 50,
                                 StartYear = 1,
                                 AnnualRecruitment = TRUE,
                                 yr_steps = 0,
                                 AssessmentPeriod = yr_steps,
                                 SampleSize = 0,
                                 SRDelay = 1,
                                 h_steepness = 0.75,
                                 autocorr = 0,
                                 recruit_cv = 0.3,
                                 sim_recruit = "independent",
                                 ref_control = NULL,
                                 rp_type = "SSB0",
                                 rseed = Sys.time()) {
  stock_assessment <- "fishblicc" 
  cat("Extracting data.\n")
  # Check fishbicc_fit is a reference point fit list and recruitment is a function or NULL
  RNGkind("L'Ecuyer-CMRG")
  set.seed(rseed) # same seed for all HCR
  factory_stream <- parallel::nextRNGStream(.Random.seed)
  sim_stream     <- parallel::nextRNGStream(factory_stream)
  .Random.seed <- factory_stream
  
  # Keep info needed from the fit
  if (fishblicc_fit$ld$NT == 1) {
    blicc_ld <- fishblicc_fit$ld
    dr_df <- fishblicc_fit$dr_df
  } else {
    # Use only last time period (change to scenario time period)
    blicc_ld <- fishblicc_fit$scenario$time_period_ld
    dr_df <- fishblicc_fit$dr_df |>
      mutate(Fk = purrr::map(Fk, \(x) x[fishblicc_fit$ld$NT]),
             SPR = purrr::map(SPR, \(x) x[fishblicc_fit$ld$NT]))
  }
  rm(fishblicc_fit)

  # Extract sufficient parameters from the fit for the n simulations
  if (nsim > nrow(dr_df)) {
    nsim <- nrow(dr_df)
  } else if (nsim < nrow(dr_df)) {
    # Select a random sample of iterations for the simulations
    dr_df <- dr_df |>
      dplyr::slice_sample(n = nsim, replace = FALSE) |>
      dplyr::arrange(.draw)
  }
  
  # Get time step
  if (yr_steps==0) {
    yr_steps <- ceiling(5*mean(vbK))    # steps are for K==0.2
    AssessmentPeriod <- yr_steps
  }
  # cohort / recruitment periods when age is incremented
  if (AnnualRecruitment) {
    incrementAge <- c(logical(yr_steps-1L), TRUE)
  } else {
    incrementAge <- rep(TRUE, yr_steps)
  }

  # expand vbK parameter
  if (length(vbK) == 1) {
    vbK <- rep(vbK, nsim)
  } else {
    if (length(vbK) > nsim) {
      vbK <- sample(vbk, size = nsim, replace = FALSE)
      warning(paste(
        "vbK has been randomly sampled to length",
        as.character(nsim)
      ))
    } else if (length(vbK) != nsim)
      stop(paste("vbK must be of length 1 or the same as nsim", as.character(nsim)))
  }
  
  # extract parameters from the fit
  Linf <- dplyr::pull(dr_df, Linf)
  Galpha <- dplyr::pull(dr_df, Galpha)
  Gbeta <- dplyr::pull(dr_df, Gbeta)
  NBphi <- dplyr::pull(dr_df, NB_phi)
  NG <- blicc_ld$NG
  Fk <- matrix(unlist(dr_df$Fk), ncol = NG, byrow = TRUE) # nsim * NG matrix
  mF <- list()
  mF[[1]] <- sweep(Fk,
                   MARGIN = 1,
                   STATS = vbK,
                   FUN = "*") # remove K
  rm(Fk)
  
  # Ages
  dAge <- 1 / sum(incrementAge)
  
  FirstAge <- min(-log(1-blicc_ld$LLB[1]/Linf)/vbK)
  LastAge <- dAge+(-log(1 - 0.975) / min(vbK))
  age <- seq(from = FirstAge, to = LastAge, by = dAge)
  Nage <- length(age)
  dt <- 1 / yr_steps 
  
  cat("Setting MSE dimensions and initial structure.\n")
 
  # Length
  Ngtg <- 50L
  
  if (SampleSize <= 0)
    SampleSize <- sum(unlist(blicc_ld$fq)) # D$samplesize
  PYN <- proj_length
  PN <- PYN * yr_steps
  
  if (is.null(ref_control)) {
    ref_control <- mF[[1]]
  } else {
    if (!is.numeric(ref_control) |
        any(dim(ref_control) != c(nsim, NG))) {
      stop(
        "parameter ref_control must be numeric vector of fishing mortality
           or catches the same length as the number of gears.\n"
      )
    }
  }
  
  # indexes
  age_idx <- rep(seq_len(Nage), each = nsim * Ngtg)         # ages
  sim_idx <- rep(rep(seq_len(nsim), each = Ngtg), Nage)   # simulations
  gtg_idx <- rep(seq_len(Ngtg), Nage * nsim)              # growth groups
  Pi <- seq_len(Nage * nsim * Ngtg)
  rec_idx <- Pi[age_idx == 1L]
  Age_from_idx <- Pi[age_idx != Nage]
  Age_to_idx <- Pi[age_idx != 1L]
  plusgp_idx <- Pi[age_idx == Nage]
  rm(Pi)
  
  # Growth groups lengths
  quartiles <- ((1:Ngtg) - 0.5) / Ngtg
  
  # Weights and weight*maturity is calculated from the model parameters in blicc_ld
  wtmat <- wt <- len <- list()
  age_inc <- 0.5 / yr_steps
  len_val <- wt_val <- rep(list(NULL), yr_steps)
  ii <- 1L
  for (i in seq_len(yr_steps)) {
    if (i>1L) {
      if (incrementAge[i-1L]) {
        age_inc <- 0.5 / yr_steps  # reset age when new recruits
        ii <- 1L
      } else {
        age_inc <- age_inc + dAge
        ii <- ii + 1L
      }
    }
    if (is.null(len_val[[ii]])) {
      len_val[[ii]] <- vapply(
        seq_len(Ngtg * nsim * Nage),
        FUN = function(j)
          qgamma(quartiles[gtg_idx[j]], shape = Galpha[sim_idx[j]], 
                 rate = Gbeta[sim_idx[j]]) *
                  (1 - exp(-vbK[sim_idx[j]] * (age[age_idx[j]] + age_inc))),
        FUN.VALUE = double(1)
      )
      wt_val[[ii]] <- blicc_ld$a * len_val[[ii]]^blicc_ld$b
    }
    len[[i]] <- len_val[[ii]]
    wt[[i]] <- wt_val[[ii]]
  }
  rm(len_val, wt_val)
  # SSB only calculated at the end of the year
  wtmat <- wt[[1]] / (1 + exp(-blicc_ld$Ls * (len[[1]] - blicc_ld$L50))) 
  rm(quartiles)
  
  LB <- blicc_ld$LLB
  LMP <- blicc_ld$LMP
  M_L <- blicc_ld$M_L
  # Increase length bins if necessary
  if (max(unlist(len)) > max(LB)) {
    lenmax <- max(unlist(len))
    LBmax <- max(LB)
    inc <- max(1, LB[length(LB)] - LB[length(LB) - 1L])
    N <- 0L
    repeat {
      LBmax <- LBmax + inc
      N <- N + 1L
      if (LBmax > lenmax)
        break
    }
    newvalues_idx <- with(blicc_ld, (NB + 1L):(NB + N))
    LB <- c(LB, seq(
      from = LB[length(LB)] + inc,
      to = LBmax,
      by = inc
    ))
    LMP <- c(LMP, seq(
      from = LMP[length(LMP)] + inc,
      to = LB[length(LB)] + 1,
      by = inc
    ))
    if (length(LMP) > length(LB))
      stop("Error: LMP and LB incompatible lengths.\n")
    LN <- length(LB)
    M_L <- blicc_ld$M_L
    WMP <- c(blicc_ld$wt_L, blicc_ld$a * LMP[newvalues_idx]^blicc_ld$b)  #Weight*Maturity for each length
    WtMat <- c(blicc_ld$ma_L, WMP[newvalues_idx] * rep(blicc_ld$ma_L[blicc_ld$NB] /
                                                         WMP[blicc_ld$NB], N))  #Weight*Maturity for each length
    ILen <- 1 / blicc_ld$LMP
    M_Lmod <- stats::lm(M_L ~ ILen)
    M_Lest <- stats::predict(M_Lmod, newdata = data.frame(ILen = LMP[newvalues_idx]))
    # Update blicc_ld
    cat(paste0(
      "Length bins from ",
      as.character(blicc_ld$NB),
      " to ",
      as.character(LN),
      "\n"
    ))
    blicc_ld$NB <- LN
    blicc_ld$LLB <- LB
    blicc_ld$LMP <- LMP
    blicc_ld$wt_L <- WMP
    blicc_ld$ma_L <- WtMat #Weight*Maturity for each length
    blicc_ld$M_L <- c(M_L, M_Lest)
    blicc_ld$fq <- lapply(blicc_ld$fq, \(freq) return(c(freq, rep(0, N))))
    rm(inc,
       M_Lest,
       M_Lmod,
       ILen,
       M_L,
       WMP,
       WtMat,
       LBmax,
       N,
       lenmax,
       newvalues_idx)
  }
  
  surv_idx <- len_idx <- list()
  for (i in 1:yr_steps) {
    len_idx[[i]] <- findInterval(len[[i]], vec = LB)  # Find length bin
    len_idx[[i]][len_idx[[i]] == 0] <- 1L
    surv_idx[[i]] <- (sim_idx - 1L) * LN + len_idx[[i]] # survival calculation index
  }
  #rm(len)
  
  #CL_N <- PL_N <- array(0, c(nsim, LN))
  SSB <- Rec <- pjIndex <- LF <- list() #Age length keys for each time step for males and females
  
  # Mortality and selectivity
  mM <- as.vector(outer(blicc_ld$M_L, dr_df$Mk * vbK * dt, "*")) # LN * nsim matrix
  
  # flatten selectivity
  selectivity <- rep(list(matrix(0, nrow = LN, ncol = nsim)), NG)
  for (rowi in 1:nrow(dr_df)) {
    sl <- fishblicc::Rselectivities(unlist(dr_df$Sm[rowi]), blicc_ld)
    for (gi in 1:NG) {
      selectivity[[gi]][, rowi] <- sl[[gi]]
    }
  }
  rm(sl)
  
  # Initial mortality
  FaL <- list() # F at length
  ZaL <- mM
  for (gi in seq_len(NG)) {
    FaL[[gi]] <- selectivity[[gi]] * rep(mF[[1]][, gi] * dt, each = LN)
    ZaL <- ZaL + FaL[[gi]]
  }
  surv <- exp(-ZaL)
  
  # Construct equilibrium initial population size / age structure
  pop1 <- fb_build_equil(Ngtg,
                         Nage,
                         yr_steps,
                         surv,
                         surv_idx,
                         rec_idx,
                         Age_to_idx,
                         Age_from_idx,
                         plusgp_idx,
                         incrementAge)
  # Get catch
  recruits <- rep(1 / (sum(incrementAge)*Ngtg), each = Ngtg) # add recruits as split gtg age 1 for each sim
  p <- (ZaL - mM) / ZaL
  catch <- double(nsim)
  for (mi in seq_len(yr_steps)) {
    pop2 <- pop1 * surv[surv_idx[[mi]]]
    mort <- pop1 - pop2
    Ca <- mort * p[surv_idx[[mi]]] # catch numbers
    catch <- catch + tapply(Ca * wt[[mi]], INDEX = sim_idx, FUN = "sum")
    if (incrementAge[mi]) {
      pop1[Age_to_idx] <- pop2[Age_from_idx]
      pop1[rec_idx] <- recruits
      pop1[plusgp_idx] <- pop2[plusgp_idx] + pop1[plusgp_idx]
    } else {
      pop1 <- pop2
    }
  } #mi
  # Approximate rescale population size based on estimated catch
  recruit <- catch_weight / catch   # separate recruit for each sim
  pop1 <- pop1 * recruit[sim_idx]
  xSSB <- tapply(pop1 * wtmat, INDEX = sim_idx, FUN = "sum")
  # Beverton-Holt stock recruitment relationship
  # R0 = = Rt*SPR/(SPR-(1-h)/(4h))
  # Rb = SBR0 R0 (1-h) / (4h)
  
  SBR0 <- vapply(seq_len(nsim), function(i) {
    with(dr_df, fishblicc::RSPR_0(Galpha[i], Gbeta[i], Mk[i], blicc_ld))
  }, numeric(1))
  SPR <- unlist(dr_df$SPR)
  
  # R0 = Rt*SPR/(SPR-(1-h)/(4h))
  steepness_increased <- FALSE
  denom <- SPR - (1 - h_steepness) / (4 * h_steepness)
  if (any(denom < 0)) {
    h_stp <- vapply(seq_along(denom), 
                    FUN = \(i) if (denom[i]<=0) return(0.01 + 1/(4*SPR[i]+1)) else return(h_steepness), 
                    FUN.VALUE = 0)
    h_steepness <- max(h_stp)
    steepness_increased <- TRUE
  }
  R0 <- recruit * SPR / (SPR - (1 - h_steepness) / (4 * h_steepness))
  
  lRd0 <- array(0, dim = c(nsim, PYN))
  # Find R0 consistent with current recruitment (i.e. catch)
  
  TargetRecruitment <- function(Rmul) {
    Rb <- Rmul * R0 * SBR0 * (1 - h_steepness) / (4 * h_steepness)
    Recruitment <- BH_recruitment(lR0 = log(Rmul * R0),
                                  Rb = Rb,
                                  lRs = recruit_cv)
    
    mean(Recruitment(xSSB, lRd0) - recruit)
    
    return(mean(Recruitment(xSSB, lRd0) - recruit))
  }
  
  s1 <- s2 <- -1
  repeat {  
    Rmul_interval <- c(0.001, 100)
    s1 <- sign(TargetRecruitment(Rmul_interval[1]))
    for (i in 1:10) {
      s2 <- sign(TargetRecruitment(Rmul_interval[2]))
      if (s1!=s2) {
        break
      } else {
        Rmul_interval[1] <- Rmul_interval[2]
        Rmul_interval[2] <- Rmul_interval[2]*2
      }
    }
    if (s1!=s2) break
    steepness_increased <- TRUE
    h_steepness <- round(h_steepness + 0.01, 2)
    if (h_steepness==1.00) break  
  }
  
  rm(s1, s2)
  if (h_steepness==1.00) {
    Recruitment <- fixed_recruitment(lR0 = log(R0), lRs = recruit_cv)
    warning(paste("No valid SR found. Recruitment fixed at ", format(R0, digits=3)))
  } else {
    res <- uniroot(TargetRecruitment, interval = Rmul_interval)
  
    if (steepness_increased)
      warning(paste("Steepness adjusted to achieve stable status: ", format(h_steepness, digits=3)))
  
    R0 <- res$root * R0
    Rb <- R0 * SBR0 * (1 - h_steepness) / (4 * h_steepness)
    R00 <- R0 * (5 * h_steepness - 1) / (4 * h_steepness) # longterm recruitment with no fishing
    Recruitment <- BH_recruitment(lR0 = log(R0),
                                Rb = Rb,
                                lRs = recruit_cv)
  
    rm(res)
  }

  cat(
    paste(
      "BH recruitment applied with (means) steepness",
      format(h_steepness),
      ", CV =",
      format(mean(recruit_cv) * 100),
      "%,",
      format(mean(autocorr)),
      "autocorrelation and log R0 =",
      format(log(mean(R0))),
      "\n"
    )
  )
  
  index_func <- switch(
    index_type,
    MeanLength = meanLength_func(blicc_ld),
    SPR = {
      SPR_func_para(blicc_ld)
    },
    SPR_F = local({
      blicc_ld1 <- blicc_ld
      blicc_ld1$poLinfm <- mean(Linf)
      blicc_ld1$polGam <- log(mean(Galpha))
      blicc_ld1$polMkm <- log(mean(dr_df$Mk))
      blicc_ld1$polSm <- log(rowMeans(simplify2array(dr_df$Sm)))
      blicc_ld1$polNB_phim <- log(mean(NBphi))
      SPR_onlyF_para(blicc_ld1)
    }),
    stop("index_type not recognised.\n")
  )

  # Determine SSB and status
  SSB[[1]] <- xSSB
  
  if (SRDelay > 1) {
    for (i in 2:SRDelay) {
      SSB[[i]] <- xSSB
    }
  }
  recruit <- Recruitment(xSSB, lRd0[, c(1, 2L)]) / (sum(incrementAge)*Ngtg)

  # Calculate first index value
  freq <- simplify2array(blicc_ld$fq)
  dim(freq) <- c(LN, 1L, NG)
  pjIndex[[1]] <- rep(index_func(freq), nsim)
  rm(freq)
  
  # Calculate indexes for computing catch and length comps
  len_sim_idx <- list()
  sim_vals <- seq_len(nsim)
  for (mi in seq_len(yr_steps)) {
    len_vals <- sort(unique(len_idx[[mi]]))
    len_sim_idx[[mi]] <- as.matrix(expand.grid(len_vals, sim_vals))
  }
  Ngtg <- as.double(Ngtg)
  
  cat("Calculating MSY reference levels.\n")
  if (index_type %in% c("SPR", "SPR_F")) {
    # Parallel processing
    old_plan <- future::plan()          # save current plan
    on.exit(future::plan(old_plan), add = TRUE)  # restore plan safely
    cores <- future::availableCores()
    future::plan(multisession, workers = cores)  # works on Windows/macOS/Linux
  }
  MSY_rp <- fb_calc_MSY_refpt(
    Ngtg,
    Nage,
    LN,
    yr_steps,
    mF[[1]],
    selectivity,
    mM,
    sim_idx,
    len_idx,
    wtmat,
    wt,
    R0,
    Rb,
    R00,
    SampleSize,
    index_func,
    incrementAge
  )
  #MSY_rp$ld <- blicc_ld
  #return(MSY_rp)
  cat("Calculating reference points.\n")
  ref_pt <- c(standard_risk_ref_pt(), MSY_rp)
  ref_pt$rp_type <- rp_type
  ref_pt$TRP <- switch(rp_type,
                       Current = 1.0,
                       MSY = 1.0,
                       SSB0 = 0.4)
  ref_F <- switch(rp_type,
                  Current = mF[[1]],
                  MSY = MSY_rp$FMSY,
                  SSB0 = mF[[1]])
  SSBtar <- switch(
    ref_pt$rp_type,
    Current = SSB[[1]] * ref_pt$TRP,
    MSY = MSY_rp$SSBMSY * ref_pt$TRP,
    #* ref_pt$TRP
    SSB0 = MSY_rp$SSB0 * ref_pt$TRP
  )
  ref_pt$B_tar <- SSBtar
  tar_rp <- fb_find_refpts(
    Ngtg,
    Nage,
    LN,
    yr_steps,
    ref_F,
    selectivity,
    mM,
    sim_idx,
    len_idx,
    wtmat,
    wt,
    R0,
    Rb,
    SampleSize,
    index_func,
    incrementAge,
    recalc_F = (ref_pt$rp_type == "SSB0"),
    SSBtar
  )
  ref_pt <- c(ref_pt, tar_rp)
  ref_pt$B_lim <- ref_pt$B_tar*ref_pt$LRP
  
  #remove everything we don't need
  rm(ref_F, tar_rp, MSY_rp, rp_type, sim_vals, len_vals, len, proj_length,
    catch, vbK,
    dr_df,
    ZaL,
    Ca,
    mort,
    LF,
    SSBtar,
    xSSB,
    p,
    pop2,
    i,
    mi,
    gi,
    rowi,
    SPR,
    SBR0,
    gtg_idx,
    R00,
    lRd0, steepness_increased
  )
  cat("MSE function completed.\n")

  # ><> #  # ><> #  # ><> #  # ><> #
  ###  GENERATED FUNCTION    # ><> #
  # ><> #  # ><> #  # ><> #  # ><> #
    
  gf <-  function(trIndex, trControl, control_type = "Effort", change_limit, ma, ctrl_pF = 1) {
    #Define the HCR: Create HCR functions
    .Random.seed <- sim_stream
    if (sim_recruit == "shared")
      lRd <- matrix(rep(rnorm(PYN + 1), each = nsim),
                    nrow = nsim,
                    ncol = PYN + 1)
    else
      lRd <- matrix(rnorm(nsim * (PYN + 1)), nrow = nsim, ncol = (PYN +
                                                                    1))
    Rec[[1]] <- Recruitment(SSB[[1]], cbind(rep(0, nsim), lRd[, 1])) # Add recruitment variation
    recruits <- rep(Rec[[1]] / (sum(incrementAge)*Ngtg), each = Ngtg)
    
    HCRIndex <- fb_create_HCR_index_func(ma, index_func)
    
    CalcControl <- fb_create_control_func(
      refControl = ref_control,
      trIndex = trIndex,
      trControl = trControl,
      change_limit = change_limit
    )
    
    if (index_type %in% c("SPR", "SPR_F")) {
      # Parallel processing
      old_plan <- future::plan()          # save current plan
      on.exit(future::plan(old_plan), add = TRUE)  # restore plan safely
      cores <- future::availableCores() - 4L
      future::plan(multisession, workers = cores)  # works on Windows/macOS/Linux
    }
    yi <- 1L
    mi <- 1L
    CW <- LF <- list()
    CaL <- array(0, dim = c(LN, nsim, NG))
    CW[[1]] <- double(nsim)
    for (ti in 1:PN) {
      ZaL <- mM
      for (gi in seq_len(NG)) {
        FaL[[gi]] <- selectivity[[gi]] * rep(mF[[yi]][, gi] * dt, each = LN)
        ZaL <- ZaL + FaL[[gi]]
      }
      surv <- exp(-ZaL)
      # Survival
      pop2 <- pop1 * surv[surv_idx[[mi]]]
      # Deaths
      mort <- pop1 - pop2
      for (gi in seq_len(NG)) {
        p <- FaL[[gi]] / ZaL
        Ca <- mort * p[surv_idx[[mi]]] # catch numbers
        CW[[yi]] <- CW[[yi]] + tapply(Ca * wt[[mi]], INDEX = sim_idx, FUN =
                                        "sum")
        # sum up in categories
        cl <- tapply(Ca,
                     INDEX = list(len_idx[[mi]], sim_idx),
                     FUN = "sum")
        cl[is.na(cl)] <- 0
        CaL[cbind(len_sim_idx[[mi]], gi)] <-
          CaL[cbind(len_sim_idx[[mi]], gi)] + cl
      }
      
      if (mi == AssessmentPeriod) {
        for (gi in seq_len(NG)) {
          cs <- colSums(CaL[, , gi])
          cs[cs == 0] <- 1   # avoid divide-by-zero
          CaL[, , gi] <- CaL[, , gi] / rep(cs, each = LN) * SampleSize
        }
        LF[[yi]] <- array(rnbinom(
          length(CaL),
          mu   = as.vector(CaL),
          size = rep(rep(NBphi, each = LN), NG)
        ), dim = dim(CaL))
        pjIndex[[yi + 1L]] <- HCRIndex(pjIndex[[yi]], LF[[yi]])
        mF[[yi + 1L]] <- CalcControl(pjIndex[[yi + 1L]], mF[[yi]])
        CaL[] <- 0
      }
      
      if (incrementAge[mi]) {
        pop1[Age_to_idx] <- pop2[Age_from_idx]
        pop1[rec_idx] <- recruits # add recruits as split gtg age 1 for each sim
        pop1[plusgp_idx] <- pop2[plusgp_idx] + pop1[plusgp_idx]
      } else {
        pop1 <- pop2
      }
      
      # Recruitment is calculated only once in a year 
      if (mi == yr_steps) {
        Rec[[yi]] <- Recruitment(SSB[[yi]], lRd[, c(yi, yi + 1L)]) # nsim recruitments
        recruits <- rep(Rec[[yi]] / (sum(incrementAge)*Ngtg), each = Ngtg)
        SSB[[yi+SRDelay]] <- tapply(pop1 * wtmat, INDEX = sim_idx, 
                                    FUN = "sum")
        mi <- 1L
        yi <- yi + 1L
        CW[[yi]] <- double(nsim)
      } else {
        mi <- mi + 1L
      }
    } #ti
    CW[[yi]] <- NULL
    mFarray <- aperm(simplify2array(mF), c(1, 3, 2)) # gear last term
    return(
      list(
        stock_assessment = stock_assessment,
        SSB = simplify2array(SSB),
        Rec = simplify2array(Rec),
        mF = mFarray,
        LF = LF,
        CW = simplify2array(CW),
        pjIndex = simplify2array(pjIndex),
        blicc_ld = blicc_ld,
        ref_pt = ref_pt,
        HCR = list(
          nsim = nsim,
          TN = 1L,
          PN = PN,
          PYN = PYN,
          StartYear = StartYear,
          trIndex = trIndex,
          trControl = trControl,
          control_type = control_type,
          index_type = index_type,
          change_limit = change_limit,
          ma = ma,
          SampleSize = SampleSize
        )
      )
    )
  } #gf
    
  return(gf)
}



