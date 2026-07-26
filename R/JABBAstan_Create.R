 ##  ><> ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>
 ##  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>
 ##  ><> ##  ><>
 ##  ><> ##    ><>     JABBAStan HCR Functions v0.96
 ##  ><> ##  ><>
 ##  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>
 ##  ><> ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>


#' Create a HCR MSE function using a Stan fit of the JABBA
#'
#' The HCR function returns simulation results from HCR parameters using a
#' jabba model fit results.
#'
#' @param stan_fit A stanfit object
#' @param stan_dat Data list used to fit the Stan model
#' @param control_type Type of control applied by HCR either "Effort" or "Catch"
#' @param proj_length   Projection length in years
#' @param nsim Number of simulations to run
#' @param ref_pt Specifies risk based reference points in a list. Otherwise
#'   defaults apply.
#' @return A function that will apply the HCR using the stan model fit and
#'   specified HCR parameters: index-control inflection point vectors, the
#'   control change limit and moving average parameter.
#' @export
#'
create_ptStan_MSE <- function(stan_fit,
                              stan_dat,
                              proj_length = 50,
                              nsim = 1000,
                              ref_pt = standard_risk_ref_pt()) {
  # candidate HCR definitions: not used in this function, but recorded for later evaluation
  # can be obtained from the function using "get"
  # Reference points are relative to B0 for Schaefer
  stock_assessment <- "JABBAstan" 
  minstatus <- 0.005
  NGears <- 1
  if (NGears==0L) NGears <- 1L
  # Extract the data components that will be used from the fit
  Par <- rstan::extract(stan_fit, permuted = TRUE, inc_warmup = FALSE, include = TRUE) |>
    tibble::as_tibble()
  
  # Dimensions
  TN <- as.integer(stan_dat$TN)
  PN <- as.integer(proj_length)
  PTN <- PN + TN
  start_year <- stan_dat$YR
  
  Avg_CPUE <- NULL #Only single gear   #create_mean_CPUE(jabba_fit)

  #if (is.null(rho)) rho <- -log(1-median(Bio$H))/1  #H=(1-e(-rho*f))
  #obserr <- sqrt(pull(jabba_fit$pars,
  #                    Median)[startsWith(rownames(jabba_fit$pars), "tau2")])

  # Extract sufficient parameters from the fit for the n simulations
  if (nsim > nrow(Par)) {
    nsim <- nrow(Par)
  } else if (nsim < nrow(Par)) {
    # Select a random sample of iterations for the simulations
    ii <- 1:nrow(Par)
    ii <- ii[Par$r < 2] # preferentially exclude unstable r
    sn <- sample(ii, size=nsim, replace=FALSE)

    sn <- sort(sn)
    Par <- Par |>
      dplyr::slice(sn)

    rm(sn, ii)
  }
  obserr <- dplyr::pull(Par, ce_cv)

  #Calculate past biomass etc.
  # Bio <- jabba_fit$kbtrj |>
  #   dplyr::select(year, iter, harvest:BB0)
  # Bio <- dplyr::mutate(Bio, tim = as.integer(year - ymin))

  Effort <- with(stan_dat, TCA_ca / ( sum(TCE_ca) / sum(TCE_ef) )) # Fixed effort
  pB <- matrix(0, ncol=PTN+1, nrow=nsim)
  Ft <- C <- matrix(0.0, nrow=nsim, ncol=PTN)

  for (i in seq_len(nsim)) {
    pvFt <- with(Par, Effort * exp(lq[i]) + dF[i,])
    pS <- exp( -pvFt )
    pB[i, 1:TN] <- with(Par, pt_BioDyn(sNut[i,], Nus[i], pS, P0[i], r[i], m[i]))
    Ft[i, 1:TN] <- pvFt
    C[i, 1:TN] <- (1-pS) * pB[i, 1:TN] * Par$Binf[i]
  }
  rm(pS, pvFt)

  # Single gear
  pvCPUE <- array(0, dim=c(TN, 1))
  pvCPUE[stan_dat$TCE_t,] <- with(stan_dat, TCE_ca / TCE_ef)

  # Model parameters
  r <- dplyr::pull(Par, r)
  lsigma <- dplyr::pull(Par, Nus)
  # The following will need adaptation if more than one gear
  q <- exp(dplyr::pull(Par, lq))
  m_1 <- pull(Par, m) - 1
  r_m_1 <- r/m_1
  Binf <- dplyr::pull(Par, Binf)
  prod_fun <- function(pBti) {
    return(pBti * (1 + r_m_1 * (1 - pBti^m_1)) * rlnorm(nsim, 0, lsigma)) }

  B_trial <- apply(sweep(pB[, 1:TN], MARGIN=1, STATS=as.array(Binf), FUN = "*"), 1, min)
  dat <- stan_dat
  # Complete the default reference points
  ref_pt$B_tar <- dplyr::pull(Par, BMSY)*ref_pt$TRP
  ref_pt$B_lim <- ref_pt$B_tar*ref_pt$LRP
  ref_pt$F_tar <- with(Par, - log((m-1) / (m - 1 + r * (1-ref_pt$B_tar^(m-1)))))
  ref_pt$rp_type <- "MSY"
  Par <- dplyr::select(Par, -BMSY)
  
  # Tidy up
  rm(stan_fit, proj_length, stan_dat, i)

  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #
  ###  GENERATED FUNCTION    # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #
  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #
  
  function(trIndex, trControl, control_type, change_limit, ma, ctrl_pF = 1) {
    #Initial Control/Index
    NCtrl <- length(control_type)
    if (length(ctrl_pF)==1L) ctrl_pF <- rep(ctrl_pF, NCtrl)

    if (NCtrl > 1) {
      if (any(NCtrl != c(length(trIndex), length(trControl), length(ctrl_pF)) |
          !c(is.list(trIndex) | !is.list(trControl)))) {
        stop("HCR: Number of controls must be the same for all arrays.")
      }}
    # if (!is.list(trIndex)) trIndex <- list(trIndex)
    # if (!is.list(trControl)) trIndex <- list(trControl)

    UpdateIndex <- pt_create_HCR_index(ma, func=Avg_CPUE, obserr)
    CalcControl <- pt_create_linear_control(trIndex, trControl,
                                            change_limit=change_limit)

    if (NCtrl == 1L)
      C_trial <- max(trControl)
    else
      C_trial <- max(trControl[[1L]])

    if (control_type[1L] == "Catch") {
      F_trial <- -log(1 - C_trial/ B_trial)
    } else {
      F_trial <- q * C_trial
    }

    pvIndex <- double(TN)
    pvControl <- matrix(0, ncol=TN, nrow=NCtrl)
    pvIndex[1] <- mean(pvCPUE, na.rm=T)   # First CPUE may be NA, so use the overall average to start
    pvControl[ , 1] <- ifelse(is.list(trControl), trControl[[1]][2], trControl[2]) # Need a better start than this...

    for (yi in 2:TN){
      pvIndex[yi] <- UpdateIndex(pvIndex[yi-1L], pvCPUE[yi,])
      pvControl[ , yi] <- drop(CalcControl(pvIndex[yi], pvControl[ , yi-1L]))
    }

    ###  PROJECTION   # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #
    pjIndex <- matrix(0, nrow=nsim, ncol=(PN+2L))
    pjControl <- array(0, dim=c(nsim, NCtrl, (PN+2L)))
    pjIndex[, 1L] <- pvIndex[TN]
    pjControl[, , 1L] <- as.array(rep(as.vector(pvControl[ , TN]), each=nsim), dim=c(nsim, NCtrl))  # sim, ctrl, year

    for (pi in seq(PN+1L)) {
      ti <- TN + pi - 1L

      pB1 <- prod_fun(pB[, ti])
      pB1[pB1 < minstatus] <- minstatus

      # Apply Control
      Ft1 <- F_trial
      for (cj in seq(NCtrl)) {
        if (control_type[cj]=="Catch") {
          exp_F <- pmax(minstatus, 1 - pjControl[ , cj, pi] / (Binf*pB1))
          F_limit <- Ft1*(1 - ctrl_pF[cj]) - log(exp_F)
        } else if (control_type[cj]=="Effort") {
          F_limit <- Ft1*(1 - ctrl_pF[cj]) + q * pjControl[ , cj, pi]
        } else { # Opportunities
          F_limit <- Ft1 * (1 - ctrl_pF[cj] + ctrl_pF[cj] * pjControl[ , cj, pi])
        }
        Ft1 <- pmin(Ft1, F_limit)
      }

      Ft[ , ti] <- Ft1
      pB[, ti+1L] <- pB1 * exp(- Ft[, ti])
      C[ , ti] <- (pB1 - pB[, ti+1L]) * Binf
      CPUE <- q * C[ , ti] / Ft[ , ti]

      pjIndex[, pi+1L] <- UpdateIndex(pjIndex[, pi], CPUE)
      pjControl[,, pi+1L] <- CalcControl(pjIndex[, pi+1L], pjControl[, , pi])
    } #ti

    return(list(stock_assessment = stock_assessment,
                pB=pB, C=C, Ft=Ft, Par=Par,
                pvIndex=pvIndex, pvControl=pvControl,
                pjIndex=pjIndex, pjControl=pjControl,
                HCR=list(nsim=nsim, TN=TN, PN=PN, PTN=PTN,
                         start_year = start_year,
                         trIndex=trIndex, trControl=trControl,
                         control_type=control_type,
                         change_limit=change_limit, ma=ma),
                dat = dat,
                ref_pt = ref_pt))
  }
}




#' Create the empirical abundance index update function.
#'
#' The index function created calculates a moving average from one or more
#' abundance indices, usually CPUE.
#'
#' @param  ma  Moving average parameter for index time series
#' @param  func Function used to combine multiple CPUE series into a single index
#' @param  se  Observation standard error to be added in to the index
#' @return A function taking a single CPUE list as a parameter that can be used
#'   to calculate the HCR index
#' @export
#'
pt_create_HCR_index <- function(ma, func=NULL, se) {
  default_factor <- 1.05  # default 5% increase to allow recovery
  if (is.null(func)) # Single gears
    return(
      function(indx_1, CPUE) {
        indx_0 <- CPUE*rlnorm(length(CPUE), 0, se)
        invalid <- is.na(indx_0) | (indx_0 <= 0)
        indx_0[invalid] <- indx_1[invalid]*default_factor
        return(ma*indx_0 + (1-ma)*indx_1)
      })
  else
    return(    # Allow for multiple gears
      function(indx_1, CPUE) {
        # CPUE can be a matrix
        if (is.matrix(CPUE)) {
          err <- matrix(rlnorm(length(CPUE), 0, se), ncol=length(se), byrow=T)
          indx_0 <- apply(CPUE*err, MARGIN=1, func)
        } else {
          err <- rlnorm(length(CPUE), 0, se)
          indx_0 <- func(array(CPUE*err, dim=c(1, length(CPUE))))  # Apply func vector
        }
        invalid <- is.na(indx_0) | (indx_0 <= 0)
        indx_0[invalid] <- indx_1[invalid]*default_factor
         return((ma*indx_0 + (1-ma)*indx_1))
        }
      )
}


#' Calculate a linear sliding control to be applied based on the index value
#'
#' The function creates a function to calculate a control based on the index,
#' previous control and HCR parameters held in the HCR list.
#'
#' @param trIndex A sorted matrix of index inflexion points from low to high.
#'   Same size as `trControl`
#' @param trControl A sorted matrix of control inflection points defining linear
#'   HCR.
#' @param change_limit An annual limit on the control change as a proportion
#' @param control_type Either "Effort" or "Catch" or "Both"
#' @return A function to calculate the new control from the index (vector) and previous
#'   year's controls (matrix)
#' @export
#'
pt_create_linear_control <- function(trIndex, trControl,
                                     change_limit=NA) {
  if (is.list(trIndex) & length(trIndex) > 1L) { #several controls
    NCtrl <- length(trIndex)
    trIndex <- lapply(trIndex, \(x) c(0, x, Inf))
    trControl <- lapply(trControl, \(x) c(x[1], x, x[length(x)]))
    if (all(is.na(change_limit)))
      return(
        function(indx, prev_con){
          # Calculate the control for each index value
          con <- array(0, dim=c(length(indx), NCtrl))
          for (cj in seq(NCtrl)) {
            ii <- findInterval(indx, trIndex[[cj]])
            con[, cj] <- trControl[[cj]][ii] +
              (trControl[[cj]][ii+1L] - trControl[[cj]][ii]) *
              (indx-trIndex[[cj]][ii])/(trIndex[[cj]][ii+1L]-trIndex[[cj]][ii])
          }
          return(con)
        }
      )
    else
      return(
        function(indx, prev_con){
          # Calculate the control for each index value
          con <- array(0, dim=c(length(indx), NCtrl))
          for (cj in seq(NCtrl)) {
            ii <- findInterval(indx, trIndex[[cj]])
            con[, cj] <- trControl[[cj]][ii] +
              (trControl[[cj]][ii+1L] - trControl[[cj]][ii]) *
              (indx-trIndex[[cj]][ii])/(trIndex[[cj]][ii+1L]-trIndex[[cj]][ii])
            if (!is.na(change_limit[cj])) {
              prop_change <- (con[, cj]-prev_con[, cj])/prev_con[, cj]
              var_limit <- (abs(prop_change) > change_limit[cj])
              if (any(var_limit)) {
                prop_change <- 1 + sign(prop_change[var_limit])*change_limit[cj]
                con[var_limit, cj ] <- prev_con[var_limit, cj] * prop_change
              }
            }
          }
          return(con)
        }
      )
  }
  # Single control
  trIndex <- c(0, unlist(trIndex), Inf)
  trControl <- c(trControl[1], unlist(trControl), trControl[length(trControl)])
  if (is.na(change_limit))
    return(
      function(indx, prev_con){
        # Calculate the control for each index value
        con <- double(length(indx))
        ii <- findInterval(indx, trIndex)
        con <- trControl[ii] +
          (trControl[ii+1L] - trControl[ii]) *
          (indx-trIndex[ii])/(trIndex[ii+1L]-trIndex[ii])
        return(as.array(con, dim=c(length(indx), 1)))
      }
    )
  else
    return(
      function(indx, prev_con){
        # Calculate the control for each index value
        con <- double(length(indx))
        ii <- findInterval(indx, trIndex)
        con <- trControl[ii] +
          (trControl[ii+1L] - trControl[ii]) *
          (indx - trIndex[ii])/(trIndex[ii+1L] - trIndex[ii])

        prop_change <- (con-prev_con)/prev_con
        var_limit <- (abs(prop_change) > change_limit)
        if (any(var_limit)) {
          prop_change <- 1 + sign(prop_change[var_limit])*change_limit
          con[var_limit] <- prev_con[var_limit] * prop_change
        }
        return(as.array(con, dim=c(length(indx), 1)))
      }
    )
}


#' Collapses multiple CPUE indices into single weighted geometric mean index
#'
#' A function is returned taking a matrix of CPUE with columns for each
#' CPUE index and convert to a single vector index. The indices are assumed to
#' be log-normally distributed, so the geometric mean is taken weighted by the
#' median JABBA observation error estimates (tau2).
#'
#' @inheritParams create_ptStan_MSE
#' @return A function taking a matrix with CPUE in columns and returning a
#'   single vector index
#' @export
#'
pt_create_mean_CPUE <- function(stanfit) {
  q <- unlist(rstan::extract(stanfit, pars="lq"))
  if (length(q)<=1) {
    return(
      NULL
      # function(vCPUE) {
      #   return(vCPUE)
      # }
  )
  } else {
    ce_cv <- unlist(rstan::extract(stanfit, pars="ce_cv"))
    # CPUE = qB e^E  E~N(0, tau2)
    weights <- 1/ce_cv^2
    weights <- weights / sum(weights)
    rm(q, tau2)
    return(
      function(CPUE) {
        is_na <- is.na(CPUE)
        if (any(is_na)) {
          not_na <- !is_na
          if (any(not_na))
            gm <- exp(sum(log(CPUE[not_na])*weights[not_na]/sum(weights[not_na])))
          else
            gm <- NA
        } else
          gm <- exp(sum(log(CPUE)*weights))
        return(gm)
      }
    )
  }
}


#' Calculates the projected biomass from input parameters
#'
#' @param nut Population standard log-normal deviates for the process error 
#'   time series
#' @param nus Population process error scale parameter (log-normal)
#' @param pS  Proportion survival of biomass for each time 
#' @param P0 Initial stock as as a proportion of Binf
#' @param r Intrinsic rate of increase
#' @param m Pella-thomlinson m parameter
#' @return A vector of biomass projected based on the biomass dynamics 
#'   parameters
#' @noRd
#'
pt_BioDyn <- function(nut, nus, pS, P0, r, m) {
  # Production model
  TN <- length(pS)
  eBt <- double(TN)
  Bdev <- exp(nut * nus)
  
  Bt <- P0 * Bdev[1]
  eBt[1] <- Bt
  
  m1 <- m - 1.0
  for (ti in 2:TN) {
    eBt[ti] <- (Bt + (r / m1) * Bt * (1.0 - Bt^m1)) * pS[ti - 1] * Bdev[ti]
    Bt <- eBt[ti]
  }
  return(eBt)
}


#' Extracts MSY reference points from a JABBAstan model fit
#'
#' This only works on a single gear. To extend it, a function will need to be
#' provided to convert catch-effort to the HCR index.
#'
#' @inheritParams create_ptStan_MSE
#' @return List of MSY reference points for different variation
#' @export
#'
pt_MSY_mean_refpt <- function(stan_fit) {
  # Effort <- with(stan_dat, TCA_ca / ( sum(TCE_ca) / sum(TCE_ef) )) # Fixed effort
  PP <- rstan::extract(stan_fit, permuted = TRUE, inc_warmup = FALSE, include = TRUE)
  msyrp <- pt_calc_MSY_refpt(PP)
  return(lapply(msyrp, mean))
}


#' Extracts MSY reference points from a JABBAstan model fit
#'
#' This only works on a single gear. To extend it, a function will need to be
#' provided to convert catch-effort to the HCR index.
#'
#' @param PP Parameter list from a JABBAstan fit
#' @return List of MSY reference points for different variation
#' @noRd
#'
pt_calc_MSY_refpt <- function(PP) {
  FMSY <- with(PP, log(1 + r * (1-BMSY^(m-1)) / (m-1)))
  CPUEMSY <- exp(PP$lq) * PP$Binf * (1-exp(-FMSY)) * PP$BMSY / FMSY
  return(  list(BMSY = mean(PP$BMSY),
                IMSY = mean(CPUEMSY),
                FMSY = mean(FMSY),
                fMSY = mean(FMSY * exp(-PP$lq)),
                MSY = mean(PP$MSY))  
  )
}



