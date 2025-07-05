
# #### MSE function creation ####

#' Create a HCR MSE function 
#' 
#' The HCR function returns simulation results from HCR parameters using a 
#' JABBA model fit results.
#' 
#' The the same state space lognormal random number is applied across 
#' simulations, so HCR performance is more comparible.
#' 
#' @param jabba_fit A list returned from package JABBA `fit_jabba` function
#' @param control_type Type of control applied by HCR either "Effort" or "Catch"
#' @param ProjLength   Projection length in years
#' @param nsim Number of simulations to run
#' @param ref_pt Specifies risk based reference points in a list. Otherwise 
#'   defaults apply.
#' @param rho Conversion factor from effort (f) to fishing mortality, so 
#'   H=(1-e(-rho*f)). Only used if control_type is "Effort". Control defaults to 
#'   fishing mortality (rho=1).
#' @return A function that will apply the HCR using the JABBA model fit and 
#'   specified HCR parameters: index-control inflection point vectors, the 
#'   control change limit and moving average parameter.
#' @export
#' 
create_HCR_MSE <- function(jabba_fit, 
                           control_type = "Effort",
                           ProjLength=50, 
                           nsim = 1000,
                           ref_pt = standard_risk_ref_pt(),
                           rho = 1) {
  # candidate HCR definitions: not used in this function, but recorded for later evaluation
  # can be obtained from the function using "get"
  # Reference points are relative to B0 for Schaefer
  minstatus <- 0.005 
  NGears <- sum(substr(rownames(jabba_fit$pars), 1, 2)=="q.")
  if (NGears==0L) NGears <- 1L
  # Extract the data components that will be used from the fit
  Par <- jabba_fit$pars_posterior 
  Bio <- jabba_fit$kbtrj |>
    dplyr::select(year, iter, harvest:BB0)
  
  ymin <- min(Bio$year) - 1L
  # Dimensions
  TN <- max(Bio$year) - ymin
  PN <- as.integer(ProjLength)
  PTN <- PN + TN
  
  Avg_CPUE <- create_mean_CPUE(jabba_fit)
  obserr <- sqrt(dplyr::pull(jabba_fit$pars,
                      Median)[startsWith(rownames(jabba_fit$pars), "tau2")])
  
  # Extract sufficient parameters from JABBA fit for the n simulations
  if (nsim > max(Bio$iter)) {
    nsim <- as.integer(max(Bio$iter))
  } else if (nsim < max(Bio$iter)) {
    # Select a random sample of iterations for the simulations
    ii <- 1:max(Bio$iter)
    ii <- ii[Par$r < 2]
    sn <- sample(ii, size=nsim, replace=FALSE)
    
    sn <- sort(sn)
    Par <- Par |>
      dplyr::slice(sn)
    
    Bio <- Bio |>
      dplyr::rename(siter = iter) |>
      dplyr::filter(siter %in% sn) |>
      dplyr::mutate(iter = match(siter, sn))
    rm(sn, ii)
  }
  Bio <- dplyr::mutate(Bio, tim = as.integer(year - ymin))
  
  pB <- matrix(0, ncol=PTN+1, nrow=nsim)
  pB[cbind(Bio$iter, Bio$tim)] <- Bio$BB0   # extract B/B0 estimates as matrix
  #NGear <- ncol(jabba_fit$inputseries$cpue) - 1L
  #Catch <- jabba_fit$inputseries$catch$catch
  pvCPUE <- as.matrix(jabba_fit$inputseries$cpue[,-1])
  
  # Model parameters   
  r <- dplyr::pull(Par, r)
  lsigma <- sqrt(dplyr::pull(Par, sigma2))
  if (NGears > 1L) 
    q <- as.matrix(dplyr::select(Par, tidyselect::starts_with("q."))) #Excludes auxiliary
  else
    q <- as.matrix(dplyr::select(Par, q)) #Excludes auxiliary
  
  if (jabba_fit$settings$model.type=="Schaefer") {
    invBMSY <- 2
    invHMSY <- 2 / r     
    prod_fun <- function(pBti) {
      pBti * (1 + r * (1 - pBti)) * rlnorm(1, 0, lsigma)}
  } else if (jabba_fit$settings$model.type %in% c("Fox", "Pella_m")) {
    m_ <- dplyr::pull(Par, m)
    m_1 <- m_ - 1
    r_m_1 <- r/(m_ - 1)
    invBMSY  <- m_^(1 / m_1) # BMSY/K = m^(-1/(m-1))
    invHMSY <- m_ / r        # HMSY = r / m
    prod_fun <- function(pBti) {
      pBti * (1 + r_m_1 * (1 - pBti^m_1)) * rlnorm(1, 0, lsigma)}
  } else {
    stop("Error: Model type not recognised.")
  }
  
  # Tidy up
  rm(jabba_fit, ProjLength, ymin, m_)
  
  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #
  ###  GENERATED FUNCTION    # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> # 
  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #
  
  function(trIndex, trControl, change_limit, ma) {
    #Initial Control/Index
    
    UpdateIndex <- create_empirical_index(ma, func=Avg_CPUE, obserr)
    CalcControl <- create_linear_control(trIndex, trControl,  
                                         change_limit=change_limit, 
                                         control_type=control_type)
    
    pvIndex <- double(TN)
    pvControl <- double(TN)
    pvIndex[1] <- mean(pvCPUE, na.rm=T)   # First CPUE may be NA, so use the overall average
    pvControl[1] <- trControl[2]
    
    for (yi in 2:TN){
      pvIndex[yi] <- UpdateIndex(pvIndex[yi-1L], pvCPUE[yi,])
      pvControl[yi] <- CalcControl(pvIndex[yi], pvControl[yi-1])
    }
    
    ###  PROJECTION   # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #  # ><> #
    
    pjIndex <- pjControl <- matrix(0, nrow=nsim, ncol=(PN+1))
    pjIndex[,1] <- pvIndex[TN]
    pjControl[,1] <- pvControl[TN]
    C <- matrix(0.0, nrow=nsim, ncol=PN)
    H <- C
    
    for (pi in 1:(PN)) {
      ti <- TN+pi-1
      #pB1 <- (pB[, ti] * (1 + r * (1 - pB[, ti])) * rlnorm(nsim, 0, lsigma)) - Schaefer
      #pB1 <- (pB[, ti] * (1 + r_m_1 * (1 - pB[, ti]^m_1)) * rlnorm(nsim, 0, lsigma) - Pella_m
      
      pB1 <- prod_fun(pB[, ti])
      pB1[pB1 < minstatus] <- minstatus
      
      if (control_type=="Catch") { 
        C[ , pi] <- pjControl[, pi] / Par$K
        H[ , pi] <- C[ , pi] / pB1
      } else {
        H[ , pi] <- (1 - exp(-pjControl[, pi] * rho))
        C[ , pi] <- H[ , pi] * pB1
      }
      
      invalid <- (pB1 <= C[, pi])
      C[invalid, pi] <- 0
      pB[, ti+1L] <- pB1 - C[, pi]
      CPUE <- apply(q, 2, function(x) x*pB1*Par$K)
      pjIndex[, pi+1L] <- UpdateIndex(pjIndex[, pi], CPUE)
      pjControl[, pi+1L] <- CalcControl(pjIndex[, pi+1L], pjControl[, pi])
    } #ti
    
    return(list(pB = pB, C = C, H = H, Par = Par, Bio = Bio, 
                pvIndex = pvIndex, pvControl = pvControl, 
                pjIndex = pjIndex, pjControl = pjControl, 
                invBMSY = invBMSY, invHMSY = invHMSY,
                HCR=list(nsim=nsim, TN=TN, PN=PN, PTN=PTN, 
                         trIndex=trIndex, trControl=trControl, 
                         control_type=control_type,
                         change_limit=change_limit, ma=ma),
                ref_pt = ref_pt))
  }
}


#' Generates combinations of index, controls, change limits and the index moving 
#' average parameter as list columns in a tibble, each row being a unique HCR.
#' 
#' A tibble is created defining full combinations of the defined parameter ranges.
#' 
#' @param  IMSY Index at MSY estimated from the JABBA model
#' @param  fMSY Control (TAC or effort) at MSY estimated from the JABBA model
#' @param  rel_index_range  (lower, upper) range for indices as proportions of 
#'   the MSY values
#' @param  rel_control_range  (lower, upper) range for controls as proportions 
#'   of the MSY values
#' @param  NBreaks The number of breaks including the ranges to generate the 
#'   test HCR
#' @param  NInflex The number of inflexion (trigger) points for the linear HCR
#' @param  change_limit Limit on the annual change of the control. NA implies 
#'   no limit
#' @param  ma  Moving average parameter for index time series
#' @return A tibble of all combinations of index, control, control change limit,
#'   and moving average parameters for the HCR. 
#' @export
#' 
define_HCR_test_range <- function(IMSY, fMSY, 
                                  rel_index_range = c(0.5, 1.0),
                                  rel_control_range = c(0.1, 1.0), 
                                  NInflex = 2L, 
                                  NBreaks = 5L, 
                                  change_limit = NA, 
                                  ma = 0.5) {
  NInflex <- as.integer(NInflex) 
  NBreaks <- as.integer(NBreaks) 
  Controls <- seq(rel_control_range[1], rel_control_range[2], length.out=NBreaks) * fMSY # Index values intervention points
  Indices <- seq(rel_index_range[1], rel_index_range[2], length.out=NBreaks) * IMSY  # Control levels to be applied
  
  trIndex <- list()
  trControl <- list()
  trCtrlFix <- list()
  Ir <- integer(NInflex)
  Ir[] <- 1L
  
  while (TRUE) {
    for (j in Ir[NInflex-1L]:NBreaks) {
      Ir[NInflex] <- j
      trIndex <- append(trIndex, list(Indices[Ir]))
      if (all(Ir==Ir[1]))
        trCtrlFix <- append(trCtrlFix, list(Controls[Ir])) #keep track of fixed controls
      else
        trControl <- append(trControl, list(Controls[Ir]))
    }
    lvl <- NInflex - 1L
    while (Ir[lvl]==NBreaks) {
      if (lvl>1) Ir[lvl] <- Ir[lvl-1L]
      lvl <- lvl - 1L
      if (lvl==0) break
    }
    if (lvl==0) break
    Ir[lvl] <- Ir[lvl] <- Ir[lvl] + 1L
  }
  
  fx_df <- tibble::tibble(trIndex = list(trIndex[[1]]), trControl = trCtrlFix, 
                          change_limit = change_limit[1], ma = ma[1])
  
  df <- tidyr::expand_grid(trIndex, trControl, change_limit, ma) |> 
    dplyr::bind_rows(fx_df) |>    # where controls are fixed, only 1 HCR is required
    dplyr::mutate(ID = dplyr::row_number()) |> 
    dplyr::select(ID, dplyr::everything())
  
  return(df)
}


#' Create the empirical abundance index update function. 
#' 
#' The index function created calculates a moving average from one or more
#' abundance indices, usually CPUE. 
#' 
#' The same lognormal error is applied across all simulations to make HCR more 
#' comparable.
#'
#' @param  ma  Moving average parameter for index time series
#' @param  func Function used to combine multiple CPUE series into a single index
#' @param  se  Observation standard error to be added in to the index
#' @return A function taking a single CPUE list as a parameter that can be used
#'   to calculate the HCR index
#' @export
#' 
create_empirical_index <- function(ma, func=NULL, se) {
  default_factor <- 1.05  # default 5% increase to allow recovery
  if (is.null(func)) # Single gears
    return(
      function(indx_1, CPUE) {
        #indx_0 <- CPUE*rlnorm(length(CPUE), 0, se)
        indx_0 <- CPUE*rlnorm(1, 0, se)
        invalid <- is.na(indx_0) | (indx_0 <= 0)
        indx_0[invalid] <- indx_1[invalid]*default_factor  
        return(ma*indx_0 + (1-ma)*indx_1)
      })
  else
    return(    # Allow for multiple gears
      function(indx_1, CPUE) {
        # CPUE can be a matrix
        if (is.matrix(CPUE)) {
          #err <- matrix(rlnorm(length(CPUE), 0, se), ncol=length(se), byrow=T)
          err <- matrix(rep(rlnorm(length(se), 0, se), length(CPUE)), ncol=length(se), byrow=T)
          indx_0 <- apply(CPUE*err, MARGIN=1, func)
        } else {
          #err <- rlnorm(length(CPUE), 0, se)
          err <- rlnorm(length(se), 0, se)
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
#' @param trIndex A sorted list of index inflexion points from low to high. 
#'   Same length as `trControl`
#' @param trControl A sorted list of control inflection points defining linear 
#'   HCR.
#' @param change_limit An annual limit on the control change as a proportion
#' @param control_type Either "Effort" or "Catch"
#' @return A function to calculate the new control from the index and previous
#'   year's control
#' @export
#' 
create_linear_control <- function(trIndex, trControl, 
                                  change_limit=NA, 
                                  control_type) {
  trIndex <- c(0, trIndex, Inf)
  trControl <- c(trControl[1], trControl, trControl[length(trControl)])
  if (is.na(change_limit))
    return(
      function(indx, prev_con){
        # Calculate the control for each index value
        con <- double(length(indx))
        ii <- findInterval(indx, trIndex)
        con <- trControl[ii] + 
          (trControl[ii+1L] - trControl[ii]) * 
          (indx-trIndex[ii])/(trIndex[ii+1L]-trIndex[ii])
        return(con)
      }
    )
  else
    return(
      function(indx, prev_con){
        # Calculate the control for each index value
        prop_change <- con <- double(length(indx))
        ii <- findInterval(indx, trIndex)
        con <- trControl[ii] + 
          (trControl[ii+1L] - trControl[ii]) * 
          (indx - trIndex[ii])/(trIndex[ii+1L] - trIndex[ii])
        tochange <- prev_con > 0
        
        prop_change[tochange] <- (con[tochange] - prev_con[tochange])/prev_con[tochange] 
        
        var_limit <- (abs(prop_change) > change_limit)
        if (any(var_limit)) {
          prop_change <- 1 + sign(prop_change[var_limit])*change_limit
          con[var_limit] <- prev_con[var_limit] * prop_change
        }
        return(con)
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
#' @inheritParams create_HCR_MSE
#' @return A function taking a matrix with CPUE in columns and returning a
#'   single vector index
#' @export
#'   
create_mean_CPUE <- function(jabbafit) {
  q <- jabbafit$pars$Median[substr(rownames(jabbafit$pars), 1, 2)=="q."]
  if (length(q)<=1) {
    return(    
      function(CPUE) return(CPUE)
    )
  } else {
    tau2 <- jabbafit$pars$Median[substr(rownames(jabbafit$pars), 1, 3)=="tau"]
    # CPUE = qB e^E  E~N(0, tau2)
    weights <- 1/tau2
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

#' Calculate the performance measures of a harvest control rule and return it as
#' a row in a data frame
#'
#' The performance measures consist of the 1) average catch, 2) the mean catch
#' range around the average (mean(abs(catch-mean_catch))), 3) the HCR Type 1
#' error rate, 4) the HCR Type 2 error rate, 5) the proportion of time that the
#' stock is below the limit reference point (usually 0.5 BMSY) and 6) the
#' proportion of time that the stock is outside the target range (default
#' 90-110% BMSY).
#'
#' The HCR error rates are scores indicating how well the HCR is responding. For
#' Type 1, When the stock is below the limit reference point and the index is
#' below the upper trigger but above the lower trigger, it scores one, whereas
#' if it is above both triggers it scores two. The mean of this value is taken.
#' Similarly for Type 2, responding by reducing harvest when the stock is at MSY
#' scores two when the reduction is minimised, one otherwise. The higher the
#' scores are, the worse the performance.
#'
#' @param HCR_sim Fishery HCR simulation result list from the MSE (see
#'   [create_HCR_MSE]).
#' @return Tibble row of performance measures
#' @export
#' 
HCR_performance <- function(HCR_sim) {
  B_tar_range <- HCR_sim$ref_pt$B_tar_range #c(0.90, 1.2)
  B_tar <- HCR_sim$ref_pt$B_tar
  B_lim <- HCR_sim$ref_pt$B_lim # 0.5
  MaxRisk <- HCR_sim$ref_pt$max_risk
  Mostly <- HCR_sim$ref_pt$mostly
  
  pBMSY <- with(HCR_sim, as.vector(apply(pB[ , (HCR$TN+1):(HCR$PTN)], 2, function(x) x*invBMSY)))
# ### MSY ####

  Index <- as.vector(HCR_sim$pjIndex[ , -1])
  
  Catch <- apply(HCR_sim$C, 2, function(x) x*HCR_sim$Par$K)
  Catch_Avg <- mean(Catch)
  Catch_Rng <- mean(abs(Catch-Catch_Avg))
  Catch_pcile <- quantile(Catch, probs=MaxRisk)
  lt_Blim <- sum(pBMSY < B_lim)/length(pBMSY)
  at_Btar <- sum((pBMSY >= B_tar_range[1]) & (pBMSY <= B_tar_range[2])) / length(pBMSY)
  gt_Btar <- sum(pBMSY > B_tar_range[2])/length(pBMSY)
  # Management response when not necessary
  Err_Type1 <- with(HCR_sim$HCR, 0.5*mean(((pBMSY > B_lim) & (Index < trIndex[1])) +
                          ((pBMSY > B_tar_range[1]) & (Index < trIndex[length(trIndex)]))))
  # No management response when it is necessary
  Err_Type2 <- with(HCR_sim$HCR, 0.5*mean(((pBMSY < B_lim) & (Index > trIndex[1])) +
                          ((pBMSY < B_lim) & (Index > trIndex[2]))))
  
  return(tibble::tibble(
    Catch_Avg = Catch_Avg, 
    Catch_Rng = Catch_Rng,
    Catch_pcile = Catch_pcile,
    lt_Blim = lt_Blim,
    at_Btar = at_Btar,
    gt_Btar = gt_Btar,
    State = at_Btar - lt_Blim,
    Err_Type1 = Err_Type1,
    Err_Type2 = Err_Type2
  ))
}


#' Run an HCR MSE based on a table of parameters and HCR MSE function 
#' 
#' The function runs the HCR MSE based on the set of HCR parameters in the 
#' dataframe and the provided MSE function. The performance indicators are 
#' added as additional columns to the dataframe.
#' 
#' @param HCR_df A data frame produced by eg [define_HCR_test_range] with 
#'   HCR parameters
#' @param HCR_MSE  HCR function produced by [create_HCR_MSE]
#' @return A data frame with the original parameters and additional columns 
#'   containing the performance indicators.
#' @export
#' 
run_HCR_MSE <- function(HCR_df, HCR_MSE) {
  HCR_df <- HCR_df |>
    dplyr::select(ID:ma) |>
    dplyr::mutate(Res = purrr::pmap(
      list(trIndex, trControl, change_limit, ma),
      \(trIndex, trControl, change_limit, ma) HCR_performance(HCR_MSE(trIndex, trControl, change_limit, ma)),
      .progress = "HCR Sim"
    )) |>
    tidyr::unnest(cols=Res)
  HCR_df <- evaluate_HCR(HCR_df, HCR_MSE)
  return(HCR_df)
}


#' Combines multiple HCR simulation results into a single tibble
#'
#' @inheritParams run_HCR_MSE
#' @param HCR_res A list of results from [HCR_performance]
#' @return A tibble containing combined results from the simulations
#' @export
#'   
combine_HCR_results <- function(HCR_res, HCR_df, HCR_MSE) {
  HCR_res <- dplyr::bind_rows(HCR_res)
  # Higher rank is better
  df <- HCR_res |>
    dplyr::group_by(ID) |>
    dplyr::summarise(Catch_Avg = mean(Catch_Avg),
                     Catch_Rng = mean(Catch_Rng),
                     Catch_pcile = min(Catch_pcile),  # this is not quite correct, but precautionary
                     lt_Blim = mean(lt_Blim),
                     at_Btar = mean(at_Btar),
                     gt_Btar = mean(gt_Btar),
                     Err_Type1 = mean(Err_Type1),
                     Err_Type2 = mean(Err_Type2)) |>
    dplyr::ungroup() |>
    dplyr::mutate(State = at_Btar - lt_Blim) |>
    dplyr::select(ID:gt_Btar, State, everything())
  
  HCR_res <- HCR_df |>
    dplyr::select(ID:ma) |>
    dplyr::left_join(df, by="ID") |>
    evaluate_HCR(HCR_MSE)
}


#' Ranks HCR based on performance indicators and evaluates whether the HCR is a 
#' candidate.
#'   
#' The resulting ranks and whether the HCR is a candidate or not is added as
#' columns to the data frame.
#' 
#' Ranks are based on catch (10th percentile) and stock state (proportion of 
#' the simulation spent around MSY level - that spent below Blim).
#' @inheritParams run_HCR_MSE
#' @return HCR data frame with added evaluation and ranks
#' @export
#' 
evaluate_HCR <- function(HCR_df, HCR_MSE) {
  ref_pt <- get("ref_pt", envir=environment(HCR_MSE))
  Ranks <- nrow(HCR_df)+1
  HCR_df$Catch_Rank <- Ranks - rank(HCR_df$Catch_pcile)   
  HCR_df$State_Rank <- Ranks - rank(HCR_df$State)    
  HCR_df <- HCR_df |>
    dplyr::mutate(Evaluation = 
                    dplyr::if_else((at_Btar >= ref_pt$mostly) & 
                                                (lt_Blim <= ref_pt$max_risk), 
                                   "Candidate", "Rejected")) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      Rank = max(Catch_Rank, State_Rank)) |>
    dplyr::ungroup()
  return(HCR_df)
}


#' Extracts MSY reference points from a JABBA model fit
#' 
#' JABBA uses the harvest rate, H, which is catch / biomass and limited to the 
#' range 0-1. For effort control, this is converted to a fishing mortality 
#' equivalent (H=1-exp(-F)), which avoids problems with boundary conditions.
#' 
#' @inheritParams run_HCR_MSE_list
#' @param ref_year Year to use for relative change in effort 
#' @return List of MSY reference points and a reference year for fishing 
#'   mortality for relative changes to that year
#' @export
#' 
MSY_refpt <- function(jabba_list, 
                      ref_year = NULL) {
  
  jabba_list <- convert_jabba2list(jabba_list)
  if (is.null(ref_year)) ref_year <- max(jabba_fit[[1]]$yr) 
  H_ref <- HMSY <- IMSY <- MSY <- m_ <- double(0)
  for (i in 1:length(jabba_list)) {
    Average_CPUE <- create_mean_CPUE(jabba_list[[i]]) # A function to combine multiple indices into a single index
    H_ref <- c(H_ref, dplyr::pull(dplyr::filter(jabba_list[[i]]$kbtrj, year==ref_year), H))
    HMSY <- c(HMSY, jabba_list[[i]]$pfunc$Hmsy)
    MSY <- c(MSY, jabba_list[[i]]$pfunc$MSY)
    IMSY = c(IMSY, Average_CPUE(CPUE_MSY(jabba_list[[i]])))
    m_ <- c(m_, jabba_list[[i]]$pars_posterior$m)
  }
  H_ref <- median(H_ref)
  IMSY <- mean(IMSY)
  HMSY <- median(HMSY)
  MSY <- median(MSY)
  m_ <- median(m_)
  BMSY <- m_^(-1 / (m_ - 1))
  
  return(list(BMSY = BMSY,         
              IMSY = IMSY,
              HMSY = HMSY,
              fMSY = -log(1-HMSY),
              ref_year = ref_year,
              f_ref = -log(1-H_ref), 
              MSY = MSY))
}


#' Defines standard risk-based reference points 
#' 
#' This maintains the reference points definition in one place and can be 
#' edited. The reference points are in a list and consist of:
#'   - BMSY = biomass at MSY as a proportion of the unexploited biomass
#'   - BMSY_range = accepted target range around BMSY
#'   - B_lim = limit reference point as a proportion of the unexploited biomass
#'   - mostly, likely, highly likely, high certainty = probabilistic definitions 
#'     of each word not currently used
#'   - max_risk = maximum acceptable risk for the candidate HCR selection. 
#'     Primarily used to exlcude HCR with this risk of falling below the limit
#'     reference point B_lim
#' 
#' @param B_tar Target reference point relative to MSY
#' @return List of reference points
#' @export
#' 
standard_risk_ref_pt <- function(B_tar = 1) {
  
  return(ref_pt <- list(
    B_tar = B_tar,  # Target relative to BMSY
    B_tar_range = B_tar*c(0.9, 1.2),
    B_lim = 0.5*B_tar,
    mostly = 0.5,
    likely = 0.7,
    highly_likely = 0.8,
    high_certainty = 0.95,
    max_risk = 0.1
  ))
}


#' Returns the median CPUE for each gear expected at MSY from a JABBA fit 
#'
#' @inheritParams create_HCR_MSE
#' @return A 1D array of each CPUE at MSY
#' @export
#' 
CPUE_MSY <- function(jabba_fit) {
  q <- jabba_fit$pars$Median[substr(rownames(jabba_fit$pars), 1, 1)=="q"]
  K <- jabba_fit$pars$Median[rownames(jabba_fit$pars)=="K"]
  m <- jabba_fit$pars$Median[rownames(jabba_fit$pars)=="m"]
  return(array(q * K * m^(-1/(m-1)), dim=c(1, length(q))))
}


#' Does minimal checks on jabba fit and ensures it is in a list 
#'
#' @inheritParams create_HCR_MSE
#' @return jabba_fit in a list 
#' @export
#' 
convert_jabba2list <- function(jabba_fit) {
  
  if ( ! is.null(jabba_fit$assessment)) jabba_fit <- list(jabba_fit) # turn jabba_fit into a list of fits if necessary
  
}


#' Run an HCR MSE based on a table of parameters and HCR MSE function over a 
#' list of JABBA fits 
#' 
#' The function runs the HCR MSE based on the set of HCR parameters in the 
#' dataframe and the provided MSE function on a list of JABBA fits and returns 
#' the combined results in a dataframe. This is the same as [run_HCR_MSE], but 
#' instead of a single JABBA fit, the function works for multiple alternative 
#' fits such as sensitivity analyses.
#' 
#' @inheritParams create_HCR_MSE
#' @inheritParams run_HCR_MSE
#' @param jabba_list a single jabba fit or a list of jabba fits
#' @param risk_rp Risk based reference levels defining targets relative to MSY 
#'   and maximum acceptable risks
#' @return A data frame with the original parameters and additional columns 
#'   containing the performance indicators.
#' @export
#' 
run_HCR_MSE_list <- function(
    jabba_list, 
    HCR_df, 
    control_type = "Effort",
    ProjLength=50, 
    nsim = 1000,
    ref_pt = standard_risk_ref_pt(),
    rho = 1) {

  jabba_list <- convert_jabba2list(jabba_list)  
  HCR_res <- list()

  # Each fit is run in the MSE HCR simulations

  for (i in 1:length(jabba_list)) {
    HCR <- create_HCR_MSE(jabba_list[[i]], 
                          control_type = control_type,
                          ProjLength = ProjLength, 
                          nsim = nsim,
                          ref_pt = ref_pt,
                          rho = rho)
    df <- run_HCR_MSE(HCR_df, HCR)
    df$fit <- names(jabba_list)[i]
    HCR_res[[i]] <- df
  }

  HCR_res <- combine_HCR_results(HCR_res, HCR_df, HCR)
  return(HCR_res)
}

