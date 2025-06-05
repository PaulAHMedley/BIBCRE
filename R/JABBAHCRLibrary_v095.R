 ##  ><> ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>
 ##  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>
 ##  ><> ##  ><> 
 ##  ><> ##    ><>     JABBA HCR Functions v0.95
 ##  ><> ##  ><> 
 ##  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>
 ##  ><> ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>  ><>


library("here")
library("dplyr")
library("tidyr")
library("ggplot2")
library("JABBA")
library("flextable")
#remotes::install_github("jabbamodel/JABBA")  

#' The functions below implement a simplified harvest control rule projection 
#' based on a fitted JABBA model. The method will work with some adaptation on 
#' any application of a Pella-Thomlinson MCMC (or related) fitted model to generate projected stock 
#' status over the projection period. It is configured for the Schaefer
#' model and the Pella-Thomlinson, so the `m` parameter is used. 
#' There is no checking that parameters are correctly specified in this version


#' Workflow
#' 
#' jabba_fits is a single fit or list of fits.
#' 
#' 1. Generate a jabba_fits
#' 2. Define the linear HCR, ma, change_limit ranges to test and generate a tibble of these values
#' 3. Run all HCR on JABBA fit info generating performance indicators
#' 4. Generate plots of HCR performance with rejected HCR in red.
#' 
  
# #### MSE function creation ####

#' Create a HCR MSE function 
#' 
#' The HCR function returns simulation results from HCR parameters using a 
#' jabba model fit results.
#' 
#' @param jabba_fit A list returned from package JABBA `fit_jabba` function
#' @param control_type Type of control applied by HCR either "Effort" or "Catch"
#' @param ProjLength   Projection length in years
#' @param nsim Number of simulations to run
#' @param ref_pt Specifies risk based reference points in a list. Otherwise 
#'   defaults apply.
#' @param rho Conversion factor from effort (f) to fishing mortality, so 
#'   H=(1-e(-rho*f)). Only used if control_type is "Effort"
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
                           rho = NULL) {
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
  
  if (is.null(rho)) rho <- -log(1-median(Bio$H))/1  #H=(1-e(-rho*f))  
  obserr <- sqrt(pull(jabba_fit$pars,
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
      slice(sn)
    
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
  r <- pull(Par, r)
  lsigma <- sqrt(dplyr::pull(Par, sigma2))
  if (NGears > 1L) 
    q <- as.matrix(dplyr::select(Par, tidyselect::starts_with("q."))) #Excludes auxiliary
  else
    q <- as.matrix(dplyr::select(Par, q)) #Excludes auxiliary
  
  if (jabba_fit$settings$model.type=="Schaefer") {
    prod_fun <- function(pBti) {
      pBti * (1 + r * (1 - pBti)) * rlnorm(nsim, 0, lsigma)}
  } else if (jabba_fit$settings$model.type %in% c("Fox", "Pella_m")) {
    m_1 <- pull(Par, m)-1
    r_m_1 <- r/m_1
    prod_fun <- function(pBti) {
      pBti * (1 + r_m_1 * (1 - pBti^m_1)) * rlnorm(nsim, 0, lsigma)}
  } else {
    stop("Error: Model type not recognised.")
  }
  # Tidy up
  rm(jabba_fit, ProjLength, ymin)
  
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
        H[ , pi] <- (1 - exp(-pjControl[, pi]*rho))
        C[ , pi] <- H[ , pi] * pB1
      }
      
      invalid <- (pB1 <= C[, pi])
      C[invalid, pi] <- 0
      pB[, ti+1L] <- pB1 - C[, pi]
      CPUE <- apply(q, 2, function(x) x*pB1*Par$K)
      pjIndex[, pi+1L] <- UpdateIndex(pjIndex[, pi], CPUE)
      pjControl[, pi+1L] <- CalcControl(pjIndex[, pi+1L], pjControl[, pi])
    } #ti
    
    return(list(pB=pB, C=C, H=H, Par=Par, Bio=Bio, 
                pvIndex=pvIndex, pvControl=pvControl, 
                pjIndex=pjIndex, pjControl=pjControl, 
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
                                  NInflex = 2, 
                                  NBreaks = 5, 
                                  change_limit = NA, 
                                  ma = 0.5) {
  Controls <- seq(rel_control_range[1], rel_control_range[2], length.out=NBreaks) * fMSY # Index values intervention points
  Indices <- seq(rel_index_range[1], rel_index_range[2], length.out=NBreaks) * IMSY  # Control levels to be applied
  
  trIndex <- list()
  trControl <- list()
  Ir <- integer(NInflex)
  Ir[] <- 1L

  while (TRUE) {
    for (j in Ir[NInflex-1L]:NBreaks) {
      Ir[NInflex] <- j
      trIndex <- append(trIndex, list(Indices[Ir]))
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
  
  return(tidyr::expand_grid(trIndex, trControl, change_limit, ma) |> 
         dplyr::mutate(ID = dplyr::row_number()) |> 
         dplyr::select(ID, dplyr::everything()))
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
create_empirical_index <- function(ma, func=NULL, se) {
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
    rm(q)
    return(    
      NULL
      # function(vCPUE) {
      #   return(vCPUE)
      # }
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

#' Calculate the performance measures of a harvest control rule and return it 
#' as a row in a data frame
#' 
#' The performance measures consist of the 1) average catch, 2) the mean catch 
#' range around the average (mean(abs(catch-mean_catch))), 3) the HCR Type 1 
#' error rate, 4) the HCR Type 2 error rate, 5) the proportion of time that the 
#' stock is below 0.5 BMSY (LRP) and 6) the proportion of time that the stock 
#' is outside 90-110% BMSY.
#' 
#' The HCR error rates are scores indicating how well the HCR is responding. 
#' For Type 1, When the stock is below the limit reference point and the index 
#' is below the upper trigger but above the lower trigger, it scores one, 
#' whereas if it is above both triggers it scores two. The mean of this value 
#' is taken. Similarly for Type 2, responding by reducing harvest when the 
#' stock is at MSY scores two when the reduction is minimised, one otherwise.
#' The higher the scores are, the worse the performance.
#' 
#' @param HCR_sim Fishery HCR simulation result list from the MSE 
#'   (see [create_HCR_MSE]).
#' @return Tibble row of performance measures
#' @export
#' 
HCR_performance <- function(HCR_sim) {
  BMSY_range <- HCR_sim$ref_pt$BMSY_range #c(0.90, 1.2)
  B_MSY <- HCR_sim$ref_pt$BMSY
  B_lim <- HCR_sim$ref_pt$B_lim # 0.5
  MaxRisk <- HCR_sim$ref_pt$max_risk
  Mostly <- HCR_sim$ref_pt$mostly
  
  BB0 <- as.vector(HCR_sim$pB[ , (HCR_sim$HCR$TN+1):(HCR_sim$HCR$PTN)])   #For Schaefer
  Index <- as.vector(HCR_sim$pjIndex[ , -1])
  
  Catch <- apply(HCR_sim$C, 2, function(x) x*HCR_sim$Par$K)
  Catch_Avg <- mean(Catch)
  Catch_Rng <- mean(abs(Catch-Catch_Avg))
  Catch_pcile <- quantile(Catch, probs=MaxRisk)
  lt_Blim <- sum(BB0 < B_lim)/length(BB0)
  at_BMSY <- sum((BB0 >= BMSY_range[1]) & (BB0 <= BMSY_range[2])) / length(BB0)
  gt_BMSY <- sum(BB0 > B_MSY)/length(BB0)
  # Management response when not necessary
  Err_Type1 <- 0.5*mean(((BB0 > B_lim) & (Index < HCR_sim$HCR$trIndex[1])) +
                          ((BB0 > BMSY_range[1]) & (Index < HCR_sim$HCR$trIndex[2])))
  # No management response when it is necessary
  Err_Type2 <- 0.5*mean(((BB0 < B_lim) & (Index > HCR_sim$HCR$trIndex[1])) +
                          ((BB0 < B_lim) & (Index > HCR_sim$HCR$trIndex[2])))
  
  return(tibble(
    Catch_Avg = Catch_Avg, 
    Catch_Rng = Catch_Rng,
    Catch_pcile = Catch_pcile,
    lt_Blim = lt_Blim,
    at_BMSY = at_BMSY,
    gt_BMSY = gt_BMSY,
    State = at_BMSY - lt_Blim,
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
  HCR_res <- bind_rows(HCR_res)
  # Higher rank is better
  df <- HCR_res |>
    dplyr::group_by(ID) |>
    dplyr::summarise(Catch_Avg = mean(Catch_Avg),
                     Catch_Rng = mean(Catch_Rng),
                     Catch_pcile = min(Catch_pcile),  # this is not quite correct, but precautionary
                     lt_Blim = mean(lt_Blim),
                     at_BMSY = mean(at_BMSY),
                     Err_Type1 = mean(Err_Type1),
                     Err_Type2 = mean(Err_Type2)) |>
    dplyr::ungroup() |>
    dplyr::mutate(State = at_BMSY - lt_Blim) |>
    dplyr::select(ID:at_BMSY, State, everything())

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
    dplyr::mutate(Evaluation = if_else((at_BMSY >= ref_pt$mostly) & (lt_Blim <= ref_pt$max_risk), "Candidate", "Rejected")) |>
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
#' @inheritParams create_HCR_MSE
#' @param ref_year Year to use for relative change in effort 
#' @return List of MSY reference points
#' @export
#' 
MSY_refpt <- function(jabba_fit, ref_year=NULL) {
  
  if ( ! is.null(jabba_fit$assessment)) jabba_fit <- list(jabba_fit) # turn jabba_fit into a list of fits if necessary
  if (is.null(ref_year)) ref_year <- max(jabba_fit[[1]]$yr) 
  H_ref <- HMSY <- IMSY <- MSY <- double(0)
  for (i in 1:length(jabba_fit)) {
    print(i)
    Average_CPUE <- create_mean_CPUE(jabba_fit[[i]]) # A function to combine multiple indices into a single index
    H_ref <- c(H_ref, dplyr::pull(filter(jabba_fit[[i]]$kbtrj, year==ref_year), H))
    HMSY <- c(HMSY, jabba_fit[[i]]$pfunc$Hmsy)
    MSY <- c(MSY, jabba_fit[[i]]$pfunc$MSY)
    IMSY = c(IMSY, Average_CPUE(CPUE_MSY(jabba_fit[[i]])))
  }
  H_ref <- median(H_ref)
  IMSY <- mean(IMSY)
  HMSY <- median(HMSY)
  MSY <- median(MSY)
  
  risk_rp <- standard_risk_ref_pt()
  
  return(list(BMSY = risk_rp$BMSY,         
              Blim = risk_rp$B_lim,     # 50%BMSY
              IMSY = IMSY,
              HMSY = HMSY,
              fMSY = -log(1-HMSY),
              ref_year = ref_year,
              f_ref = -log(1-H_ref), 
              MSY = MSY))
}


#' Defines standard risk-based reference points.
#' 
#' This is purely to maintain reference points definition in one place. 
#' Can be edited.
#' 
#' @return List of reference points
#' @export
#' 
standard_risk_ref_pt <- function(m=2) {
  BMSY <- m^(-1/(m-1))
  return(ref_pt <- list(
      BMSY = BMSY,  
      BMSY_range = BMSY*c(0.9, 1.2),
      B_lim = 0.5*BMSY,
      mostly = 0.5,
      likely = 0.7,
      highly_likely = 0.8,
      high_certainty = 0.95,
      max_risk = 0.1
    ))
}


#' Gets the median CPUE expected at MSY from a JABBA fit
#'
#' @inheritParams create_HCR_MSE
#' @return A 1D array of each CPUE at MSY
#' @export
#' 
CPUE_MSY <- function(jabba_fit) {
  q <- jabba_fit$pars$Median[substr(rownames(jabba_fit$pars), 1, 1)=="q"]
  K <- jabba_fit$pars$Median[rownames(jabba_fit$pars)=="K"]
  return(array(0.5*q*K, dim=c(1, length(q))))  #Schaefer
}


# #### Individual Sims: Graphs and Tables ####


#' Plot the stock status and harvest rate for a single HCR test simulation.
#'   
#' A ggplot is returned with the biomass as a proportion of BMSY and H as a 
#' proportion of HMSY over time.
#'   
#' @inheritParams HCR_performance
#' @return A probability density plot of biomass dynamics projection under a HCR
#' @export
#' 
graph_sim_BMSY_FMSY <- function(HCR_sim) {

  ## Fishing Mortality
  catch_mat <- with(HCR_sim, matrix(Bio$Catch, nrow=HCR$nsim, ncol=HCR$TN))
  invHMSY <- with(HCR_sim, 4 / (Par$r*Par$K))
  obHMSY <- with(HCR_sim, apply(catch_mat, 2, function(x) x*invHMSY))
  invHMSY <- with(HCR_sim, 4 / Par$r)     # At MSY C/B = 0.5 r for Schaefer
  pjHMSY  <- with(HCR_sim, 
                  apply(C, 2, function(x) x*invHMSY))
  HMSY  <- as.vector(cbind(obHMSY, pjHMSY))
  rsl_H <- with(HCR_sim, tibble::tibble(year = rep(Bio$year[1]:(Bio$year[1] + HCR$PTN-1), each=HCR$nsim),
                  Var = "HMSY", Val = HMSY))
  
  ## Biomass    
  rsl_B <- with(HCR_sim, tibble::tibble(year = rep(Bio$year[1]:(Bio$year[1] + HCR$PTN), each=HCR$nsim),
                                      Var = "BMSY",
                                      Val = 2*as.vector(pB)))
  
  rbind(rsl_B, rsl_H) |>
    dplyr::filter(Val > 0, Val < 2.4, !is.na(Val)) |>
    ggplot2::ggplot(aes(x=year, y=Val)) +
    ggplot2::geom_bin2d(bins=HCR_sim$HCR$PTN-2) +
    #geom_line(aes(colour=Type)) +
    ggplot2::scale_y_continuous(limits = c(0, 2.4))+
    ggplot2::geom_vline(xintercept=(HCR_sim$Bio$year[1] + HCR_sim$HCR$TN-1))+
    ggplot2::geom_hline(yintercept=1.0)+
    ggplot2::geom_hline(yintercept=0.5)+
    ggplot2::scale_fill_gradient(low = "#96B1F7", high = "#030B43") +
    ggplot2::facet_wrap(vars(Var), nrow=2) +
    ggplot2::labs(x="Year", y="Biomass / harvest rate relative to MSY")
  
}


#' Plot the catches for a single HCR test simulation.
#'
#' @inheritParams HCR_performance
#' @return graph showing binned probability for observed and simulated catch 
#'   time series
#' @export
#' 
graph_sim_catch <- function(HCR_sim) {
  catch_df <- dplyr::select(HCR_sim$Bio, year, Catch)

  catch_pj <- with(HCR_sim, apply(C, 2, function(x) x*Par$K))
  rsl_C <- with(HCR_sim, 
                tibble::tibble(year = rep((Bio$year[1] + HCR$TN):(Bio$year[1] + HCR$PTN-1), each=HCR$nsim),
                               Catch = as.vector(catch_pj))) |>
           rbind(catch_df)
  
  ggplot2::ggplot(rsl_C, aes(x=year, y=Catch)) +
    ggplot2::geom_bin2d(bins=HCR_sim$HCR$PTN-2) +
    ggplot2::scale_fill_gradient(low = "#96B1F7", high = "#030B43") +
    ggplot2::labs(x="Year", y="Catch(t)")
} 


#' Table the proportion of each state-response for HCR decisions
#'
#' @inheritParams HCR_performance
#' @return flextable object containing the state-response proportions for the
#'   HCR
#' @export
#' 
table_sim_decision <- function(HCR_sim) {
  Lim <- HCR_sim$ref_pt$B_lim
  Trig <- HCR_sim$ref_pt$BMSY_range[1]
  lvlBt <- c("B<Limit",   #Limit
             "B<Trigger", #Trigger
             "B~Target")  #Target
  lvlIt <- c("It<Limit",  #Limit
             "It<Trigger", #Trigger
             "It~Target")  #Target
  BB0 <- as.vector(HCR_sim$pB[ , (HCR_sim$HCR$TN+1):(HCR_sim$HCR$PTN)])
  Index <- as.vector(HCR_sim$pjIndex[ , -1])
  HCR_r <- tibble::tibble(State=lvlBt[3 - ((BB0 < Lim) + (BB0 < Trig))],
                    Response=lvlIt[3 - ((Index < HCR_sim$HCR$trIndex[1]) + 
                                      (Index < HCR_sim$HCR$trIndex[2]))]) |>
    # mutate(Response=ifelse(Response %in% c("It~Target"), "It=Target", Response),
    #        State = ifelse(State %in% c("B<Target", "B>Target"), "B>=Target", State)) |>
    dplyr::mutate(Response = factor(Response, levels=lvlIt),
           State = factor(State, levels = lvlBt))
  
  tbl <- tibble::as_tibble(table(HCR_r))
  tbl$n <- tbl$n/sum(tbl$n)
  tbl <- tbl |> 
    tidyr::pivot_wider(names_from=State, values_from=n)
  
  tbl <- flextable::flextable(tbl) |>
    flextable::color(i=3, j=2:3, color="red", part="body") |>
    flextable::color(i=1:2, j=4, color="blue", part="body") |>
    flextable::colformat_double(digits=3)
  return(tbl)
}


#' Table the simulation proportion of time that the stock spends in each status 
#' interval.
#'
#' @inheritParams HCR_performance
#' @return `flextable` containing the simulations status proportions for the HCR
#' @export
#' 
table_sim_status <- function(HCR_sim) {
  std_brk <- c(0.0, 0.5, 0.9, 1.2, 100.0) # Appropriate for the Schaefer model
  c_std_brk <- paste0(format(std_brk[-length(std_brk)],nsmall=1), "-",
                      format(std_brk[-1], nsmall=1))
  c_std_brk[length(c_std_brk)] <- paste0(">", format(std_brk[length(std_brk)-1L],nsmall=1))
  
  invHMSY <- with(HCR_sim, 2 / rep(Par$r, HCR$PTN-HCR$TN))

  # Projection statistics
  BMSY <- with(HCR_sim, 2 * as.vector(pB[ ,(HCR$TN+1):HCR$PTN]))
  HMSY  <- with(HCR_sim, as.vector(C) * invHMSY)
  
  cnt <- tibble(
    `B/BMSY` = c_std_brk,
    `B%` = hist(BMSY, breaks=std_brk, plot=FALSE)$counts,
    `H%` = hist(HMSY, breaks=std_brk, plot=FALSE)$counts,
    #IMSY = hist(as.vector(HCR_sim$pjIndex[ , -1])*invIMSY, breaks=std_brk, plot=FALSE)$counts
  )
  cnt$`B%` <- 100*cnt$`B%`/sum(cnt$`B%`)
  cnt$`H%` <- 100*cnt$`H%`/sum(cnt$`H%`)
  return(
    flextable::flextable(cnt) |>
      flextable::colformat_double(digits=1)
  )
}


#' Table the standard performance measures for a HCR simulation
#'
#' @inheritParams HCR_performance
#' @return `flextable` containing the HCR performance measures
#' @export
#' 
table_sim_performance <- function(HCR_sim) {
  tbl <- HCR_performance(HCR_sim)
  
  tbl |>
    flextable::flextable() |>
    flextable::set_header_labels(values = c("Mean Catch", "Mean Catch Range", "Catch Lower Percentile", 
                                            "Below LRP", "In Target Range", "Above Target", "Good Fishery State", 
                                            "HCR False Positive", "HCR False Negative")) |>
    flextable::colformat_double(j=1:3, digits=0, big.mark="") |>
    flextable::colformat_double(j=4:9, digits=3, big.mark="")
}

  
# #### HCR tests: Graphs and Tables ####


#' Plot harvest control rules with the stock status index on the x-axis and 
#' control on the y-axis.
#'   
#' A ggplot is returned with the HCR plotted as a line between 
#' the index of stock status and control.
#'   
#' @inheritParams run_HCR_MSE
#' @inheritParams graph_BMSY_FMSY
#' @param HCR_ID Logical - whether to include the ID as a factor in the plot. 
#'   Will only work for 12 or fewer HCR.
#' @return A ggplot object plotting the HCR's
#' @export
#' 
graph_linear_HCR <- function(HCR_df,
                             HCR_sim = NULL,
                             HCR_ID = TRUE) {
  line_df <- HCR_df |>
    dplyr::select(ID, trIndex, trControl) |>
    tidyr::unnest(c(trIndex, trControl)) |>
    dplyr::mutate(ID = factor(ID))
  
  maxIndex <- 1.2 * max(line_df$trIndex)
  
  lo_df <- line_df |>
    dplyr::group_by(ID) |>
    dplyr::summarise(trIndex = 0, trControl = min(trControl)) |>
    dplyr::ungroup()
  hi_df <- line_df |>
    dplyr::group_by(ID) |>
    dplyr::summarise(trIndex = maxIndex, trControl = max(trControl)) |>
    dplyr::ungroup()
  
  line_df <- dplyr::bind_rows(line_df, lo_df, hi_df)
  
  gp <- ggplot2::ggplot(line_df, aes(x = trIndex, y = trControl, group = ID)) +
    ggplot2::geom_line() +
    ggplot2::labs(y = "Control", x = "HCR Index") +
    ggplot2::coord_cartesian( y = c(0, NA))

  if (HCR_ID & nrow(HCR_df) <= 12) gp <- gp + facet_wrap(vars(ID))
  
  if (!is.null(HCR_sim)) {
    pv_df <- tibble::tibble(pvIndex = pvIndex, pvControl = pvControl)
    tr_df <- tibble::tibble(trIndex = trIndex, trControl = trControl)
    gp <- gp +
      ggplot2::geom_point(pv_df, aes(x = pvIndex, y = pvControl)) +
      ggplot2::geom_line(tr_df, aes(x = pvIndex, y = pvControl), color = "blue")
  }
  return(gp)
}


#' Plot the average catch and mean catch range.
#'   
#' A ggplot is returned with the average catch and mean catch range on the 
#' x- and y-axes respectively. Candidate HCR are plotted as black points and 
#' rejected HCR as red points. A linear fit is shown going through the origin
#' indicating the exchange between increasing catch with increasing catch range. 
#' HCR below this line indicate a better exchange rate between these two
#' measures.
#'   
#' @inheritParams graph_linear_HCR
#' @return A ggplot object plotting the HCR performance catches
#' @export
#' 
graph_HCR_catches <- function(HCR_df) {
  
  # Add dotted lines for means
  
  Mean_Catch <- mean(HCR_df$Catch_Avg)
  Var_Catch <- mean(HCR_df$Catch_Rng)
  
  ggplot2::ggplot(HCR_df, aes(x=Catch_Avg, y=Catch_Rng)) +
    ggplot2::geom_point(aes(color=Evaluation)) +
    ggplot2::scale_color_manual(values = c("Candidate" = "black", "Rejected" = "red")) +
    ggplot2::geom_hline(yintercept=Var_Catch, linetype="dotted") +
    ggplot2::geom_vline(xintercept=Mean_Catch, linetype="dotted") +
    ggplot2::geom_smooth(method="lm", formula = y~x+0, se=FALSE, linetype="solid", alpha=0.5)
    
}


#' Plot the HCR stock status relative to the BMSY (TRP) and Blim (LRP)
#'   
#' A ggplot is returned with the the proportion of the simulation time across 
#' all runs for an HCR below the limit (x-axis) and around the target (y-axis). 
#' Candidate HCR are plotted as black points and rejected HCR as red points. 
#' Risk based reference points are indicated by dotted vertical and horizontal 
#' lines, which are set at default values if no alternative candidate 
#' definition. is provided.
#'   
#' @inheritParams evaluate_HCR

#' @return A ggplot object plotting the HCR performance catches
#' @export
#' 
graph_HCR_status <- function(HCR_df, HCR_MSE) {
  
  ref_pt <- get("ref_pt", envir=environment(HCR_MSE))
  ggplot2::ggplot(HCR_df, aes(x=lt_Blim, y=at_BMSY)) +
    ggplot2::geom_point(aes(color=Evaluation)) +
    ggplot2::scale_color_manual(values = c("Candidate" = "black", "Rejected" = "red")) +
    ggplot2::geom_hline(yintercept=ref_pt$mostly, linetype="dotted") +
    ggplot2::geom_vline(xintercept=ref_pt$max_risk, linetype="dotted") +
    ggplot2::labs(x = "Proportion Below B_lim", y = "Proportion at MSY")
  
}

#' Plot the HCR stock status relative to the BMSY (TRP) and Blim (LRP)
#'   
#' A ggplot is returned with the the proportion of the simulation time across 
#' all runs for an HCR below the limit (x-axis) and around the target (y-axis). 
#' Candidate HCR are plotted as black points and rejected HCR as red points. 
#' Risk based reference points are indicated by dotted vertical and horizontal 
#' lines, which are set at default values if no alternative candidate 
#' definition. is provided.
#'   
#' @inheritParams graph_linear_HCR
#' @return A ggplot object plotting the HCR performance catches
#' @export
#' 
graph_HCR_status_catches <- function(HCR_df) {
    
  HCR_df |>
    ggplot2::ggplot(aes(x=State, y=Catch_pcile, colour=factor(change_limit))) +
    ggplot2::geom_point() +
    ggplot2::labs(x="Stock State", y="Catch - lower 10%_ile", color = "Change Limit") 
}


#' Plot the HCR decision performance in relation to mistakes made by the 
#' decision rule.
#'   
#' A ggplot is returned with the the proportion of Type 1 (false positive) and 
#' Type 2 (false negative) errors. For Type 1, the decision rule results in 
#' management intervention reducing harvest when the stock in reality is at or 
#' above its target level. For Type 2 (false negative), the decision rule does 
#' not harvest when the stock is in reality below the limit reference point, 
#' which in this case is considered worse as it is less precautionary. Higher 
#' error rates indicate poorer HCR performance and might be improved by 
#' adjusting the index trigger points.
#'  
#' Candidate HCR are plotted as black points and rejected HCR as red points. 
#' Risk based reference points are indicated by dotted vertical and horizontal 
#' lines, which are set at default values if no alternative candidate 
#' definition. is provided.
#'   
#' @inheritParams graph_linear_HCR
#' @return A ggplot object plotting the HCR performance catches
#' @export
#' 
graph_HCR_decision <- function(HCR_df, ref_pt=NULL) {
  
  # Add dotted lines for means
  if (is.null(ref_pt)) {
    # defaults
    MaxRisk <- 0.1
  } else {
    MaxRisk <- ref_pt$max_risk
  }
  
  ggplot2::ggplot(HCR_df, aes(x=Err_Type1, y=Err_Type2)) +
    ggplot2::geom_point(aes(color=Evaluation)) +
    ggplot2::scale_color_manual(values = c("Candidate" = "black", "Rejected" = "red")) +
    ggplot2::geom_hline(yintercept=MaxRisk, linetype="dotted") +
    ggplot2::geom_vline(xintercept=MaxRisk, linetype="dotted") +
    ggplot2::labs(x = "Type 1 (False Positive)", y = "Type 2 (False Negative)")
}



#' Table the standard performance measures for a set of candidate HCRs
#'
#' @inheritParams graph_linear_HCR
#' @return A `flextable` containing the HCR performance measures
#' @export
#' 
table_HCR_performance <- function(HCR_df) {
  HCR_df |>
    dplyr::select(ID, change_limit:Err_Type2, State, Catch_Rank, State_Rank) |>
    flextable::flextable() |>
    flextable::set_header_labels(values = c("ID", "Change Limit", "Index Smoother",
                                            "Mean Catch", "Mean Catch Range", "Catch Lower Percentile", 
                                            "Below LRP", "In Target Range", "Above TRP", 
                                            "Good Fishery State", 
                                            "HCR False Positive", "HCR False Negative", 
                                            "Catch Rank", "State Rank")) |>
    flextable::colformat_int(j=1, big.mark="") |>
    flextable::colformat_double(j=4:6, digits=0, big.mark="") |>
    flextable::colformat_double(j=7:11, digits=4) |>
    flextable::colformat_double(j=12:13, digits=0, big.mark="") |>
    flextable::autofit() 
}


# #### Miscellaneous functions ####


#' Plot correlation between simulated biomass and index using different moving
#' average parameters. 
#'   
#' A ggplot is returned with the correlation for different moving average 
#' parameters over range 0-1. Higher correlation is better. When the moving 
#' average is zero, it will not change, whereas 1 indicates previous index 
#' values are ignored. The moving average reduces the index noise. This 
#' graph may not show a clear optimum, in which case a value around 0.5 is 
#' recommended
#'   
#' @inheritParams create_HCR_MSE
#' @return A ggplot plotting moving average parameter against index/biomass 
#'   correlation
#' @export
#' 
graph_ma <- function(jabba_fit) {
  minstatus <- 0.005 
  nsim <- 1000L
  PN <- 40
  
  Average_CPUE <- Create_MeanCPUEFunction(jabba_fit) # A function to combine multiple indices into a single index
  IMSY <- Average_CPUE(CPUE_MSY(jabba_fit))
  HMSY <- jabba_fit$refpts$fmsy[1]
  obserr <- sqrt(pull(jabba_fit$pars,
                     Median)[startsWith(rownames(jabba_fit$pars), "tau2")])
  
  # Extract the data components that will be used from the fit
  Par <- jabba_fit$pars_posterior 
  Bio <- jabba_fit$kbtrj |>
    dplyr::select(year, iter, harvest:BB0)
  ymin <- min(Bio$year) - 1L
  
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
  }
  
  Bio <- dplyr::mutate(Bio, tim = as.integer(year - ymin))
  
  # Dimensions
  TN <- max(Bio$year) - ymin
  PTN <- PN + TN
  
  Index <- double(TN)
  pB <- matrix(0, ncol=PTN, nrow=nsim)
  pB[cbind(Bio$iter, Bio$tim)] <- Bio$BB0   # extract B/B0 estimates as matrix
  NGear <- ncol(jabba_fit$inputseries$cpue) - 1L
  Catch <- jabba_fit$inputseries$catch$catch
  CPUE <- as.matrix(jabba_fit$inputseries$cpue[,-1])
  UpdateIndex <- Create_EmpiricalIndexFunction(0.5, Average_CPUE, obserr)
  Index[1] <- mean(CPUE, na.rm=T)
  for (ti in 2:TN) {
    Index[ti] <- UpdateIndex(Index[ti-1L], CPUE[ti-1L,])
  }
  
  pjIndex <- matrix(0, nrow=nsim, ncol=(PN+1))
  pjIndex[,1] <- Index[TN]
  C <- matrix(0.0, nrow=nsim, ncol=PN)
  H <- C
  
  # Parameters   
  Hexp <- HMSY*0.9
  r <- pull(Par, r)
  lsigma <- sqrt(pull(Par, sigma2))
  q <- as.matrix(dplyr::select(Par, starts_with("q.")))
  
  for (pi in 1:PN) {
    ti <- TN+pi-1
    pB1 <- (pB[, ti] * (1 + r * (1 - pB[, ti])) * rlnorm(nsim, 0, lsigma))
    pB1[pB1 < minstatus] <- minstatus
    
    H[ , pi] <- Hexp*rlnorm(nsim, 0, 0.1)
    C[ , pi] <- H[ , pi] * pB1
    
    invalid <- (pB1 <= C[, pi])
    C[invalid, pi] <- 0
    pB[, ti+1L] <- pB1 - C[, pi]
  } #ti
  
  B <- apply(pB[,(TN):(PTN)], 2, function(x) x*Par$K)
  
  cor_func <- function(ma){
    UpdateIndex <- Create_EmpiricalIndexFunction(ma, Average_CPUE, obserr)
    for (pi in 1:(PN)) {
      CPUE <- apply(q, 2, function(x) x*B[, pi])
      pjIndex[, pi+1L] <- UpdateIndex(pjIndex[, pi], CPUE)
    }
    return(cor(as.vector(B), as.vector(pjIndex), use="complete.obs"))
  }
  
  bm_ma <- double(20)
  p_ma <- seq(0, 1.0, length.out=20)
  for (i in seq(p_ma)) {
    bm_ma[i] <- cor_func(p_ma[i]) 
  }
  
  tibble::tibble(x = p_ma, y = bm_ma) |>
    ggplot2::ggplot(aes(x=x, y=y)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::labs(x="HCR Moving Average Parameter", y="Correlation (B : HCR Index)")
}



# #### Database Access tests ####

# This code needs testing and will only work on a Access database that is local network or machine


#' Test procedure to pull data directly from MS Access using RODBC in Windows. 
#'   
#' This is only a test procedure that would use a query maintained in MS Access
#' to extract a table of data suitable for the stock assessment analyses. The 
#' query would need to be set up in the database. The records would be 
#' filtered by the species name. Separate functions should be created for each 
#' data type (catch, catch-effort, length frequency).
#' 
#' As well as this function, there should be a manual procedure to save data in 
#' Excel files for manual transfer. This will likely be necessary if the 
#' database is not on a LAN.
#'   
#' @param SpName The species name to filter the data
#' @return A data frame of the required data
#' @export
#' 
GetCatchEffortDataQuery <- function(SpName) {
  # This works on more recent versions of 
  ## Set up driver info and database path - note hard address for the database location
  DRIVERINFO <- "Driver={Microsoft Access Driver (*.mdb, *.accdb)};"
  MDBPATH <- "D:/paula/Dropbox/Americas/Seabob/DataCollection/GuyanaSeabob/GSeabobData_20240203.mdb"
  PATH <- paste0(DRIVERINFO, "DBQ=", MDBPATH)
  qry <- "EffortByYear"  #"qryCatchEffort"
  
  ## Establish connection
  channel <- RODBC::odbcDriverConnect(PATH)

  ## Load data into R dataframe
  # df <- sqlQuery(channel,
  #                "SELECT [Year], [CountOfID], [PlantID] FROM [EffortByYear] WHERE [Year] > 2006;",
  #                stringsAsFactors = FALSE)
  
  # This is not necessary, but allows for changes in column names. 
  # "SELECT * " could also be tried, but some users have reported problems
  qry_vars <- RODBC::sqlColumns(channel, qry)["COLUMN_NAME"]
  ## Add brackets to each variable (ie. [variable]) to maintain ACCESS syntax
  qry_vars$COLUMN_NAME <- paste0("[", qry_vars$COLUMN_NAME, "]")
  ## Transform dataframe column into string separated by comma
  cols <- paste0(qry_vars[1:nrow(qry_vars), ], collapse=",")
  ## Create ORDER BY string

  ## Extract table of interest as dataframe
  df <- RODBC::sqlQuery(channel,
                        paste0("SELECT ", cols, " FROM [", qry, "] WHERE Species = ", SpName, ";"),
               stringsAsFactors=FALSE)
  
  ## Close and remove channel
  RODBC::odbcClose(channel)
  rm(channel)
  return(df)
}
  

