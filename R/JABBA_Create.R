
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
#' @param proj_length   Projection length in years
#' @param nsim Number of simulations to run
#' @param ref_pt Specifies risk based reference points in a list. Otherwise 
#'   defaults apply.
#' @return A function that will apply the HCR using the JABBA model fit and 
#'   specified HCR parameters: index-control inflection point vectors, the 
#'   control change limit and moving average parameter.
#' @export
#' 
create_JABBA_MSE <- function(jabba_fit, 
                           proj_length=50, 
                           nsim = 1000,
                           ref_pt = standard_risk_ref_pt()) {
  stock_assessment <- "JABBA" 
  # candidate HCR definitions: not used in this function, but recorded for later evaluation
  # can be obtained from the function using "get"
  # Reference points are relative to B0 for Schaefer
  minstatus <- 0.005 
  NGears <- sum(substr(rownames(jabba_fit$pars), 1, 2)=="q.")
  if (NGears==0L) NGears <- 1L
  # Extract the data components that will be used from the fit
  
  # Standardise names to JABBAstan
  Par <- jabba_fit$pars_posterior 
  names(Par)[1] <- "Binf"
  names(Par)[4] <- "P0"
  names(Par)[length(Par)] <- "Nus"
  Par$Nus <- sqrt(Par$Nus)
  Par$BMSY <- jabba_fit$refpts_posterior$BmsyK
  
  Bio <- jabba_fit$kbtrj |>
    dplyr::select(year, iter, harvest:BB0) |>
    dplyr::mutate(Ft = -log(1-H))
  
  Par <- Par |>
    dplyr::mutate(
      B_BMSY = Bio |> 
        dplyr::slice_max(year, n = 1) |> 
        dplyr::pull(BB0) / BMSY
    )
  
  start_year <- min(Bio$year)
  # Dimensions
  TN <- as.integer(max(Bio$year) - start_year + 1)
  PN <- as.integer(proj_length)
  PTN <- PN + TN
  
  Avg_CPUE <- JABBA_mean_CPUE(jabba_fit)
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
  Bio <- dplyr::mutate(Bio, tim = as.integer(year - start_year + 1))
  
  pB <- matrix(0, ncol=PTN+1L, nrow=nsim)
  pB[cbind(Bio$iter, Bio$tim)] <- Bio$BB0   # extract B/B0 estimates as matrix
  Ft <- C <- matrix(0.0, nrow=nsim, ncol=PTN)
  C[cbind(Bio$iter, Bio$tim)] <- Bio$Catch
  Ft[cbind(Bio$iter, Bio$tim)] <- Bio$Ft
  pvCPUE <- as.matrix(jabba_fit$inputseries$cpue[,-1])
  
  # Model parameters   
  r <- dplyr::pull(Par, r)
  Binf <- dplyr::pull(Par, Binf)
  lsigma <- dplyr::pull(Par, Nus)
  B_trial <- apply(sweep(pB[, 1:TN], MARGIN=1, STATS=as.array(Binf), FUN = "*"), 1, min)
  
  if (NGears > 1L) 
    q <- as.matrix(dplyr::select(Par, tidyselect::starts_with("q."))) #Excludes auxiliary
  else
    q <- as.matrix(dplyr::select(Par, q)) #Excludes auxiliary
  
  if (jabba_fit$settings$model.type=="Schaefer") {
    m_ <- 2
    m_1 <- m_ - 1
    r_m_1 <- r/(m_ - 1)
    prod_fun <- function(pBti) {
      pBti * (1 + r * (1 - pBti)) * rlnorm(1, 0, lsigma)}
  } else if (jabba_fit$settings$model.type %in% c("Fox", "Pella_m")) {
    m_ <- dplyr::pull(Par, m)
    m_1 <- m_ - 1
    r_m_1 <- r/(m_ - 1)
    prod_fun <- function(pBti) {
      pBti * (1 + r_m_1 * (1 - pBti^m_1)) * rlnorm(1, 0, lsigma)}
  } else {
    stop("Error: Model type not recognised.")
  }
  dat <- list(YR = start_year,
              TN = nrow(jabba_fit$inputseries$catch),
              TCA_ca = jabba_fit$inputseries$catch$catch,
              Catch_cv = jabba_fit$settings$catch.cv,
              inputseries = jabba_fit$inputseries
              )
  ref_pt$B_tar <- dplyr::pull(Par, BMSY)*ref_pt$TRP
  ref_pt$B_lim <- ref_pt$B_tar*ref_pt$LRP
  ref_pt$F_tar <- with(Par, - log((m-1) / (m - 1 + r * (1-ref_pt$B_tar^(m-1)))))

  #F_tar <- with(Par, - log((m-1) / (m - 1 + r * (1-tmp1$ref_pt$B_tar^(m-1)))))
  #Btar <- with(tmp1$Par, tmp1$ref_pt$B_tar*(1 + r * (1-tmp1$ref_pt$B_tar^(m-1)) / (m-1))*(exp(-F_tar)))
  
  
  ref_pt$rp_type <- "MSY"
  Par <- dplyr::select(Par, -BMSY)
  # "YR"       "TN"       "TCA_ca"   "TCEN"     "TCE_t"    "TCE_ca"   "TCE_ef"   "polP0m"   "polP0s"   "polrm"    "polrs"   
  # "polMm"    "polMs"    "polqm"    "polqs"    "polBinfm" "polBinfs" "poNus"    "poNuss"   "poce_cvm" "poce_cvs" "Catch_cv"  # Tidy up
  rm(jabba_fit, proj_length, m_, Bio)
  
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
      F_trial <- C_trial
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
          F_limit <- Ft1*(1 - ctrl_pF[cj]) + pjControl[ , cj, pi]
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

#' Collapses multiple CPUE indices into single weighted geometric mean index
#'
#' A function is returned taking a matrix of CPUE with columns for each
#' CPUE index and convert to a single vector index. The indices are assumed to
#' be log-normally distributed, so the geometric mean is taken weighted by the
#' median JABBA observation error estimates (tau2).
#'   
#' @inheritParams create_JABBA_MSE
#' @return A function taking a matrix with CPUE in columns and returning a
#'   single vector index
#' @export
#'   
JABBA_mean_CPUE <- function(jabbafit) {
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



#' Extracts MSY reference points from a JABBA model fit
#' 
#' JABBA uses the harvest rate, H, which is catch / biomass and limited to the 
#' range 0-1. For effort control, this is converted to a fishing mortality 
#' equivalent (H=1-exp(-F)), which avoids problems with boundary conditions.
#' 
#' @inheritParams run_JABBA_MSE_list
#' @param ref_year Year to use for relative change in effort 
#' @return List of MSY reference points and a reference year for fishing 
#'   mortality for relative changes to that year
#' @export
#' 
JABBA_MSY_refpt <- function(jabba_list, 
                      ref_year = NULL) {
  
  jabba_list <- JABBA_convert2list(jabba_list)
  if (is.null(ref_year)) ref_year <- max(jabba_fit[[1]]$yr) 
  H_ref <- HMSY <- IMSY <- MSY <- m_ <- double(0)
  for (i in 1:length(jabba_list)) {
    Average_CPUE <- JABBA_mean_CPUE(jabba_list[[i]]) # A function to combine multiple indices into a single index
    H_ref <- c(H_ref, dplyr::pull(dplyr::filter(jabba_list[[i]]$kbtrj, year==ref_year), H))
    HMSY <- - log(1-c(HMSY, jabba_list[[i]]$pfunc$Hmsy))
    MSY <- c(MSY, jabba_list[[i]]$pfunc$MSY)
    IMSY = c(IMSY, Average_CPUE(JABBA_CPUE_MSY(jabba_list[[i]])))
    m_ <- c(m_, jabba_list[[i]]$pars_posterior$m)
  }
  F_ref <- - log( 1 - median(H_ref) )
  IMSY <- mean(IMSY)
  FMSY <- - log( 1 - median(HMSY) )
  MSY <- median(MSY)
  m_ <- median(m_)
  BMSY <- m_^(-1 / (m_ - 1))
  
  return(list(BMSY = BMSY,         
              IMSY = IMSY,
              FMSY = FMSY,
              ref_year = ref_year,
              F_ref = F_ref, 
              MSY = MSY))
}


#' Returns the median CPUE for each gear expected at MSY from a JABBA fit 
#'
#' @inheritParams create_JABBA_MSE
#' @return A 1D array of each CPUE at MSY
#' @export
#' 
JABBA_CPUE_MSY <- function(jabba_fit) {
  q <- jabba_fit$pars$Median[substr(rownames(jabba_fit$pars), 1, 1)=="q"]
  K <- jabba_fit$pars$Median[rownames(jabba_fit$pars)=="K"]
  m <- jabba_fit$pars$Median[rownames(jabba_fit$pars)=="m"]
  return(array(q * K * m^(-1/(m-1)), dim=c(1, length(q))))
}


#' Does minimal checks on jabba fit and ensures it is in a list 
#'
#' @inheritParams create_JABBA_MSE
#' @return jabba_fit in a list 
#' @export
#' 
JABBA_convert2list <- function(jabba_fit) {
  
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
#' @inheritParams create_JABBA_MSE
#' @inheritParams run_HCR_MSE
#' @param jabba_list a single jabba fit or a list of jabba fits
#' @param risk_rp Risk based reference levels defining targets relative to MSY 
#'   and maximum acceptable risks
#' @return A data frame with the original parameters and additional columns 
#'   containing the performance indicators.
#' @export
#' 
run_JABBA_MSE_list <- function(
    jabba_list, 
    HCR_df, 
    control_type = "Effort",
    proj_length=50, 
    nsim = 1000,
    ref_pt = standard_risk_ref_pt()) {

  jabba_list <- JABBA_convert2list(jabba_list)  
  HCR_res <- list()

  # Each fit is run in the MSE HCR simulations

  for (i in 1:length(jabba_list)) {
    HCR <- create_JABBA_MSE(jabba_list[[i]], 
                          proj_length = proj_length, 
                          nsim = nsim,
                          ref_pt = ref_pt)
    df <- run_HCR_MSE(HCR_df, HCR)
    df$fit <- names(jabba_list)[i]
    HCR_res[[i]] <- df
  }

  HCR_res <- combine_HCR_results(HCR_res, HCR_df, HCR)
  return(HCR_res)
}

