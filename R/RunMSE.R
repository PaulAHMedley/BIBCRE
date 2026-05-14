

#' Run an HCR MSE based on a table of parameters and HCR MSE function
#'
#' The function runs the HCR MSE based on the set of HCR parameters in the
#' dataframe and the provided MSE function. The performance indicators are
#' added as additional columns to the dataframe.
#' 
#' The data frame HCR_df must have the fields: ID, trIndex, trControl,
#' control_type, change_limit, ma and ctrl_pF.
#'
#' @param HCR_df A data frame produced by eg [define_HCR_test_range] with
#'   HCR parameters
#' @param HCR_MSE  HCR function produced by [create_fishblicc_MSE], 
#'   [create_JABBA_MSE] or [create_ptStan_MSE]
#' @return A data frame with the original parameters and additional columns
#'   containing the performance indicators.
#' @export
#'
#'
run_HCR_MSE <- function(HCR_df, HCR_MSE) {
  HCR_df <- HCR_df |>
    dplyr::select(ID, trIndex, trControl, control_type, change_limit, ma, ctrl_pF) |>
    dplyr::mutate(Res = purrr::pmap(
      list(trIndex, trControl, control_type, change_limit, ma, ctrl_pF),
      \(trIndex, trControl, change_limit, ma, control_type, ctrl_pF)
        HCR_performance(HCR_MSE(trIndex, trControl, change_limit, ma, control_type, ctrl_pF)),
      .progress = "HCR Sim"
    )) |>
    tidyr::unnest(cols=Res)
  HCR_df <- evaluate_HCR(HCR_df, HCR_MSE)
  return(HCR_df)
}


#' Parallel version of [run_HCR_MSE]
#'
#' The function runs the HCR MSE based on the set of HCR parameters in the
#' dataframe and the provided MSE function. The performance indicators are
#' added as additional columns to the dataframe. This should be faster than the 
#' serial version when running several HCR. The function uses the package 
#' `furrr`.
#'
#' @inheritParams run_HCR_MSE
#' @return A data frame with the original parameters and additional columns
#'   containing the performance indicators.
#' @export
#'
run_HCR_MSE_para <- function(HCR_df, HCR_MSE) {
  old_plan <- future::plan()          # save current plan
  on.exit(future::plan(old_plan), add = TRUE)  # restore plan safely
  cores <- future::availableCores()
  future::plan(multisession, workers = cores)  # works on Windows/macOS/Linux
  RNGkind("L'Ecuyer-CMRG")
  rseed <- get("rseed", envir=environment(HCR_MSE))
  set.seed(rseed) # same seed for all HCR
  factory_stream <- parallel::nextRNGStream(.Random.seed)
  sim_stream     <- parallel::nextRNGStream(factory_stream)
  .Random.seed <- factory_stream
  
  HCR_df <- HCR_df |>
    dplyr::select(ID:ma) |>
    dplyr::mutate(Res = furrr::future_pmap(
      list(trIndex, trControl, change_limit, ma),
      \(trIndex, trControl, control_type, change_limit, ma, ctrl_pF)
      HCR_performance(HCR_MSE(trIndex, trControl, control_type, change_limit, ma, ctrl_pF)),
      .options = furrr_options(seed = TRUE),
      .progress = TRUE
    )) |>
    tidyr::unnest(cols=Res)
  HCR_df <- evaluate_HCR(HCR_df, HCR_MSE)
  return(HCR_df)
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
    dplyr::mutate(Evaluation = ifelse((at_Btar >= ref_pt$mostly) & 
                                        (lt_Blim <= ref_pt$max_risk), 
                                      "Candidate", "Rejected")) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      Rank = max(Catch_Rank, State_Rank)) |>
    dplyr::ungroup()
  return(HCR_df)
}


