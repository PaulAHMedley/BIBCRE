
# #### HCR tests: Graphs and Tables ####


#' JABBA version Calculate the performance measures of a harvest control rule and return it
#' as a row in a data frame
#'
#' The performance measures consist of the 1) average catch, 2) the mean catch
#' range around the average (mean(abs(catch-mean_catch))), 3) the HCR Type 1
#' error rate, 4) the HCR Type 2 error rate, 5) the proportion of time that the
#' stock is below 0.5 Btar (LRP) and 6) the proportion of time that the stock
#' is outside 90-110% Btar.
#'
#' The HCR error rates are scores indicating how well the HCR is responding.
#' For Type 1, When the stock is below the limit reference point and the index
#' is below the upper trigger but above the lower trigger, it scores one,
#' whereas if it is above both triggers it scores two. The mean of this value
#' is taken. Similarly for Type 2, responding by reducing harvest when the
#' stock is at target scores two when the reduction is minimised, one otherwise.
#' The higher the scores are, the worse the performance.
#'
#' @param HCR_sim Fishery HCR simulation result list from the MSE
#'   (see [create_fishblicc_MSE], [create_JABBA_MSE], [create_ptStan_MSE]).
#' @return Tibble row of performance measures
#' @export
#'
HCR_performance <- function(HCR_sim) {
  if (HCR_sim$stock_assessment == "fishblicc")
    PjN <- HCR_sim$HCR$PYN + 1L
  else
    PjN <- with(HCR_sim$HCR, PN + 1L)
  B_tar <- with(HCR_sim, rep(ref_pt$B_tar, PjN))
  B_lim <- with(HCR_sim, rep(ref_pt$B_lim, PjN)) # 0.5
  MaxRisk <- HCR_sim$ref_pt$max_risk
  Mostly <- HCR_sim$ref_pt$mostly
  
  B_tar_lower <- with(HCR_sim, rep(ref_pt$B_tar_range[1]*ref_pt$B_tar, PjN)) #c(0.90, 1.2)
  B_tar_higher <- with(HCR_sim, rep(ref_pt$B_tar_range[2]*ref_pt$B_tar, PjN)) #c(0.90, 1.2)
  
  if (HCR_sim$stock_assessment == "fishblicc") {
    Bt <- with(HCR_sim, as.vector(SSB))
    Index <- as.vector(HCR_sim$pjIndex)
    Catch <- HCR_sim$CW
  } else {  
    Bt <- with(HCR_sim, as.vector(pB[ , (HCR$TN+1L):(HCR$PTN+1L)]))
    Index <- as.vector(HCR_sim$pjIndex[ , -1L])
    Catch <- with(HCR_sim, C[ , (HCR$TN+1L):HCR$PTN]) #with(HCR_sim, sweep(C[ , (HCR$TN+1):HCR$PTN], MARGIN=1, STATS=Par$Binf, FUN="*"))
  }
  
  Catch_Avg <- mean(Catch)
  Catch_Rng <- mean(abs(Catch-Catch_Avg))
  Catch_pcile <- quantile(Catch, probs=MaxRisk)
  lt_Blim <- sum(Bt < B_lim)/length(Bt)
  at_Btar <- sum((Bt >= B_tar_lower) & (Bt <= B_tar_higher)) / length(Bt)
  gt_Btar <- sum(Bt > B_tar)/length(Bt)
  # Management response when not necessary
  lo_trigger <- min(unlist(HCR_sim$HCR$trIndex))
  hi_trigger <- max(unlist(HCR_sim$HCR$trIndex))
  Err_Type1 <- 0.5*mean(((Bt > B_lim) & (Index < lo_trigger)) +
                          ((Bt > B_tar_lower) & (Index < hi_trigger)))
  # No management response when it is necessary
  Err_Type2 <- 0.5*mean(((Bt < B_lim) & (Index > lo_trigger)) +
                          ((Bt < B_lim) & (Index > hi_trigger)))
  
  return(tibble(
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


#' Ranks HCR based on performance indicators and evaluates whether the HCR is a 
#' candidate.
#'   
#' The resulting ranks and whether the HCR is a candidate or not is added as
#' columns to the data frame.
#' 
#' Ranks are based on catch (10th percentile) and stock state (proportion of 
#' the simulation spent around the target level - that spent below Blim).
#' @inheritParams run_HCR_MSE
#' @return HCR data frame with added evaluation and ranks
#' @export
#' 
rank_HCR <- function(HCR_df, HCR_MSE) {
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
  
  ggplot2::ggplot(HCR_df, ggplot2::aes(x=Catch_Avg, y=Catch_Rng)) +
    ggplot2::geom_point(ggplot2::aes(color=Evaluation)) +
    ggplot2::scale_color_manual(values = c("Candidate" = "black", "Rejected" = "red")) +
    ggplot2::geom_hline(yintercept=Var_Catch, linetype="dotted") +
    ggplot2::geom_vline(xintercept=Mean_Catch, linetype="dotted") +
    ggplot2::geom_smooth(method="lm", formula = y~x+0, se=FALSE, linetype="solid", alpha=0.5) + 
    ggplot2::labs(x="Average Catch", y = "Catch Range")
  
}


#' Plot the HCR stock status relative to the Btar (TRP) and Blim (LRP)
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
  ggplot2::ggplot(HCR_df, ggplot2::aes(x=lt_Blim, y=at_Btar)) +
    ggplot2::geom_point(ggplot2::aes(color=Evaluation)) +
    ggplot2::scale_color_manual(values = c("Candidate" = "black", "Rejected" = "red")) +
    ggplot2::geom_hline(yintercept=ref_pt$mostly, linetype="dotted") +
    ggplot2::geom_vline(xintercept=ref_pt$max_risk, linetype="dotted") +
    ggplot2::labs(x = "Proportion Below B_lim", y = "Proportion in target range")
  
}


#' Plot the HCR performance with respect to stock status and catch
#'   
#' A ggplot is returned with the the proportion of the simulation time across 
#' all runs for an HCR stock status (x-axis) and lower 10%_ile catch (y-axis). 
#' HCR change limits indicated by colour. 
#'   
#' @inheritParams graph_linear_HCR
#' @return A ggplot object plotting the HCR performance catches
#' @export
#' 
graph_HCR_status_catches <- function(HCR_df) {
  
  HCR_df |>
    ggplot2::ggplot(ggplot2::aes(x=State, y=Catch_pcile, colour=factor(change_limit))) +
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
  
  ggplot2::ggplot(HCR_df, ggplot2::aes(x=Err_Type1, y=Err_Type2)) +
    ggplot2::geom_point(ggplot2::aes(color=Evaluation)) +
    ggplot2::scale_color_manual(values = c("Candidate" = "black", "Rejected" = "red")) +
    ggplot2::geom_hline(yintercept=MaxRisk, linetype="dotted") +
    ggplot2::geom_vline(xintercept=MaxRisk, linetype="dotted") +
    ggplot2::labs(x = "Type 1 (False Positive)", y = "Type 2 (False Negative)")
}



#' JABBAstan Table the standard performance measures for a HCR simulation
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



#' Table the standard performance measures for a set of candidate HCRs
#' 
#' 
#'
#' @inheritParams graph_linear_HCR
#' @return A `flextable` containing the HCR performance measures
#' @export
#' 

table_HCR_performance <- function(HCR_df) {
  HCR_df |>
    dplyr::select(ID, change_limit, ma, Catch_Avg, Catch_Rng, Catch_pcile, 
                  lt_Blim, at_Btar, gt_Btar, State, 
                  Err_Type1, Err_Type2, Catch_Rank, State_Rank) |>
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


#' JABBAstan Combines multiple HCR simulation results into a single tibble
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
                     at_Btar = mean(at_Btar),
                     gt_Btar = mean(at_Btar),
                     Err_Type1 = mean(Err_Type1),
                     Err_Type2 = mean(Err_Type2)) |>
    dplyr::ungroup() |>
    dplyr::mutate(State = at_Btar - lt_Blim) |>
    dplyr::select(ID:at_Btar, State, everything())
  
  HCR_res <- HCR_df |>
    dplyr::select(ID, trIndex, trControl, control_type, change_limit, ma, ctrl_pF) |>
    dplyr::left_join(df, by="ID") |>
    evaluate_HCR(HCR_MSE)
  return(HCR_res)
}



