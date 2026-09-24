
# #### HCR tests: Graphs and Tables ####


#' Calculate the performance measures of a harvest control rule and return it
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
  else {
    PjN <- with(HCR_sim$HCR, PN + 1L) }
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
    CPUE <- with(HCR_sim, apply(CW/mF, MARGIN=c(1,2), FUN=sum))
    CPUE <- sweep(CPUE, MARGIN=1, STATS=CPUE[,1], FUN="/")
    CPUE_tar <- NA
  } else {  
    # if (retro) {
    #   Bt <- with(HCR_sim, as.vector(pB[ , 1L:(HCR$TN)]))
    #   Index <- as.vector(HCR_sim$pvIndex)
    #   Catch <- HCR_sim$dat$TCA_ca
    # } else {
      Bt <- with(HCR_sim, as.vector(pB[ , (HCR$TN+1L):(HCR$PTN+1L)]))
      Index <- as.vector(HCR_sim$pjIndex[ , -1L])
      Catch <- with(HCR_sim, C[ , (HCR$TN+1L):HCR$PTN]) #with(HCR_sim, sweep(C[ , (HCR$TN+1):HCR$PTN], MARGIN=1, STATS=Par$Binf, FUN="*"))
      CPUE <- as.vector(HCR_sim$CPUE)
      #    }
  }
  
  Catch_Avg <- mean(Catch)
  Catch_Rng <- mean(abs(Catch-Catch_Avg))
  Catch_pcile <- quantile(Catch, probs=MaxRisk)
  CPUE_Avg <- mean(CPUE)
  lt_Blim <- sum(Bt < B_lim)/length(Bt)
  at_Btar <- sum((Bt >= B_tar_lower) & (Bt <= B_tar_higher)) / length(Bt)
  gt_Btar <- sum(Bt > B_tar_higher)/length(Bt)
  # Management response when not necessary
  lo_trigger <- min(unlist(HCR_sim$HCR$trIndex))
  hi_trigger <- max(unlist(HCR_sim$HCR$trIndex))
  Err_Type1 <- 0.5*mean(((Bt > B_lim) & (Index < lo_trigger)) +
                          ((Bt > B_tar_lower) & (Index < hi_trigger)))
  # No management response when it is necessary
  Err_Type2 <- 0.5*mean(((Bt < B_lim) & (Index > lo_trigger)) +
                          ((Bt < B_lim) & (Index > hi_trigger)))
  
  return(tibble::tibble(
    Catch_Avg = Catch_Avg,
    Catch_Rng = Catch_Rng,
    Catch_pcile = Catch_pcile,
    CPUE_Avg = CPUE_Avg,
    lt_Blim = lt_Blim,
    at_Btar = at_Btar,
    gt_Btar = gt_Btar,
    State = at_Btar - lt_Blim,
    Err_Type1 = Err_Type1,
    Err_Type2 = Err_Type2
  ))
}


#' Calculate retrospective performance from fitted model
#'
#' Calculate the performance measures for observed past performance for the 
#' fitted model from the MSE function. See [HCR_performance] for interpretation
#' of the performance measures.
#'
#' @inheritParams run_HCR_MSE
#' @return reference point list updated with past performance values
#' @export
#'
past_performance <- function(HCR_MSE) {
  stock_assessment <- get("stock_assessment", envir=environment(HCR_MSE))
  if (stock_assessment == "fishblicc")
    return(NULL)  # Not relevant
  pB <- get("pB", envir=environment(HCR_MSE))
  dat <- get("dat", envir=environment(HCR_MSE))
  ref_pt <- get("ref_pt", envir=environment(HCR_MSE))

  PjN <- with(dat, TN)
  B_tar <- rep(ref_pt$B_tar, PjN)
  B_lim <- rep(ref_pt$B_lim, PjN) # 0.5
  MaxRisk <- ref_pt$max_risk
  
  B_tar_lower <- rep(ref_pt$B_tar_range[1]*ref_pt$B_tar, PjN) 
  B_tar_higher <- rep(ref_pt$B_tar_range[2]*ref_pt$B_tar, PjN) 
  Bt <- as.vector(pB[ , seq_len(dat$TN)])

  Catch_Avg <- mean(dat$TCA_ca)
  Catch_Rng <- mean(abs(dat$TCA_ca-Catch_Avg))
  Catch_pcile <- quantile(dat$TCA_ca, probs=MaxRisk)
  CPUE_Avg <- mean(with(dat, TCE_ca/TCE_ef), na.rm=TRUE)
  lt_Blim <- sum(Bt < B_lim)/length(Bt)
  at_Btar <- sum((Bt >= B_tar_lower) & (Bt <= B_tar_higher)) / length(Bt)
  gt_Btar <- sum(Bt > B_tar_higher)/length(Bt)
  # No HCR
  Err_Type1 <- ref_pt$max_risk
  Err_Type2 <- ref_pt$max_risk
  
  return(tibble::tibble(
    Catch_Avg = Catch_Avg,
    Catch_Rng = Catch_Rng,
    Catch_pcile = Catch_pcile,
    CPUE_Avg = CPUE_Avg,
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
#'
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


#' Plot the harvest control rule performance indicators
#'
#' An XY plot is returned with the requested performance indicators. 
#' 
#' If requested (eval==TRUE) candidate HCR are plotted as black points and 
#' rejected HCR as red points. 
#' 
#' A linear fit is shown for plots covering catches and CPUE. The line goes 
#' through the origin indicating the exchange between increasing 
#' catch with increasing catch range. HCR above or below this line may indicate 
#' a better "exchange rate" between measures where they conflict. 
#' 
#' If requested (ref_pt==TRUE) and available, vertical 
#' and horizontal dotted lines represent the reference points for each 
#' indicator. Reference points are either the performance actually achieved 
#' historically (proj_length=0 for the production models only) or based on the 
#' average performance across HCR being examined (all catch based performance 
#' indicators), the CPUE at MSY for the CPUE_Avg indicator, or the risk-based 
#' reference points set in the projection function for the stock state 
#' indicators and error type indicators.
#' 
#' The variables that can be plotted are:
#' Catch_Avg : Average catch
#' Catch_Rng : Catch range
#' Catch_pcile : Lower catch 10 percentile 
#' CPUE_Avg : Average CPUE
#' lt_Blim : Proportion of the simulations below the limit reference point 
#' at_Btar : Proportion of the simulations in the target region
#' State : at_Btar - lt_Blim
#' gt_Btar : Greater than the target region
#' Err_Type1 : Proportion of the simulations with a false positive
#' Err_Type2 : Proportion of the simulations with a false negative
#'
#' The reference points and lower catch 10 percentile can be adjusted in the 
#' model's reference point list. 
#' 
#' For the proportion of Type 1 (false positive), the decision rule results in 
#' management intervention reducing harvest when the stock in reality is at or 
#' above its target level. For Type 2 (false negative), the decision rule does 
#' not harvest when the stock is in reality below the limit reference point, 
#' which in this case is considered worse as it is less precautionary. Higher 
#' error rates indicate poorer HCR performance and might be improved by 
#' adjusting the index trigger points.
#' 
#' @inheritParams graph_linear_HCR
#' @param PIX Performance indicator to go on the X axis
#' @param PIY Performance indicator to go on the Y axis
#' @param ref_pt Whether to plot the performance reference points as 
#'   vertical/horizontal dotted lines
#' @param eval Whether to plot the evaluation for candidate vs rejected HCR
#' @param avg_line whether to include a line regression through the performance
#'   indicators
#' @return An XY point plot of the HCR performance indicators
#' @export
#'
graph_HCR_performance <- function(HCR_df, PIX, PIY, 
                                  ref_pt = TRUE, 
                                  eval = FALSE,
                                  avg_line = FALSE) {
  pi_text <- c("Catch_Avg", "Catch_Rng", "Catch_pcile", "CPUE_Avg", "lt_Blim", 
               "at_Btar", "State", "gt_Btar", "Err_Type1", "Err_Type2")
  lbl_text <- c("Average Catch", "Catch Range", "Catch Lower Percentile", 
                "Average CPUE", "Biomass Less than Limit", 
                "Biomass At Target", "Stock State", "Biomass Above Target", 
                "Error Type 1", "Error Type 2")
  names(lbl_text) <- pi_text
  
  PIX_enq <- rlang::enquo(PIX)
  PIY_enq <- rlang::enquo(PIY)
  # PIX_text <- rlang::quo_text(PIX_enq)
  # PIY_text <- rlang::quo_text(PIY_enq)
  PIX_text <- rlang::as_name(PIX_enq)
  PIY_text <- rlang::as_name(PIY_enq)

  if (!(PIY_text %in% pi_text & PIX_text %in% pi_text)) {
    stop(paste("Error: No such column name: Valid names are ", paste(pi_text, collapse=", ")))
  }
  
  gg <- ggplot2::ggplot(dplyr::filter(HCR_df, ID > 0), 
                        ggplot2::aes(x = !!PIX_enq, y = !!PIY_enq)) + 
    ggplot2::labs(x = lbl_text[PIX_text], y = lbl_text[PIY_text])
  
  if (eval) {
    gg <- gg +
    ggplot2::geom_point(ggplot2::aes(color = Evaluation)) + 
      ggplot2::scale_color_manual(values = c(Candidate = "black", 
                                             Rejected = "darkred"))  
  } else {
    gg <- gg + ggplot2::geom_point() 
  }

  if (avg_line) {
    # Plot line if there is an 'exchange rate'
    gg <- gg +
      ggplot2::geom_smooth(method = "lm", 
                           formula = y ~ x + 0, se = FALSE, linetype = "solid", 
                           alpha = 0.4)
  }
  if (PIX_text %in% c("Err_Type1", "Err_Type2")) {
    gg <- gg +
      ggplot2::coord_cartesian(xlim = c(0, NA))
  }
  if (PIY_text %in% c("Err_Type1", "Err_Type2")) {
    gg <- gg +
      ggplot2::coord_cartesian(ylim = c(0, NA))
  }
  
  if (ref_pt) {
    df <- dplyr::filter(HCR_df, ID == 0)
    if (nrow(df) == 1L) {
      rp_Xaxis <- dplyr::pull(df, !!PIX_enq)
      rp_Yaxis <- dplyr::pull(df, !!PIY_enq)
      gg <- gg + 
        ggplot2::geom_hline(yintercept = rp_Yaxis, linetype = "dotted") + 
        ggplot2::geom_vline(xintercept = rp_Xaxis, linetype = "dotted")
    }}
  return(gg)
}


#' Table the standard performance measures for a set of candidate HCRs
#'
#'
#' @inheritParams graph_linear_HCR
#' @return A `flextable` containing the HCR performance measures
#' @export
#'
table_HCR_performance <- function(HCR_df) {
  HCR_df |>
    dplyr::filter(ID > 0) |>
    dplyr::select(ID, change_limit, ma, Catch_Avg, Catch_Rng, Catch_pcile, 
                  CPUE_Avg, lt_Blim, at_Btar, gt_Btar, State, 
                  Err_Type1, Err_Type2, Catch_Rank, State_Rank) |>
    flextable::flextable() |>
    flextable::set_header_labels(values = c("ID", "Change Limit", "Index Smoother",
                                            "Mean Catch", "Mean Catch Range", 
                                            "Catch Lower Percentile", "Mean CPUE",
                                            "Below LRP", "In Target Range", "Above TRP", 
                                            "Good Fishery State", 
                                            "HCR False Positive", "HCR False Negative", 
                                            "Catch Rank", "State Rank")) |>
    flextable::colformat_int(j=1, big.mark="") |>
    flextable::colformat_double(j=4:6, digits=0, big.mark="") |>
    flextable::colformat_double(j=7:12, digits=3) |>
    flextable::colformat_double(j=13:14, digits=0, big.mark="") |>
    flextable::autofit() 
}


#' JABBAstan Combines multiple HCR simulation results into a single tibble
#'
#' @inheritParams run_HCR_MSE
#' @param HCR_res_list A list of results from [HCR_performance]
#' @return A tibble containing combined results from the simulations
#' @export
#'
combine_HCR_results <- function(HCR_res_list, HCR_MSE) {
  HCR_res_list <- dplyr::bind_rows(HCR_res_list) |>
    dplyr::filter(ID > 0)

  HCR_df <- dplyr::group_by(HCR_res_list, ID) |>
    dplyr::summarise(trIndex = list(dplyr::first(trIndex)), 
                     trControl = list(dplyr::first(trControl)), 
                     control_type = list(dplyr::first(control_type)), 
                     change_limit = dplyr::first(change_limit), 
                     ma = dplyr::first(ma), 
                     ctrl_pF = list(dplyr::first(ctrl_pF)), 
                     .groups = "drop")
  
  # Assume each model is equally weighted so use mean performance
  df <- HCR_res_list |>
    dplyr::group_by(ID) |>
    dplyr::summarise(Catch_Avg = mean(Catch_Avg),
                     Catch_Rng = mean(Catch_Rng),
                     Catch_pcile = min(Catch_pcile),  # this is not quite correct, but precautionary
                     CPUE_Avg = mean(CPUE_Avg),  
                     lt_Blim = mean(lt_Blim),
                     at_Btar = mean(at_Btar),
                     gt_Btar = mean(at_Btar),
                     Err_Type1 = mean(Err_Type1),
                     Err_Type2 = mean(Err_Type2)) |>
    dplyr::ungroup() |>
    dplyr::mutate(State = at_Btar - lt_Blim) |>
    dplyr::select(ID:at_Btar, State, everything())
  
  HCR_res <- HCR_df |>
    dplyr::left_join(df, by="ID") |>
    evaluate_HCR(HCR_MSE)
  return(HCR_res)
}



