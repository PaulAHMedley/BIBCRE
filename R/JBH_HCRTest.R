
# #### HCR tests: Graphs and Tables ####


#' Plot harvest control rules with the stock status index on the x-axis and 
#' control on the y-axis.
#'   
#' A ggplot is returned with the HCR plotted as a line between 
#' the index of stock status and control.
#'   
#' @inheritParams HCR_performance
#' @inheritParams run_HCR_MSE
#' @param HCR_ID Logical - whether to include the HCR ID as a factor in the 
#'   plot. Will only work for 12 or fewer HCR.
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
  
  if (HCR_ID & nrow(HCR_df) <= 12) gp <- gp + ggplot2::facet_wrap(vars(ID))
  
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
