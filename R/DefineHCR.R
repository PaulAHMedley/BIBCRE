
#' Generates combinations of index, controls, change limits and the index moving 
#' average parameter as list columns in a tibble, each row being a unique HCR. 
#' Only works with single controls.
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
#' @param  ctrl_pF Effectiveness of the control (what proportion the harvest 
#'   activity is adjusted)
#' @param  ma  Moving average parameter for index time series
#' @return A tibble of all combinations of index, control, control change limit,
#'   and moving average parameters for the HCR. 
#' @export
#' 
define_HCR_test_range <- function(IMSY, fMSY, 
                                  control_type,
                                  rel_index_range = c(0.5, 1.0),
                                  rel_control_range = c(0.1, 1.0), 
                                  NInflex = 2L, 
                                  NBreaks = 5L, 
                                  change_limit = NA, 
                                  ma = 0.5,
                                  ctrl_pF = 1) {
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
    dplyr::mutate(ID = dplyr::row_number(),
                  control_type = control_type,
                  ctrl_pF = 1) |>
    dplyr::select(ID, dplyr::everything())
  
  return(df)
}


#' Defines standard risk-based reference points
#'
#' This maintains the reference points definition in one place and can be
#' edited. The reference points are in a list and consist of:
#'   - type = Basis for biomass target. e.g. 'MSY' or 'SPR40'
#'   - B_tar = Biomass at the target as a proportion of the unexploited biomass
#'   - B_tar_range = accepted target range around the target biomass point
#'   - B_lim = limit reference point as a proportion of the unexploited biomass
#'   - mostly, likely, highly likely, high certainty = probabilistic definitions
#'     of each word not currently used
#'   - max_risk = maximum acceptable risk for the candidate HCR selection.
#'     Primarily used to exlcude HCR with this risk or greater of falling below the limit
#'     reference point B_lim
#'
#' @param B_tar Target reference point relative to type
#' @return List of reference points
#' @export
#'
standard_risk_ref_pt <- function(B_tar = 1) {
  LRP <- 0.5
  return(ref_pt <- list(
    TRP = 1,
    LRP = LRP,
    B_tar = B_tar,
    B_tar_range = c(0.9, 1.2),
    B_lim = LRP*B_tar,
    mostly = 0.5,
    likely = 0.7,
    highly_likely = 0.8,
    high_certainty = 0.95,
    max_risk = 0.1
  ))
}




#### JABBAStan ####

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
define_HCR_test_range1 <- function(IMSY, fMSY,
                                  rel_index_range = c(0.5, 1.0),
                                  rel_control_range = c(0.1, 1.0),
                                  NInflex = 2,
                                  NBreaks = 5,
                                  control_type = "Effort",
                                  change_limit = NA,
                                  ma = 0.5,
                                  ctrl_pF = 1) {
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
           dplyr::mutate(ID = dplyr::row_number(),
                         control_type=control_type, ctrl_pF = ctrl_pF) |>
           dplyr::select(ID, dplyr::everything()))
}



#' Generates infection points for index/control pairs of values for a list of 
#' controls being applied that make up HCR
#' 
#' @param vInd A vector of index values to apply
#' @param vCon A vector of control values to apply
#' @param rel_index_range A list of index ranges relative to vInd
#' @param rel_control_range A list of control ranges relative to vCon
#' @param NBreaks The number of inflection points chosen for each line
#' @export
#' 
define_inflection_points <- function(vInd, vCon,
                                     rel_index_range = list(c(0.5, 1.5)),
                                     rel_control_range = list(c(0.5, 1.5)),
                                     NBreaks = 5) {
  if (any(length(vCon) != c(length(vInd), length(rel_index_range), length(rel_control_range))))
    stop("Number of controls must be apply to the length of all vectors / lists.")
  Index <- list()
  Controls <- list()
  for (i in seq_along(vCon)) {
    Index[[i]] <- seq(rel_index_range[[i]][1], rel_index_range[[i]][2], length.out=NBreaks)*vInd[i]
    Controls[[i]] <- seq(rel_control_range[[i]][1], rel_control_range[[i]][2], length.out=NBreaks)*vCon[i]
  }
  return(list(Index=Index, Controls=Controls))
}


#' Generates valid combinations of values corresponding to index or controls
#' 
#' The function returns a list of values with inflection points forming 
#' straight lines, with all valid different combinations of the values vector. 
#' 
#' @param values A vector of alternative controls or indices relative to some 
#'   value, so 1.0 would be that value
#' @param NInfex The number of inflection points chosen for each line
#' @export
#' 
generate_HCR_range <- function(values = c(0.5, 1.0, 1.5),
                               NInflex = 2) {
  NBreaks <- length(values)
  trValues <- list()
  Ir <- rep(1L, NInflex)
  
  while (TRUE) {
    for (j in Ir[NInflex - 1L]:NBreaks) {
      Ir[NInflex] <- j
      trValues <- append(trValues, list(values[Ir]))
    }
    lvl <- NInflex - 1L
    while (Ir[lvl] == NBreaks) {
      if (lvl > 1)
        Ir[lvl] <- Ir[lvl - 1L]
      lvl <- lvl - 1L
      if (lvl == 0)
        break
    }
    if (lvl == 0)
      break
    Ir[lvl] <- Ir[lvl] + 1L
  }
  return(trValues)
}



#' Plot harvest control rules with the stock status index on the x-axis and
#' control on the y-axis. 
#'
#' A ggplot is returned with the HCR plotted as a line between
#' the index of stock status and control.
#'
#' @inheritParams run_HCR_MSE
#' @inheritParams graph_sim_Btar_Ftar
#' @param HCR_ID Logical - whether to include the ID as a factor in the plot.
#'   Will only work for 12 or fewer HCR.
#' @return A ggplot object plotting the HCR's
#' @export
#'
graph_linear_HCR <- function(HCR_df,
                             ctrl_index = 1L,
                             HCR_sim = NULL,
                             HCR_ID = TRUE) {
  listify <- function(x) if (is.list(x)) x else list(x)
  HCRorder <- dplyr::select(dplyr::mutate(HCR_df, ID = ID, order=factor(dplyr::row_number())), ID, order)
  
  line_df <- HCR_df |>
    dplyr::ungroup() |>
    dplyr::mutate(trIndex = purrr::map(trIndex, ~ purrr::pluck(listify(.x), ctrl_index)),
                  trControl = purrr::map(trControl, ~ purrr::pluck(listify(.x), ctrl_index))) |>
    dplyr::select(ID, trIndex, trControl) |>
    tidyr::unnest(c(trIndex, trControl)) 
  
  maxIndex <- 1.2 * max(line_df$trIndex)
  
  lo_df <- line_df |>
    dplyr::group_by(ID) |>
    dplyr::summarise(trIndex = 0, trControl = min(trControl)) |>
    dplyr::ungroup()
  hi_df <- line_df |>
    dplyr::group_by(ID) |>
    dplyr::summarise(trIndex = maxIndex, trControl = max(trControl)) |>
    dplyr::ungroup()
  
  line_df <- dplyr::bind_rows(line_df, lo_df, hi_df) |>
    dplyr::left_join(HCRorder, by="ID") |>
    dplyr::arrange(order, trIndex)
  
  gp <- ggplot2::ggplot(line_df, 
                        ggplot2::aes(x = trIndex, y = trControl, group = order)) +
    ggplot2::geom_line() +
    ggplot2::labs(y = "Control", x = "HCR Index") +
    ggplot2::coord_cartesian( y = c(0, NA))
  
  if (HCR_ID & nrow(HCR_df) <= 12) {
    id_labels <- with(HCRorder, setNames(as.character(ID), as.character(order)))
    gp <- gp + ggplot2::facet_wrap(ggplot2::vars(order), 
                                   labeller=ggplot2::labeller(order = id_labels))
  }
  if (!is.null(HCR_sim)) {
    pv_df <- tibble::tibble(pvIndex = pvIndex, pvControl = pvControl)
    tr_df <- tibble::tibble(trIndex = trIndex, trControl = trControl)
    gp <- gp +
      ggplot2::geom_point(pv_df, 
                          ggplot2::aes(x = pvIndex, y = pvControl)) +
      ggplot2::geom_line(tr_df, 
                         ggplot2::aes(x = pvIndex, y = pvControl), color = "blue")
  }
  return(gp)
}

