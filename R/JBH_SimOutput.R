
# #### Individual Sims: Graphs and Tables ####


#' Plot the stock status and harvest rate for a single HCR test simulation.
#'   
#' A ggplot is returned with the biomass as a proportion of BMSY and H as a 
#' proportion of HMSY over time.
#'   
#' @inheritParams HCR_performance
#' @param type Either "both", "bmsy" or "hmsy" to indicate which plots are 
#'  provided 
#' @return A probability density plot of biomass dynamics projection under a HCR
#' @export
#' 
graph_sim_BMSY_FMSY <- function(HCR_sim, type = "both") {
  rsl_B <- rsl_H <- NULL
  
  if (type %in% c("both", "hmsy")) {
    catch_mat <- with(HCR_sim, matrix(Bio$Catch, nrow=HCR$nsim, ncol=HCR$TN))
    #invHMSY <- with(HCR_sim$Par, m / r)     
    obHMSY <- with(HCR_sim, apply(catch_mat, 2, function(x) x*invHMSY/Par$K))
    pjHMSY  <- with(HCR_sim, apply(C, 2, function(x) x*invHMSY))
    HMSY  <- as.vector(cbind(obHMSY, pjHMSY))
    
    rsl_H <- with(HCR_sim, tibble::tibble(year = rep(Bio$year[1]:(Bio$year[1] + HCR$PTN-1), each=HCR$nsim),
                                          Var = "HMSY", Val = HMSY))
  }
  
  if (type %in% c("both", "bmsy")) {
    #invBMSY  <- with(HCR_sim$Par, m^(1/(m-1)))
    pBMSY <- with(HCR_sim, apply(pB, 2, function(x) x*invBMSY))
    rsl_B <- with(HCR_sim, tibble::tibble(year = rep(Bio$year[1]:(Bio$year[1] + HCR$PTN), each=HCR$nsim),
                                        Var = "BMSY",
                                        Val = as.vector(pBMSY))) 
  }
  
  df <- switch(type,
         both = rbind(rsl_B, rsl_H),
         bmsy = rsl_B,
         hmsy = rsl_H
  )
  if (is.null(df)) stop("Error: type must be a string: 'both', 'bmsy' or 'hmsy'.")
  
  df <- dplyr::mutate(df, year=year+0.5)
  
  ggp <- dplyr::filter(df, Val > 0, Val < 2.4, !is.na(Val)) |>
    ggplot2::ggplot(ggplot2::aes(x=year, y=Val)) +
    ggplot2::geom_bin2d(bins=HCR_sim$HCR$PTN-1, na.rm = TRUE) +
    #geom_line(aes(colour=Type)) +
    ggplot2::scale_y_continuous(limits = c(0, 2.4))+
    ggplot2::geom_vline(xintercept=(HCR_sim$Bio$year[1] + HCR_sim$HCR$TN-1))+
    ggplot2::geom_hline(yintercept=1.0)+
    ggplot2::geom_hline(yintercept=0.5)+
    ggplot2::scale_fill_gradient(low = "#96B1F7", high = "#030B43") 

  if (type == "both") 
    return(ggp +
      ggplot2::facet_wrap(dplyr::vars(Var), nrow=2) +
      ggplot2::labs(x="Year", y="Biomass / harvest rate relative to MSY"))
  
  if (type == "hmsy")
    return(ggp +
             ggplot2::labs(x="Year", y="Harvest rate relative to MSY"))
  if (type == "bmsy")
    return(ggp +
             ggplot2::labs(x="Year", y="Biomass relative to MSY"))
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
  
  ggplot2::ggplot(rsl_C, ggplot2::aes(x=year, y=Catch)) +
    ggplot2::geom_bin2d(bins=HCR_sim$HCR$PTN-2, na.rm = TRUE) +
    ggplot2::scale_fill_gradient(low = "#96B1F7", high = "#030B43") +
    ggplot2::labs(x="Year", y="Catch(t)")
} 


#' Table the proportion of each state-response for HCR decisions
#'
#' The decisions compare stock status versus response, where an initial response
#' is indicated by the HCR index being below the highest trigger point in the
#' HCR or below the lowest trigger in response to being below Blim.
#'
#' @inheritParams HCR_performance
#' @return flextable object containing the state-response proportions for the
#'   HCR
#' @export
#' 
table_sim_decision <- function(HCR_sim) {
  Lim <- HCR_sim$ref_pt$B_lim
  Trig <- HCR_sim$ref_pt$B_tar_range[1]
  lvlBt <- c("B<Limit",   #Limit
             "B<Trigger", #Trigger
             "B~Target")  #Target
  lvlIt <- c("It<Limit",  #Limit
             "It<Trigger", #Trigger
             "It~Target")  #Target
  
  #BB0 <- as.vector(HCR_sim$pB[ , (HCR_sim$HCR$TN+1):(HCR_sim$HCR$PTN)])
  pBMSY <- with(HCR_sim, as.vector(apply(pB[ , (HCR$TN+1):(HCR$PTN)], 
                                         2, function(x) x*invBMSY)))
  
  Index <- as.vector(HCR_sim$pjIndex[ , -1])
  HCR_r <- with(HCR_sim$HCR, tibble::tibble(State = lvlBt[3 - ((pBMSY < Lim) + (pBMSY < Trig))],
                          Response=lvlIt[3 - ((Index < trIndex[1]) + 
                                                (Index < trIndex[length(trIndex)]))])) |>
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
  std_brk <- with(HCR_sim$ref_pt, c(0.0, B_lim, B_tar_range, 100.0))
  c_std_brk <- paste0(format(std_brk[-length(std_brk)],nsmall=1), "-",
                      format(std_brk[-1], nsmall=1))
  c_std_brk[length(c_std_brk)] <- paste0(">", format(std_brk[length(std_brk)-1L],nsmall=1))
  
  # Projection statistics
#  BMSY <- with(HCR_sim, 2 * as.vector(pB[ ,(HCR$TN+1):HCR$PTN]))
  HMSY  <- with(HCR_sim, as.vector(C) * invHMSY)
  pBMSY <- with(HCR_sim, apply(pB[ , (HCR$TN+1):(HCR$PTN)], 
                               2, function(x) x*invBMSY))
  
  cnt <- tibble::tibble(
    `B/BMSY` = c_std_brk,
    `B%` = hist(pBMSY, breaks=std_brk, plot=FALSE)$counts,
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
