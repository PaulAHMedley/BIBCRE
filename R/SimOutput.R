
# Individual Sims: Graphs and Tables ####


#' Plot the stock status and harvest rate for a single HCR test simulation.
#'
#' A ggplot is returned with the biomass as a proportion of BMSY and H as a
#' proportion of HMSY over time.
#'
#' @inheritParams HCR_performance
#' @return A probability density plot of biomass dynamics projection under a HCR
#' @export
#'
graph_sim_Btar_Ftar <- function(HCR_sim, type = "both") {
  rsl_B <- rsl_F <- NULL

  if (HCR_sim$stock_assessment == "fishblicc") {
    if (type %in% c("both", "f")) {
      mF <- HCR_sim$mF
      nsim <- HCR_sim$HCR$nsim
      PYN <- HCR_sim$HCR$PYN
      F_tar <- HCR_sim$ref_pt$F_tar
      NG <- ncol(F_tar)
      pFtar <- matrix(0, nrow=nsim, ncol=PYN+1L)
      for (i in seq_len(NG))
        pFtar <- pFtar + sweep(mF[ ,, i], 1, STATS=NG * F_tar[,i], FUN="/")
      
      rsl_F <- with(HCR_sim$HCR,
                    tibble::tibble(year = rep(StartYear - 1 + seq_len(PYN+1L), each=nsim),
                                   Var = "F",
                                   Val = as.vector(pFtar)))
    }
    
    if (type %in% c("both", "b")) {
      pBtar  <- with(HCR_sim, sweep(SSB, 1, STATS=ref_pt$B_tar, FUN="/"))
      rsl_B <- with(HCR_sim$HCR,
                    tibble::tibble(year = rep(StartYear - 1 + seq_len(PYN+1L), each=nsim),
                                   Var = "SSB",
                                   Val = as.vector(pBtar)))
    }
  } else {
    
    ## Fishing Mortality
    if (type %in% c("both", "f")) {
      pFtar <- with(HCR_sim, sweep(Ft+0.0001, MARGIN=1, STATS=ref_pt$F_tar, FUN="/"))
      
      years <- with(HCR_sim, HCR$StartYear - 1 + 1:HCR$PTN)
      rsl_F <- with(HCR_sim, tibble::tibble(year = rep(years, each=HCR$nsim),
                                            Var = "F", 
                                            Val = as.vector(pFtar)))
    }
    
    ## Biomass
    if (type %in% c("both", "b")) {
      years <- with(HCR_sim, HCR$StartYear - 1 + 1:(HCR$PTN+1L))
      pBtar <- with(HCR_sim, sweep(pB, MARGIN=1, STATS=ref_pt$B_tar, FUN="/"))
      rsl_B <- with(HCR_sim, tibble::tibble(year = rep(years, each=HCR$nsim),
                                            Var = "B",
                                            Val = as.vector(pBtar)))
    }  
  }
  
  df <- switch(type,
               both = rbind(rsl_B, rsl_F),
               b = rsl_B,
               f = rsl_F
  )
  if (is.null(df))
    stop("Error: type must be a string: 'both', 'b' or 'f'.")
  
  rp_basis <- paste0(HCR_sim$ref_pt$rp_type, ": ", format(HCR_sim$ref_pt$TRP))
  maxvalue <- pmax(2.4, quantile(df$Val, 0.99, names = FALSE))
  binwidth_x <- 1   # one time step per bin
  binwidth_y <- diff(quantile(df$Val, c(0.01, 0.99))) / 50  # example
  
  ggp <- df |>
    dplyr::filter(!is.na(Val), Val > 0, Val < maxvalue) |>
    ggplot2::ggplot(ggplot2::aes(x=year, y=Val)) +
    ggplot2::geom_bin2d(binwidth = c(binwidth_x, binwidth_y), na.rm = TRUE) +
    ggplot2::scale_y_continuous(limits = c(0, maxvalue)) +
    ggplot2::geom_vline(xintercept=(HCR_sim$HCR$StartYear + HCR_sim$HCR$TN))+
    ggplot2::geom_hline(yintercept=1.0)+
    ggplot2::geom_hline(yintercept=0.5)+
    ggplot2::scale_fill_gradient(low = "#96B1F7", high = "#030B43")

  if (type == "both")
    return(ggp +
             ggplot2::facet_wrap(ggplot2::vars(Var), nrow=2, scales="free_y") +
             ggplot2::labs(x="Year",
                           y=paste0("SSB / F relative to target ", rp_basis)))
  
  if (type == "f")
    return(ggp +
             ggplot2::labs(x="Year",
                           y=paste0("Fishing mortality relative to target ",
                                    rp_basis)))
  if (type == "b")
    return(ggp +
             ggplot2::labs(x="Year",
                           y=paste0("Spawning biomass relative to target ",
                                    rp_basis)))
}


#' Plot the catches for a single HCR test simulation.
#'
#' @inheritParams HCR_performance
#' @return graph showing binned probability for observed and simulated catch
#'   time series
#' @export
#'
graph_sim_catch <- function(HCR_sim) {
  #catch_df <- with(HCR_sim, sweep(C, MARGIN=1, STATS=Par$Binf, FUN="*"))
  if (HCR_sim$stock_assessment == "fishblicc") {
    rsl_C <- with(HCR_sim,
                  tibble::tibble(year = rep(HCR$StartYear - 1 + seq_len(HCR$PYN), each=HCR$nsim),
                                 Catch = as.vector(CW)))
  } else {
    rsl_C <- with(HCR_sim,
                  tibble::tibble(year = rep(HCR$StartYear - 1 + seq_len(HCR$PTN), each=HCR$nsim),
                                 Catch = as.vector(HCR_sim$C)))
  }
  binwidth_x <- 1   # one time step per bin
  binwidth_y <- diff(quantile(rsl_C$Catch, c(0.01, 0.99))) / 50  
  ggplot2::ggplot(rsl_C, ggplot2::aes(x=year, y=Catch)) +
    ggplot2::geom_bin2d(binwidth = c(binwidth_x, binwidth_y)) +
    ggplot2::geom_vline(xintercept=(HCR_sim$HCR$StartYear + HCR_sim$HCR$TN)) +
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
  Trig <- HCR_sim$ref_pt$B_tar_range[1]
  lvlBt <- c("B<Limit",   #Limit
             "B<Trigger", #Trigger
             "B~Target")  #Target
  lvlIt <- c("It<Limit",  #Limit
             "It<Trigger", #Trigger
             "It~Target")  #Target

  if (HCR_sim$stock_assessment == "fishblicc") {
    pBtar  <- with(HCR_sim, sweep(SSB, 1, STATS=ref_pt$B_tar, FUN="/"))
    pBlim  <- with(HCR_sim, sweep(SSB, 1, STATS=ref_pt$B_lim, FUN="/"))
    Index <- as.vector(HCR_sim$pjIndex)
  } else {
    pBtar <- with(HCR_sim, sweep(pB[ , (HCR_sim$HCR$TN+1L):(HCR_sim$HCR$PTN+1L)],
                               MARGIN=1, STATS=ref_pt$B_tar, FUN="/"))
    pBlim <- with(HCR_sim, sweep(pB[ , (HCR_sim$HCR$TN+1L):(HCR_sim$HCR$PTN+1L)],
                                 MARGIN=1, STATS=ref_pt$B_lim, FUN="/"))
    Index <- as.vector(HCR_sim$pjIndex[ , -1])
  }
  
  lo_trigger <- min(unlist(HCR_sim$HCR$trIndex))
  hi_trigger <- max(unlist(HCR_sim$HCR$trIndex))
  
  HCR_r <- tibble::tibble(State=lvlBt[3 - ((pBlim < 1) + (pBtar < Trig))],
                          Response=lvlIt[3 - ((Index < lo_trigger) +
                                                (Index < hi_trigger))]) |>
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
#' interval relative to MSY.
#'
#' @inheritParams HCR_performance
#' @return `flextable` containing the simulations status proportions for the HCR
#' @export
#'
table_sim_status <- function(HCR_sim) {
  std_brk <- with(HCR_sim$ref_pt, c(0.0, LRP, B_tar_range, 200.0))
  c_std_brk <- paste0(format(std_brk[-length(std_brk)], nsmall=1), "-",
                      format(std_brk[-1], nsmall=1))
  c_std_brk[length(c_std_brk)] <- paste0(">", format(std_brk[length(std_brk)-1L],
                                                     nsmall=1))
  
  # Projection statistics
  if (HCR_sim$stock_assessment == "fishblicc") {
    Ftar <- with(HCR_sim$HCR, matrix(0, nrow=nsim, ncol=PYN+1L))
    F_tar <- HCR_sim$ref_pt$F_tar
    NG <- ncol(F_tar)
    for (i in seq_len(NG))
      Ftar <- Ftar + sweep(HCR_sim$mF[ ,, i], 1, STATS=NG * F_tar[,i], FUN="/")
    pBtar  <- with(HCR_sim, sweep(SSB, 1, STATS=ref_pt$B_tar, FUN="/"))
  } else {
    #FMSY <- with(HCR_sim$Par, log(1 + r * (1-BMSY^(m-1)) / (m-1)))
    Ftar <- with(HCR_sim, as.vector(sweep(Ft[ ,(HCR$TN+1):HCR$PTN], MARGIN=1, STATS=ref_pt$F_tar, FUN="/")))
    pBtar <- with(HCR_sim, as.vector(sweep(pB[ ,(HCR$TN+1):HCR$PTN+1L], MARGIN=1, STATS=ref_pt$B_tar, FUN="/")))
  }
  
  cnt <- tibble::tibble(
    `B/B_target` = c_std_brk,
    `B%` = hist(pBtar, breaks=std_brk, plot=FALSE)$counts,
    `F%` = hist(Ftar, breaks=std_brk, plot=FALSE)$counts,
  )
  cnt$`B%` <- 100*cnt$`B%`/sum(cnt$`B%`)
  cnt$`F%` <- 100*cnt$`F%`/sum(cnt$`F%`)
  return(
    flextable::flextable(cnt) |>
      flextable::colformat_double(digits=1)
  )
}
