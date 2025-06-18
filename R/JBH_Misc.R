



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
  obserr <- sqrt(dplyr::pull(jabba_fit$pars,
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

