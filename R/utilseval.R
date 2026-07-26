
# helpers for fitted model evaluation in gceval.R

## get fixed coefficients according to model type
get_fixed_coef <- function(fit) {
  
  if (inherits(fit, "nlme")) {
    
    beta <- nlme::fixed.effects(fit)
    
  } else if (
    inherits(fit, "gnls") ||
    inherits(fit, "nls")
  ) {
    
    beta <- stats::coef(fit)
    
  } else {
    
    stop("ERROR: Unsupported fitted model class")
  }
  
  
  betaNames <- names(beta)
  beta <- as.numeric(beta)
  names(beta) <- betaNames
  
  return(beta)
}


## get fixed covariance matrix according to model type
get_fixed_vcov <- function(fit, beta_trans) {
  
  covm <- try(
    stats::vcov(fit),
    silent = TRUE
  )
  
  if (inherits(covm, "try-error")) {
    return(NULL)
  }
  
  covm <- as.matrix(covm)
  
  required <- names(beta_trans)
  
  if (!all(required %in% rownames(covm))) {
    return(NULL)
  }
  
  covm <- covm[
    required,
    required,
    drop = FALSE
  ]
  
  return(covm)
}


## get fixed y_hats according to model type
get_fitted_response <- function(fit) {
  
  if (inherits(fit, "nlme")) {
    
    return(
      as.numeric(
        stats::fitted(
          fit,
          level = 1
        )
      )
    )
  }
  
  return(
    as.numeric(
      stats::fitted(fit)
    )
  )
}


## get fixed residuals according to model type
get_response_residuals <- function(fit) {
  
  if (inherits(fit, "nlme")) {
    
    return(
      as.numeric(
        stats::residuals(
          fit,
          level = 1,
          type = "response"
        )
      )
    )
  }
  
  if (inherits(fit, "gnls")) {
    
    return(
      as.numeric(
        stats::residuals(
          fit,
          type = "response"
        )
      )
    )
  }
  
  return(
    as.numeric(
      stats::residuals(fit)
    )
  )
}


## get fixed norm. residuals according to model type
get_normalized_residuals <- function(fit) {
  
  if (inherits(fit, "nlme")) {
    
    return(
      as.numeric(
        stats::residuals(
          fit,
          level = 1,
          type = "normalized"
        )
      )
    )
  }
  
  if (inherits(fit, "gnls")) {
    
    return(
      as.numeric(
        stats::residuals(
          fit,
          type = "normalized"
        )
      )
    )
  }
  
  res <- as.numeric(
    stats::residuals(fit)
  )
  
  sigmaHat <- try(
    stats::sigma(fit),
    silent = TRUE
  )
  
  if (
    inherits(sigmaHat, "try-error") ||
    !is.finite(sigmaHat) ||
    sigmaHat <= 0
  ) {
    
    return(
      rep(
        NA_real_,
        length(res)
      )
    )
  }
  
  return(res / sigmaHat)
}


## provide pooled ACF for autocorrelation diagnostic 
pooled_subject_acf <- function(
    residuals,
    time,
    id,
    lag_max = NULL
) {
  
  dat <- data.frame(
    residual = residuals,
    time = time,
    id = id
  )
  
  dat <- dat[
    is.finite(dat$residual) &
      is.finite(dat$time) &
      !is.na(dat$id),
    ,
    drop = FALSE
  ]
  
  splittedData <- split(dat, dat$id)
  
  subjectAcf <- lapply(splittedData, function(d) {
    
    d <- d[order(d$time), , drop = FALSE]
    
    if (nrow(d) < 4L)
      return(NULL)
    
    localLagMax <- if (is.null(lag_max)) {
      min(10L, nrow(d) - 1L)
    } else {
      min(lag_max, nrow(d) - 1L)
    }
    
    acfResult <- stats::acf(
      d$residual,
      lag.max = localLagMax,
      plot = FALSE,
      na.action = stats::na.pass
    )
    
    data.frame(
      lag = as.numeric(acfResult$lag),
      acf = as.numeric(acfResult$acf)
    )
  })
  
  subjectAcf <- Filter(Negate(is.null), subjectAcf)
  
  if (length(subjectAcf) == 0L)
    return(NULL)
  
  allAcf <- do.call(
    rbind,
    Map(
      function(x, subject) {
        x$subject <- subject
        x
      },
      subjectAcf,
      names(subjectAcf)
    )
  )
  
  aggregate(
    acf ~ lag,
    data = allAcf,
    FUN = mean,
    na.rm = TRUE
  )
}

