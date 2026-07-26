
# functions for model evaluation in gccurve.b.R

## dataframe with estimated parameters and error in original scale
fitted_parameters <- function (
  fit,
  mspec
) {
  
  
  ## Fixed effects parameters in transformed scale
  beta_trans <- get_fixed_coef(fit)
  
  ## Fixed effects parameters in natural scale
  beta_hat <- mspec$backtransform(beta_trans)
  
  
  ## Covariance matrix in transformed scale
  covm_trans <- get_fixed_vcov(
    fit = fit,
    beta_trans = beta_trans
  )
  
  
  ## Critical value
  crit <- stats::qnorm(
    1 - 0.05 / 2
  )
  
  ## Standard Error for parameters
  if (is.null(covm_trans)) {
    
    se_trans <- rep(
      NA_real_,
      length(beta_trans)
    )
    names(se_trans) <- names(beta_trans)
    
    se_vals <- rep(
      NA_real_,
      length(beta_hat)
    )
    names(se_vals) <- names(beta_hat)
    
    lower_ci <- rep(
      NA_real_,
      length(beta_hat)
    )
    names(lower_ci) <- names(beta_hat)
    
    upper_ci <- lower_ci
    
  } else {
    
    ### Jacobian of the back-transformation
    J <- mspec$jacobian(beta_trans)
    
    J <- J[
      names(beta_hat),
      names(beta_trans),
      drop = FALSE
    ]
    
    ### Natural-scale covariance matrix
    covm <- J %*%
      covm_trans %*%
      t(J)
    
    ### Standard errors
    se_trans <- sqrt(
      pmax(
        diag(covm_trans),
        0
      )
    )
    
    se_vals <- sqrt(
      pmax(
        diag(covm),
        0
      )
    )
    
    
    ### Confidence limits in transformed scale
    lower_trans <- beta_trans - crit * se_trans
    upper_trans <- beta_trans + crit * se_trans
    
    ### Back-transformed confidence limits
    lower_ci <- mspec$backtransform(lower_trans)
    upper_ci <- mspec$backtransform(upper_trans)
    
  }
  
  
  list(
    pardf   = data.frame(
      row.names = names(beta_hat),
      parvalue  = as.numeric(beta_hat),
      parerror  = as.numeric(se_vals),
      parlower  = as.numeric(lower_ci),
      parupper  = as.numeric(upper_ci)
    )
  )
}


# calculate requested GoF and Error metrics for a model
fittedmodel_evaluation <- function (
  op,
  y,
  fit,
  npar
) {
  
  
  ## Original-scale fitted values
  yhat <- get_fitted_response(fit)
  
  ## Original-scale fitted residuals
  res <- get_response_residuals(fit)
  resnorm <- get_normalized_residuals(fit)
  
  ## Number of observations
  n <- stats::nobs(fit)
  
  
  ## GoF
  AIC <- NA_real_
  BIC <- NA_real_
  AICc <- NA_real_
  R2 <- NA_real_
  R2adj <- NA_real_
  
  if (any(c("r2", "r2adj") %in% op$gofMetrics)) {
    
    ## Residual degrees of freedom
    df2 <- n - npar
    
    ## Sum of squares
    SSt <- sum(
      (y - mean(y))^2,
      na.rm = TRUE
    )
    
    SSe <- sum(
      res^2,
      na.rm = TRUE
    )
    
    ## Descriptive goodness-of-fit
    R2 <- 1 - SSe / SSt
    
    if ("r2adj" %in% op$gofMetrics) {
      if (df2 > 0 && n > 1) {
        R2adj <- 1 -
          (SSe / df2) /
          (SSt / (n - 1))
      } else {
        R2adj <- NA_real_
      }
    }
    
  }

  
  ## Information criteria
  if ("bic" %in% op$gofMetrics) {
    BIC <- stats::BIC(fit)
  }
  
  if (any(c("aic", "aicc") %in% op$gofMetrics)) {
    
    AIC <- stats::AIC(fit)
    
    if ("aicc" %in% op$gofMetrics) {
      
      ## Likelihood
      nparTotal <- attr(
        stats::logLik(fit), "df"
      )
      
      if (
        is.finite(nparTotal) &&
        n > nparTotal + 1
      ) {
        AICc <- AIC +
          (
            2 *
              nparTotal *
              (nparTotal + 1)
          ) /
          (
            n -
              nparTotal -
              1
          )
      } else {
        AICc <- NA_real_
      }
    
    }
  }
  
  
  ## Error metrics
  RMSE <- NA_real_
  MAE <- NA_real_
  MedAE <- NA_real_
  sMAPE <- NA_real_
  RRMSE <- NA_real_
  
  if (any(c("rmse", "rrmse") %in% op$errorMetrics)) {
    
    MSE <- mean(
      res^2,
      na.rm = TRUE
    )
    RMSE <- sqrt(MSE)
    
    if ("rrmse" %in% op$errorMetrics) {
      RRMSE <- RMSE /
        mean(y, na.rm = TRUE)
    }
    
  }
  
  if ("mae" %in% op$errorMetrics) {
    MAE <- mean(
      abs(res),
      na.rm = TRUE
    )
  }
  
  if ("medae" %in% op$errorMetrics) {
    MedAE <- stats::median(
      abs(res),
      na.rm = TRUE
    )
  }
  
  if ("smape" %in% op$errorMetrics) {
    sMAPE <- mean(
      2 * abs(res) /
        pmax(
          abs(y) + abs(yhat),
          1e-8
        ),
      na.rm = TRUE
    ) * 100
  }
  
  
  list(
    gof = list(
      AIC   = AIC,
      BIC   = BIC,
      AICc  = AICc,
      R2    = R2,
      R2adj = R2adj
    ),
    em = list(
      RMSE  = RMSE,
      MAE   = MAE,
      MedAE = MedAE,
      sMAPE = sMAPE,
      RRMSE = RRMSE
    )
  )
}


## Model residuals by fit type for diagnostics
residual_diagnostics <- function(
  fit, 
  time, 
  id = NULL
) {
  
  stopifnot(length(time) > 0L)
  
  
  ## NLME: conditional diagnostics
  if (inherits(fit, "nlme")) {
    
    fittedValues <- as.numeric(
      fitted(fit, level = 1)
    )
    
    residualRaw <- as.numeric(
      residuals(
        fit,
        type = "response",
        level = 1
      )
    )
    
    residualNorm <- as.numeric(
      residuals(
        fit,
        type = "normalized",
        level = 1
      )
    )
    
    fitType <- "nlme"
    
  ## GNLS
  } else if (inherits(fit, "gnls")) {
    
    fittedValues <- as.numeric(fitted(fit))
    
    residualRaw <- as.numeric(
      residuals(fit, type = "response")
    )
    
    residualNorm <- as.numeric(
      residuals(fit, type = "normalized")
    )
    
    fitType <- "gnls"
    
  ## NLS
  } else if (inherits(fit, "nls")) {
    
    fittedValues <- as.numeric(fitted(fit))
    residualRaw <- as.numeric(residuals(fit))
    
    dfResidual <- tryCatch(
      stats::df.residual(fit),
      error = function(e) {
        length(residualRaw) -
          length(stats::coef(fit))
      }
    )
    
    sigmaHat <- sqrt(
      sum(residualRaw^2, na.rm = TRUE) /
        max(dfResidual, 1L)
    )
    
    residualNorm <- residualRaw / sigmaHat
    fitType <- "nls"
    
  } else {
    stop(
      ERROR: "Unsupported fitted model class: ",
      paste(class(fit), collapse = ", ")
    )
  }
  
  
  n <- min(
    length(time),
    length(fittedValues),
    length(residualRaw),
    length(residualNorm)
  )
  
  output <- data.frame(
    time = as.numeric(time)[seq_len(n)],
    fittedValues = fittedValues[seq_len(n)],
    residualNorm = residualNorm[seq_len(n)]
  )
  
  if (!is.null(id)) {
    output$id <- as.factor(id[seq_len(n)])
  }
  
  output
}

