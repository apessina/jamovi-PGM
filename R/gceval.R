
# functions for model evaluation in gccurve.b.R


## Natural-scale parameter estimates and fixed-effect contrasts
fitted_parameters <- function(
    fit,
    mspec,
    fixedTerms,
    data
) {
  
  ## Fixed coefficients in estimation scale
  betaTrans <- get_fixed_coef(fit)
  
  ## Covariance matrix in estimation scale
  covmTrans <- get_fixed_vcov(fit = fit, betaTrans = betaTrans)
  
  ## Fixed-effects inference
  infer <- get_fixed_inference(fit)
  
  ## Structural parameters present in this model
  transPars <- names(fixedTerms)
  
  ## Natural names
  naturalNames <- c(logA="A", logK="K", Ti="Ti", eta0="W0", lambda="d")
  
  ## Find intercept coefficient for each structural parameter
  interceptNames <- setNames(
    vapply(
      transPars,
      function(p) {
        
        if (p %in% names(betaTrans))
          return(p)
        
        intercept <- paste0(p, ".(Intercept)")
        
        if (!intercept %in% names(betaTrans)) {
          stop(
            "Could not identify the intercept for parameter ",
            p,
            ".",
            call. = FALSE
          )
        }
        
        intercept
      },
      character(1)
    ),
    transPars
  )
  
  ## Reference-profile structural parameters
  betaBase <- setNames(betaTrans[interceptNames], transPars)
  
  ## Natural-scale reference parameters
  betaHat <- mspec$backtransform(betaBase)
  
  ## 95% CI
  crit <- stats::qnorm(0.975)
  
  ## Evaluation context for parameter-table helpers
  evalContext <- list(
    betaTrans = betaTrans,
    covmTrans = covmTrans,
    infer = infer,
    mspec = mspec,
    transPars = transPars,
    interceptNames = interceptNames,
    naturalNames = naturalNames,
    fixedTerms = fixedTerms,
    data = data,
    crit = crit
  )
  
  ## Build parameter table
  rows <- list()
  
  for (p in transPars) {
    
    ## Reference parameter
    rows[[length(rows) + 1L]] <- make_row(p = p, context = evalContext)
    
    ## Fixed-effect contrasts
    prefix <- paste0(p, ".")
    
    effects <- names(betaTrans)[startsWith(names(betaTrans), prefix)]
    
    effects <- setdiff(effects, interceptNames[[p]])
    
    for (effect in effects) {
      rows[[length(rows) + 1L]] <- make_row(p, effect, context = evalContext)
    }
  }
  
  pardf <- do.call(rbind, rows)
  rownames(pardf) <- make.unique(pardf$Parameter)
  pardf$Parameter <- NULL
  
  list(
    pardf = pardf,
    beta = betaHat,
    betaTrans = betaBase,
    nfixed = length(betaTrans)
  )
}


# calculate requested GoF and Error metrics for a model
fittedmodel_evaluation <- function (
  op,
  fit,
  npar
) {
  
  
  ## Original-scale fitted values
  yhat <- get_fitted_response(fit)
  
  ## Original-scale fitted residuals
  res <- get_response_residuals(fit)
  resnorm <- get_normalized_residuals(fit)
  
  ## Observed response aligned with fitted/residual vectors
  yobs <- yhat + res
  
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
    SSt <- sum((yobs - mean(yobs))^2, na.rm = TRUE)
    SSe <- sum(res^2, na.rm = TRUE)
    
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
        mean(yobs, na.rm = TRUE)
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
          abs(yobs) + abs(yhat),
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
    
    ## Innermost grouping level
    diagLevel <- fit$dims$Q
    
    fittedValues <- as.numeric(fitted(fit, level = diagLevel))
    
    residualRaw <- as.numeric(
      residuals(fit, type = "response", level = diagLevel)
    )
    
    residualNorm <- as.numeric(
      residuals(fit, type = "normalized", level = diagLevel)
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
      "ERROR: Unsupported fitted model class: ",
      paste(class(fit), collapse = ", "),
      call. = FALSE
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


## Random-effects variance components
fitted_random_effects <- function(
    fit,
    randomGroups,
    randomGroupLabels,
    mspec,
    betaTrans
) {
  
  if (!inherits(fit, "nlme"))
    return(NULL)
  
  reStruct <-
    fit$modelStruct$reStruct
  
  if (
    is.null(reStruct) ||
    length(reStruct) == 0L
  ) {
    return(NULL)
  }
  
  
  naturalNames <- c(
    logA = "A",
    logK = "K",
    Ti = "Ti",
    eta0 = "W0",
    lambda = "d"
  )
  
  logPars <- c(
    "logA",
    "logK",
    "lambda"
  )
  
  
  ## Natural-scale reference profile
  naturalBase <-
    mspec$backtransform(
      betaTrans
    )
  
  ## Jacobian at reference profile
  Jbase <-
    mspec$jacobian(
      betaTrans
    )
  
  
  rows <- list()
  
  
  for (g in names(reStruct)) {
    
    groupLabel <- pgm_group_label(group = g, labels = randomGroupLabels)
    
    pars <- names(randomGroups)[
      vapply(
        randomGroups,
        function(x)
          g %in% x,
        logical(1)
      )
    ]
    
    if (!length(pars))
      next
    
    
    vc <- as.matrix(
      nlme::VarCorr(
        reStruct[[g]],
        sigma = fit$sigma
      )
    )
    
    n <- min(
      length(pars),
      nrow(vc)
    )
    
    
    for (i in seq_len(n)) {
      
      p <- pars[i]
      
      tau <- as.numeric(
        vc[i, "StdDev"]
      )
      
      natural <-
        naturalNames[[p]]
      
      
      ## log(A), log(K), lambda:
      ## exact induced lognormal SD
      if (p %in% logPars) {
        
        baseEstimate <-
          naturalBase[[natural]]
        
        sdNatural <-
          baseEstimate *
          exp(
            0.5 * tau^2
          ) *
          sqrt(
            expm1(
              tau^2
            )
          )
        
        
        ## eta0 = logit(W0/A):
        ## delta-method SD on W0 scale
      } else if (p == "eta0") {
        
        derivative <- as.numeric(
          Jbase[
            "W0",
            "eta0"
          ]
        )
        
        sdNatural <-
          abs(derivative) *
          tau
        
        
        ## Ti is already on natural scale
      } else {
        
        sdNatural <- tau
      }
      
      
      rows[[length(rows) + 1L]] <-
        data.frame(
          group = groupLabel,
          parameter = natural,
          sd = sdNatural,
          stringsAsFactors = FALSE
        )
    }
  }
  
  
  if (!length(rows))
    return(NULL)
  
  do.call(
    rbind,
    rows
  )
}

