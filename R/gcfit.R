
# functions for model fit in gccurve.b.R

## Fit first a nlsLM model and provide initial estimations for subsequent 
## models, when needed or requested
pgm_fit <- function(
  op,
  t,
  y,
  iid = NULL,
  mexpr,
  mformula,
  parinit,
  tparinit,
  evartype,
  ecortype
) {
  
  
  ## Structured data for modeling
  fitData <- data.frame(
    t = as.numeric(t),
    y = as.numeric(y)
  )
  if (!is.null(iid)) fitData$iid <- as.factor(iid)
  
  
  ## Check if initial values return finite output
  startEnv <- c(list(t = fitData$t), parinit)
  startVal <- eval(mexpr, envir = startEnv)
  
  if (any(!is.finite(startVal))) {
    stop("ERROR: Initial values generate non-finite model predictions")
  }
  
  
  ## Standard estimation methods
  estInfo <- "Nonlinear Least Squares"
  optInfo <- "Levenberg-Marquardt"
  randomEff <- NULL
  evarEq <- "Homoscedastic errors"
  ecorEq <- "Independent errors"
  
  
  ## Initial fit using Levenberg-Marquardt
  ## with transformed parameters
  fitLM <- tryCatch(
    
    minpack.lm::nlsLM(
      formula = mformula,
      data = fitData,
      start = tparinit
    ),
    
    error = function(e) {
      stop(
        paste0(
          "Failed to fit initial nonlinear model. ",
          "Consider checking the data structure. Original error: ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
    
  )
  
  
  ## Error weights
  if (evartype == "power") {
    
    evarfunc <- nlme::varPower(
      form = ~ fitted(.)
    )
    
    evarEq <- "Var(εᵢ) = σ²|ŷᵢ|^(2θ)."
    
  } else if (evartype == "exp") {
    
    evarfunc <- nlme::varExp(
      form = ~ fitted(.)
    )
    
    evarEq <- "Var(εᵢ) = σ²exp(2θŷᵢ)."
    
  } else if (evartype == "cpower") {
    
    evarfunc <- nlme::varConstPower(
      form = ~ fitted(.)
    )
    
    evarEq <- "Var(εᵢ) = σ²[θ₁+|ŷᵢ|^(θ₂)]²."
    
  } else {
    
    evarfunc <- NULL
  }
  
  ## Generalized Nonlinear Least Squares
  if (is.null(iid) & !is.null(evarfunc)) {
    
    estInfo <- "Nonlinear Least Squares + Generalized NLS"
    optInfo <- "Levenberg-Marquardt + Gauss-Newton"
    
    ## Fit gnls
    fit <- tryCatch(
      
      nlme::gnls(
        model = mformula,
        data = fitData,
        start = coef(fitLM),
        weights = evarfunc
      ),
      
      error = function(e) {
        stop(
          paste0(
            "Failed to fit generalized nonlinear model. ",
            "Consider changing variance function. Original error: ",
            conditionMessage(e)
          ),
          call. = FALSE
        )
      }
      
    )
    
    
  ## Mixed Effects with option AR
  ## when iid is provided
  } else if (!is.null(iid)) {
    
    estInfo <- "Nonlinear Least Squares + NL Mixed Effects"
    optInfo <- "Levenberg-Marquardt + Gauss-Newton"
    
    ## Errors autocorrelation
    if (ecortype=="ar1") {
      
      ecorfunc <- nlme::corAR1(
        form = ~ t | iid
      )
      
      ecorEq <- "AR(1): Corr(εᵢt, εᵢ,t+h) = φ^|h|."
      
    } else if (ecortype=="car1") {
      
      ecorfunc <- nlme::corCAR1(
        form = ~ t | iid
      )
      
      ecorEq <- "CAR(1): Corr(εᵢj, εᵢk) = φ^(|t_ij - t_ik|)."
      
    } else {
      
      ecorfunc <- NULL
    }
    
    
    ## Mixed Effect on the Asymptote
    fit <- nlme::nlme(
      model = mformula,
      data = fitData,
      fixed =
        logA +
        logK +
        Ti +
        logdm105 ~ 1,
      random =
        logA ~ 1 | iid,
      start = coef(fitLM),
      weights = evarfunc,
      correlation = ecorfunc,
      method = "ML"
    )
    
    randomEff <- "All Fixed + 'A' Random"
    
  } else {
    
    fit <- fitLM
  }
  

  list(
    fit = fit,
    estInfo = estInfo,
    optInfo = optInfo,
    randomEff = randomEff,
    evarEq = evarEq,
    ecorEq = ecorEq
  )
}

