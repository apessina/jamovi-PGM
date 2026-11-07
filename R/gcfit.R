
# functions for model fit in gccurve.b.R

## Fit first a nlsLM model and provide initial estimations for subsequent 
## models, when needed or requested
pgm_fit <- function(
    op,
    fitData,
    mexpr,
    mformula,
    tparinit,
    evartype,
    ecortype
) {

  ## Standard estimation methods
  estInfo <- "Nonlinear Least Squares"
  optInfo <- "Levenberg-Marquardt"
  randomInfo <- NULL
  evarEq <- "Homoscedastic errors"
  ecorEq <- "Independent errors"
  
  
  ## Initial fit using Levenberg-Marquardt -----
  fitLM <- tryCatch(
    
    minpack.lm::nlsLM(
      formula = mformula,
      data = fitData,
      start = tparinit # we use log scale for A and K
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
  
  
  ## Error weights -----
  if (evartype == "power") {
    evarfunc <- nlme::varPower(form = ~ fitted(.))
    evarEq <- "Var(εᵢ) = σ²|ŷᵢ|^(2θ)."
    
  } else if (evartype == "exp") {
    evarfunc <- nlme::varExp(form = ~ fitted(.))
    evarEq <- "Var(εᵢ) = σ²exp(2θŷᵢ)."
    
  } else if (evartype == "cpower") {
    evarfunc <- nlme::varConstPower(form = ~ fitted(.))
    evarEq <- "Var(εᵢ) = σ²[θ₁+|ŷᵢ|^(θ₂)]²."
    
  } else {
    evarfunc <- NULL
  }
  
  
  ## Model structures -----
  paramNames <- c(
    "logA",
    "logK",
    if (op$modelParam == "W0") "eta0" else "Ti",
    "lambda"
  )
  naturalNames <- setNames(
    c("A", "K", if (op$modelParam == "W0") "W0" else "Ti", "d"), paramNames
  )
  
  fixedTerms <- setNames(
    lapply(pgm_parameter_terms(op$fixedTerms), pgm_normalize_terms),
    paramNames
  )
  randomTerms <- setNames(
    lapply(pgm_parameter_terms(op$randomGroups), pgm_normalize_terms),
    paramNames
  )
  
  ## Keep only parameters present in this setup
  fixedTerms <- fixedTerms[
    intersect(names(fixedTerms), names(tparinit))
  ]
  randomTerms <- randomTerms[
    intersect(names(randomTerms), names(tparinit))
  ]
  
  randomPrepared <-
    pgm_prepare_random_groups(data = fitData, randomTerms = randomTerms)
  
  fitData <- randomPrepared$data
  randomGroups <- randomPrepared$groups
  randomGroupLabels <- randomPrepared$labels
  randomGroupComponents <- randomPrepared$components
  ## Match factor levels to the data effectively used by nlme/gnls
  fitData <- droplevels(fitData)
  pgm_validate_fixed_terms(data = fitData, fixedTerms = fixedTerms)
  
  corrGroups <- unique(
    as.character(unlist(op$corrGroups, use.names = FALSE))
  )
  
  hasFixed <- any(lengths(fixedTerms) > 0L)
  hasRandom <- any(lengths(randomGroups) > 0L)
  hasVariance <- evartype != "none"
  hasCorrelation <- ecortype != "none"
  
  if (hasCorrelation && !length(corrGroups))
    stop("AR/CAR requires at least one correlation grouping variable.")
  
  ### Fixed effects formula -----
  fixedFormula <- setNames(
    lapply(names(fixedTerms), function(p) {
      pgm_fixed_formula(parameter = p, terms = fixedTerms[[p]])
    }),
    names(fixedTerms)
  )
  pgm_validate_fixed_design(formulas = fixedFormula, data = fitData)
  fixedStart <- pgm_fixed_start(
    formulas = fixedFormula, data = fitData, baseStart = coef(fitLM)
  )
  
  fixedInfo <- paste(
    vapply(
      names(fixedTerms),
      function(p) {
        terms <- pgm_compose_terms(fixedTerms[[p]])
        rhs <- if (length(terms)) paste(terms, collapse = " + ") else "1"
        paste0(naturalNames[[p]], " ~ ", rhs)
      },
      character(1)
    ),
    collapse = "; "
  )
  
  if (hasRandom) {
    grouping <- pgm_group_hierarchy(
      data = fitData,
      randomGroups = randomGroups,
      corrGroups = if (hasCorrelation)
        corrGroups
      else
        character(0),
      groupComponents = randomGroupComponents
    )
    groupVars <- grouping$random
  } else {
    grouping <- NULL
    groupVars <- character(0)
  }
  
  
  ## Errors autocorrelation -----
  if (hasCorrelation) {
    
    corGroupVars <- if (hasRandom)
      grouping$correlation
    else
      corrGroups
    
    corGroupExpr <- paste(
      vapply(
        corGroupVars,
        function(x) deparse(as.name(x)),
        character(1)
      ),
      collapse = "/"
    )
    
    corFormula <- as.formula(
      paste("~ t |", corGroupExpr)
    )
    
    corGroupDisplay <- paste(
      vapply(
        corGroupVars,
        pgm_group_label,
        character(1),
        labels = randomGroupLabels
      ),
      collapse = "/"
    )
    
    if (ecortype == "ar1") {
      ecorfunc <- nlme::corAR1(form = corFormula)
      ecorEq <- "AR(1): Corr(εᵢt, εᵢ,t+h) = φ^|h|."
      
    } else {
      ecorfunc <- nlme::corCAR1(form = corFormula)
      ecorEq <- "CAR(1): Corr(εᵢj, εᵢk) = φ^(|t_ij - t_ik|)."
    }
    
  } else {
    ecorfunc <- NULL
  }
  
  
  ## Mixed effects -----
  if (hasRandom) {
    
    estInfo <- "Nonlinear Least Squares + NL Mixed Effects"
    optInfo <- "Levenberg-Marquardt + Gauss-Newton"
    
    ### Random effects structure -----
    randomStruct <- setNames(
      lapply(groupVars, function(g) {
        
        pars <- names(randomGroups)[
          vapply(randomGroups, function(x) g %in% x, logical(1))
        ]
        
        nlme::pdDiag(
          as.formula(paste(paste(pars, collapse = " + "), "~ 1"))
        )
      }),
      groupVars
    )
    
    ### Fit nlme -----
    fit <- tryCatch(
      
      nlme::nlme(
        model = mformula,
        data = fitData,
        fixed = fixedFormula,
        random = randomStruct,
        start = fixedStart,
        weights = evarfunc,
        correlation = ecorfunc,
        method = "ML"
      ),
      
      error = function(e) {
        stop(
          paste0(
            "Failed to fit nonlinear mixed effects model. ",
            "Consider changing terms or residuals structure. Original error: ",
            conditionMessage(e)
          ),
          call. = FALSE
        )
      }
      
    )
    
    randomInfo <- paste0(
      paste(
        vapply(
          groupVars,
          function(g) {
            pars <- names(randomGroups)[
              vapply(randomGroups, function(x) g %in% x, logical(1))
            ]
            groupLab <- pgm_group_label(group = g, labels = randomGroupLabels)
            paste0(
              paste(naturalNames[pars], collapse = " + "), " ~ 1 | ", groupLab)
          },
          character(1)
        ),
        collapse = "; "
      ), 
    " [diagonal covariance]")
    
    
    
  ## Generalized Nonlinear Least Squares -----
  } else if (hasFixed || hasVariance || hasCorrelation) {
    
    estInfo <- "Nonlinear Least Squares + Generalized NLS"
    optInfo <- "Levenberg-Marquardt + Gauss-Newton"
    
    ### Fit gnls -----
    fit <- tryCatch(
      
      nlme::gnls(
        model = mformula,
        data = fitData,
        params = fixedFormula,
        start = fixedStart,
        weights = evarfunc,
        correlation = ecorfunc
      ),
      
      error = function(e) {
        stop(
          paste0(
            "Failed to fit generalized nonlinear model. ",
            "Consider changing terms or residuals structure. Original error: ",
            conditionMessage(e)
          ),
          call. = FALSE
        )
      }
      
    )

    
  } else {
    
    fit <- fitLM
  }
  

  fmt <- function(x) {
    format(
      signif(as.numeric(x), 5),
      trim = TRUE,
      scientific = FALSE
    )
  }
  
  ## Residual scale
  sigmaHat <- tryCatch(
    as.numeric(stats::sigma(fit)),
    error = function(e) {
      if (!is.null(fit$sigma))
        as.numeric(fit$sigma)
      else
        NA_real_
    }
  )
  
  ## Variance-function estimates
  if (hasVariance) {
    
    varCoef <- coef(
      fit$modelStruct$varStruct,
      unconstrained = FALSE,
      allCoef = TRUE
    )
    
    if (evartype == "power") {
      
      evarEq <- paste0(
        evarEq,
        " σ = ", fmt(sigmaHat),
        "; θ = ", fmt(varCoef[1])
      )
      
    } else if (evartype == "exp") {
      
      evarEq <- paste0(
        evarEq,
        " σ = ", fmt(sigmaHat),
        "; θ = ", fmt(varCoef[1])
      )
      
    } else if (evartype == "cpower") {
      
      evarEq <- paste0(
        evarEq,
        " σ = ", fmt(sigmaHat),
        "; θ₁ = ", fmt(varCoef[1]),
        "; θ₂ = ", fmt(varCoef[2])
      )
    }
    
  } else if (is.finite(sigmaHat)) {
    
    evarEq <- paste0(
      evarEq,
      "; σ = ",
      fmt(sigmaHat)
    )
  }
  
  
  if (hasCorrelation) {
    
    corCoef <- coef(
      fit$modelStruct$corStruct,
      unconstrained = FALSE
    )
    
    ecorEq <- paste0(
      ecorEq,
      " Grouping: ~ t | ",
      corGroupDisplay,
      "; φ = ",
      fmt(corCoef[1])
    )
  }
  
  
  list(
    fit = fit,
    fitData = fitData,
    fixedTerms = fixedTerms,
    randomGroups = randomGroups,
    randomGroupLabels = randomGroupLabels,
    fixedInfo = fixedInfo,
    estInfo = estInfo,
    optInfo = optInfo,
    randomInfo = randomInfo,
    evarEq = evarEq,
    ecorEq = ecorEq
  )
}

