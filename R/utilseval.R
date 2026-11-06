
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
get_fixed_vcov <- function(fit, betaTrans) {
  
  covm <- try(
    stats::vcov(fit),
    silent = TRUE
  )
  
  if (inherits(covm, "try-error")) {
    return(NULL)
  }
  
  covm <- as.matrix(covm)
  
  required <- names(betaTrans)
  
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


## Fixed-effects inference according to model type
get_fixed_inference <- function(fit) {
  
  ## NLME
  if (inherits(fit, "nlme")) {
    
    tab <- summary(fit)$tTable
    
    return(
      data.frame(
        estimate  = as.numeric(tab[, 1]),
        se        = as.numeric(tab[, 2]),
        df        = as.numeric(tab[, 3]),
        statistic = as.numeric(tab[, 4]),
        p         = as.numeric(tab[, 5]),
        row.names = rownames(tab)
      )
    )
  }
  
  
  ## GNLS
  if (inherits(fit, "gnls")) {
    
    tab <- summary(fit)$tTable
    
    df <- fit$dims[["N"]] -
      fit$dims[["p"]]
    
    return(
      data.frame(
        estimate  = as.numeric(tab[, 1]),
        se        = as.numeric(tab[, 2]),
        df        = rep(
          as.numeric(df),
          nrow(tab)
        ),
        statistic = as.numeric(tab[, 3]),
        p         = as.numeric(tab[, 4]),
        row.names = rownames(tab)
      )
    )
  }
  
  
  ## NLS / nlsLM
  if (inherits(fit, "nls")) {
    
    tab <- summary(fit)$coefficients
    
    return(
      data.frame(
        estimate  = as.numeric(tab[, 1]),
        se        = as.numeric(tab[, 2]),
        df        = rep(
          stats::df.residual(fit),
          nrow(tab)
        ),
        statistic = as.numeric(tab[, 3]),
        p         = as.numeric(tab[, 4]),
        row.names = rownames(tab)
      )
    )
  }
  
  
  stop(
    "ERROR: Unsupported fitted model class: ",
    paste(class(fit), collapse = ", "),
    call. = FALSE
  )
}


## get fixed y_hats according to model type
get_fitted_response <- function(fit) {
  
  if (inherits(fit, "nlme")) {
    level <- fit$dims$Q
    return(as.numeric(stats::fitted(fit, level = level)))
  }
  
  return(as.numeric(stats::fitted(fit)))
}


## get fixed residuals according to model type
get_response_residuals <- function(fit) {
  
  if (inherits(fit, "nlme")) {
    level <- fit$dims$Q
    return(as.numeric(stats::residuals(fit, level = level, type = "response")))
  }
  
  if (inherits(fit, "gnls")) {
    return(as.numeric(stats::residuals(fit, type = "response")))
  }
  
  return(as.numeric(stats::residuals(fit)))
}


## get fixed norm. residuals according to model type
get_normalized_residuals <- function(fit) {
  
  if (inherits(fit, "nlme")) {
    level <- fit$dims$Q
    return(as.numeric(stats::residuals(fit, level=level, type="normalized")))
  }
  
  if (inherits(fit, "gnls")) {
    return(as.numeric(stats::residuals(fit, type = "normalized")))
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
    lagMax = NULL
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
    
    localLagMax <- if (is.null(lagMax)) {
      min(10L, nrow(d) - 1L)
    } else {
      min(lagMax, nrow(d) - 1L)
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


## 
reference_profile <- function(
    p,
    context
) {
  
  vars <- pgm_term_variables(
    context$fixedTerms[[p]]
  )
  
  data <- context$data
  
  if (!length(vars))
    return(data.frame())
  
  values <- setNames(
    vector(
      "list",
      length(vars)
    ),
    vars
  )
  
  for (v in vars) {
    
    x <- data[[v]]
    
    if (is.factor(x)) {
      
      lev <- levels(x)
      
      if (!length(lev)) {
        stop(
          "Factor predictor has no levels: ",
          v,
          call. = FALSE
        )
      }
      
      values[[v]] <- factor(
        lev[1],
        levels = lev,
        ordered = is.ordered(x)
      )
      
    } else {
      
      values[[v]] <- 0
    }
  }
  
  as.data.frame(
    values,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}


##
profile_settings <- function(profile) {
  
  if (!ncol(profile))
    return(character(0))
  
  setNames(
    vapply(
      names(profile),
      function(v)
        as.character(
          profile[[v]][1]
        ),
      character(1)
    ),
    names(profile)
  )
}


##
effect_design_profile <- function(
    p,
    effect,
    context
) {
  
  fixedTerms <- context$fixedTerms
  data <- context$data
  
  formula <- pgm_fixed_formula(
    parameter = p,
    terms = fixedTerms[[p]]
  )
  
  termObj <- stats::delete.response(
    stats::terms(formula)
  )
  
  ref <- reference_profile(
    p = p,
    context = context
  )
  
  Xref <- stats::model.matrix(
    termObj,
    data = ref
  )
  
  mmName <- substring(
    effect,
    nchar(p) + 2L
  )
  
  column <- match(
    mmName,
    colnames(Xref)
  )
  
  if (is.na(column)) {
    stop(
      "Could not identify model-matrix column for fixed effect ",
      effect,
      ".",
      call. = FALSE
    )
  }
  
  assign <- attr(
    Xref,
    "assign"
  )[[column]]
  
  termLabels <- attr(
    termObj,
    "term.labels"
  )
  
  if (
    assign < 1L ||
    assign > length(termLabels)
  ) {
    stop(
      "Could not identify model term for fixed effect ",
      effect,
      ".",
      call. = FALSE
    )
  }
  
  components <- jmvcore::decomposeTerm(
    termLabels[[assign]]
  )
  
  candidates <- setNames(
    lapply(
      components,
      function(v) {
        
        x <- data[[v]]
        
        if (is.factor(x)) {
          
          lev <- levels(x)
          
          if (length(lev) < 2L)
            return(character(0))
          
          return(lev[-1])
        }
        
        1
      }
    ),
    components
  )
  
  if (any(lengths(candidates) == 0L)) {
    stop(
      "Could not construct predictor profile for fixed effect ",
      effect,
      ".",
      call. = FALSE
    )
  }
  
  grid <- expand.grid(
    candidates,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  selectedProfile <- NULL
  selectedDesign <- NULL
  
  for (i in seq_len(nrow(grid))) {
    
    profile <- ref
    
    for (v in components) {
      
      x <- data[[v]]
      value <- grid[[v]][i]
      
      if (is.factor(x)) {
        
        profile[[v]] <- factor(
          value,
          levels = levels(x),
          ordered = is.ordered(x)
        )
        
      } else {
        
        profile[[v]] <- as.numeric(
          value
        )
      }
    }
    
    X <- stats::model.matrix(
      termObj,
      data = profile
    )
    
    target <- X[
      1,
      mmName
    ]
    
    if (
      is.finite(target) &&
      abs(target) >
      sqrt(.Machine$double.eps)
    ) {
      
      selectedProfile <- profile
      selectedDesign <- X[
        1,
        ,
        drop = TRUE
      ]
      
      break
    }
  }
  
  if (is.null(selectedProfile)) {
    stop(
      "Could not construct predictor profile for fixed effect ",
      effect,
      ".",
      call. = FALSE
    )
  }
  
  list(
    design = selectedDesign,
    profile = selectedProfile
  )
}


##
build_profile <- function(
    p,
    effect = NULL,
    context
) {
  
  transPars <- context$transPars
  betaTrans <- context$betaTrans
  interceptNames <- context$interceptNames
  
  M <- matrix(
    0,
    nrow = length(transPars),
    ncol = length(betaTrans),
    dimnames = list(
      transPars,
      names(betaTrans)
    )
  )
  
  ## Reference intercept for every structural parameter
  for (q in transPars) {
    M[
      q,
      interceptNames[[q]]
    ] <- 1
  }
  
  profile <- reference_profile(
    p = p,
    context = context
  )
  
  if (!is.null(effect)) {
    
    effectProfile <-
      effect_design_profile(
        p = p,
        effect = effect,
        context = context
      )
    
    design <- effectProfile$design
    profile <- effectProfile$profile
    
    designNames <- vapply(
      names(design),
      function(x) {
        
        if (identical(
          x,
          "(Intercept)"
        )) {
          interceptNames[[p]]
          
        } else {
          paste0(
            p,
            ".",
            x
          )
        }
      },
      character(1)
    )
    
    if (!all(
      designNames %in%
      names(betaTrans)
    )) {
      stop(
        "Fixed-effect design does not match fitted coefficients for parameter ",
        p,
        ".",
        call. = FALSE
      )
    }
    
    M[p, ] <- 0
    
    M[
      p,
      designNames
    ] <- as.numeric(
      design
    )
  }
  
  x <- as.numeric(
    M %*% betaTrans
  )
  
  names(x) <- transPars
  
  list(
    x = x,
    M = M,
    settings =
      profile_settings(profile)
  )
}


##
reference_levels <- function(
    p,
    context
) {
  
  profile_settings(
    reference_profile(
      p = p,
      context = context
    )
  )
}


##
make_label <- function(
    p,
    settings = NULL,
    context
) {
  
  natural <-
    context$naturalNames[[p]]
  
  if (is.null(settings)) {
    settings <- reference_levels(
      p = p,
      context = context
    )
  }
  
  if (!length(settings))
    return(natural)
  
  paste0(
    natural,
    " [",
    paste(
      names(settings),
      settings,
      sep = " = ",
      collapse = ", "
    ),
    "]"
  )
}


## Build one natural-scale table row
make_row <- function(
    p,
    effect = NULL,
    context
) {
  
  betaTrans <- context$betaTrans
  covmTrans <- context$covmTrans
  infer <- context$infer
  mspec <- context$mspec
  transPars <- context$transPars
  naturalNames <- context$naturalNames
  crit <- context$crit
  
  profile <- build_profile(
    p = p,
    effect = effect,
    context = context
  )
  
  x <- profile$x
  M <- profile$M
  
  naturalName <- naturalNames[[p]]
  
  natural <- mspec$backtransform(
    x
  )
  
  estimate <- as.numeric(
    natural[[naturalName]]
  )
  
  se <- lower <- upper <- NA_real_
  
  if (!is.null(covmTrans)) {
    
    J <- mspec$jacobian(
      x
    )
    
    gBeta <- as.numeric(
      J[
        naturalName,
        transPars,
        drop = FALSE
      ] %*% M
    )
    
    names(gBeta) <- colnames(M)
    
    varNatural <- as.numeric(
      t(gBeta) %*%
        covmTrans %*%
        gBeta
    )
    
    se <- sqrt(
      pmax(
        varNatural,
        0
      )
    )
    
    
    ## W0 depends jointly on A and eta0
    if (naturalName == "W0") {
      
      if (
        is.finite(estimate) &&
        estimate > 0
      ) {
        
        gLog <- gBeta /
          estimate
        
        varLog <- as.numeric(
          t(gLog) %*%
            covmTrans %*%
            gLog
        )
        
        seLog <- sqrt(
          pmax(
            varLog,
            0
          )
        )
        
        lower <- estimate *
          exp(-crit * seLog)
        
        upper <- estimate *
          exp(crit * seLog)
      }
      
      
    } else {
      
      L <- M[
        p,
        ,
        drop = TRUE
      ]
      
      varEta <- as.numeric(
        t(L) %*%
          covmTrans %*%
          L
      )
      
      seEta <- sqrt(
        pmax(
          varEta,
          0
        )
      )
      
      xLower <- x
      xUpper <- x
      
      xLower[[p]] <-
        x[[p]] -
        crit * seEta
      
      xUpper[[p]] <-
        x[[p]] +
        crit * seEta
      
      lower <- as.numeric(
        mspec$backtransform(
          xLower
        )[[naturalName]]
      )
      
      upper <- as.numeric(
        mspec$backtransform(
          xUpper
        )[[naturalName]]
      )
    }
  }
  
  
  ## p-values only for fixed-effect contrasts
  pvalue <- NA_real_
  
  if (
    !is.null(effect) &&
    effect %in% rownames(infer)
  ) {
    
    pvalue <- infer[
      effect,
      "p"
    ]
  }
  
  
  data.frame(
    Parameter = make_label(
      p = p,
      settings = profile$settings,
      context = context
    ),
    parvalue = as.numeric(estimate),
    parerror = as.numeric(se),
    parlower = as.numeric(lower),
    parupper = as.numeric(upper),
    parp = as.numeric(pvalue),
    stringsAsFactors = FALSE
  )
}

