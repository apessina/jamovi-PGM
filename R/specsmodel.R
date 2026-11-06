
# build models equations. provide specs and complementary functions

## Ti U-Richards - Tjørve and Tjørve, 2010 (Eq. 9)
spec_urichards <- function(
    t,
    y,
    family = "richards",
    parameterization = "Ti",
    phases = 1L,
    fixValue = NULL
) {
  fam <- ur_family(family)
  
  lambdaFixed <- fam$lambda
  dFree <- is.null(lambdaFixed)
  
  ### Order observations
  ord <- order(t)
  t <- t[ord]
  y <- y[ord]
  
  
  ### Initial asymptote
  fixValue <- trimws(if (is.null(fixValue)) "" else fixValue)
  fixA <- nzchar(fixValue)
  
  if (fixA) {
    fixValue <- suppressWarnings(as.numeric(fixValue))
    if (!is.finite(fixValue) || fixValue <= 0)
      stop("The fixed asymptote must be a positive number")
    
    A0 <- fixValue
    
  } else {
    
    A0 <- max(y, na.rm = TRUE) * 1.10
  }
  
  
  ### Initial inflection time
  ss <- smooth.spline(t, y)
  
  dy <- predict(
    ss,
    x = t,
    deriv = 1
  )$y
  
  iMax <- which.max(dy)
  
  Ti0 <- t[iMax]
  
  
  ### Initial relative growth rate
  K0 <- max(
    max(dy, na.rm = TRUE) / A0,
    1e-4
  )
  
  
  ### Initial shape parameter
  lambda0 <- if (is.null(lambdaFixed)) {
    0
  } else {
    lambdaFixed
  }
  d0 <- exp(lambda0)
  
  
  ### 
  p0 <- pgm_R(lambda0, 1)
  c0 <- K0 * exp(lambda0) / p0
  
  rhoFirst <- y[1] / A0
  rhoFirst <- min(
    max(rhoFirst, 1e-8),
    1 - 1e-8
  )
  
  hFirst <- ur_h(
    lambda0,
    rhoFirst
  )
  
  h00 <- hFirst * exp(
    c0 * t[1]
  )
  
  a0 <- expm1(lambda0)
  
  if (a0 < 0) {
    hmax <- -1 / a0
    
    h00 <- min(
      h00,
      hmax * (1 - 1e-8)
    )
  }
  
  rho0 <- pgm_R(
    lambda0,
    h00
  )
  
  rho0 <- min(
    max(rho0, 1e-8),
    1 - 1e-8
  )
  
  W00 <- A0 * rho0
  eta00 <- qlogis(rho0)
  
  ### Model expression in natural parameterization
  ANat <- if (fixA) fixValue else quote(A)
  lambdaNat <- if (dFree) quote(log(d)) else lambdaFixed
  
  if (parameterization == "Ti") {
    rhsNat <- bquote(
      ur_value(
        t,
        .(ANat),
        K,
        .(lambdaNat),
        "Ti",
        Ti = Ti
      )
    )
  } else {
    rhsNat <- bquote(
      ur_value(
        t,
        .(ANat),
        K,
        .(lambdaNat),
        "W0",
        rho = W0 / .(ANat)
      )
    )
  }
  mexpr <- as.expression(rhsNat)
  
  ### Initial values for transformed fit
  if (parameterization == "Ti") {
    tparinit <- list(logA = log(A0), logK = log(K0), Ti = Ti0)
  } else {
    tparinit <- list(logA = log(A0), logK = log(K0), eta0 = eta00)
  }
  if (is.null(lambdaFixed)) {
    tparinit$lambda <- lambda0
  }
  
  ### Transformed parameters restrictions
  tparlower <- setNames(rep(-Inf, length(tparinit)), names(tparinit))
  tparupper <- setNames(rep(Inf, length(tparinit)), names(tparinit))
  
  if ("logA" %in% names(tparinit)) {
    tparlower["logA"] <- log(max(y, na.rm = TRUE))
    tparupper["logA"] <- log(max(y, na.rm = TRUE) * 10)
  }
  tparlower["logK"] <- log(1e-8)
  
  if ("Ti" %in% names(tparinit)) {
    tparlower["Ti"] <- min(t, na.rm = TRUE)
    tparupper["Ti"] <- max(t, na.rm = TRUE)
  }
  
  if (fixA) {
    tparinit$logA <- NULL
    tparlower <- tparlower[names(tparlower) != "logA"]
    tparupper <- tparupper[names(tparupper) != "logA"]
  }
  
  ### Transformed formula
  AFit <- if (fixA) fixValue else quote(exp(logA))
  
  ## Richards: lambda is estimated
  if (dFree) {
    if (parameterization == "Ti") {
      rhsFit <- bquote(
        .(AFit) *
          ifelse(
            (lambda + 0 * t) == 0,
            ## Gompertz limit
            exp(
              -exp(
                -exp(logK + 1) *
                  (t - Ti)
              )
            ),
            ## Richards, lambda != 0
            exp(
              -log1p(
                exp(
                  -exp(logK) *
                    (t - Ti) *
                    exp(
                      lambda +
                        lambda / expm1(lambda)
                    )
                ) *
                  expm1(lambda)
              ) /
                expm1(lambda)
            )
          )
      )
    } else {
      rhsFit <- bquote(
        .(AFit) *
          ifelse(
            (lambda + 0 * t) == 0,
            ## Gompertz W0 limit
            exp(
              -log1p(exp(-eta0)) *
                exp(
                  -exp(logK + 1) * t
                )
            ),
            ## Richards W0, lambda != 0
            exp(
              -log1p(
                expm1(
                  expm1(lambda) *
                    log1p(exp(-eta0))
                ) *
                  exp(
                    -exp(logK) *
                      t *
                      exp(
                        lambda +
                          lambda / expm1(lambda)
                      )
                  )
              ) /
                expm1(lambda)
            )
          )
      )
    }
  ## Gompertz: lambda = 0 fixed
  } else if (lambdaFixed == 0) {
    if (parameterization == "Ti") {
      rhsFit <- bquote(
        .(AFit) *
          exp(
            -exp(
              -exp(logK + 1) *
                (t - Ti)
            )
          )
      )
    } else {
      rhsFit <- bquote(
        .(AFit) *
          exp(
            -log1p(exp(-eta0)) *
              exp(
                -exp(logK + 1) * t
              )
          )
      )
    }
  ## Logistic / von Bertalanffy: lambda fixed and nonzero
  } else {
    aFixed <- expm1(lambdaFixed)
    cFactorFixed <- exp(lambdaFixed + lambdaFixed / aFixed)
    if (parameterization == "Ti") {
      rhsFit <- bquote(
        .(AFit) *
          exp(
            -log1p(
              .(aFixed) *
                exp(
                  -exp(logK) *
                    (t - Ti) *
                    .(cFactorFixed)
                )
            ) /
              .(aFixed)
          )
      )
    } else {
      rhsFit <- bquote(
        .(AFit) *
          exp(
            -log1p(
              expm1(
                .(aFixed) *
                  log1p(exp(-eta0))
              ) *
                exp(
                  -exp(logK) *
                    t *
                    .(cFactorFixed)
                )
            ) /
              .(aFixed)
          )
      )
    }
  }
  mformula <- as.formula(call("~", quote(y), rhsFit), env = environment())
  
  
  ### Back-transformation of transformed coef.s
  backtransform <- function(coef) {
    
    A <- if (fixA) fixValue else exp(unname(coef["logA"]))
    K <- exp(unname(coef["logK"]))
    
    out <- numeric()
    if (!fixA) out["A"] <- A
    out["K"] <- K
    
    if (parameterization == "Ti") {
      out["Ti"] <- unname(coef["Ti"])
    } else {
      rho <- plogis(unname(coef["eta0"]))
      out["W0"] <- A * rho
    }
    
    if (dFree) out["d"] <- exp(unname(coef["lambda"]))
    
    out
  }
  
  jacobian <- function(coef) {
    
    betaOriginal <- backtransform(coef)
    
    J <- matrix(
      0,
      nrow = length(betaOriginal),
      ncol = length(coef),
      dimnames = list(
        names(betaOriginal),
        names(coef)
      )
    )
    
    if (!fixA) J["A", "logA"] <- betaOriginal["A"]
    
    J["K", "logK"] <- betaOriginal["K"]
    
    if (parameterization == "Ti") {
      J["Ti", "Ti"] <- 1
    } else {
      A <- if (fixA) fixValue else betaOriginal["A"]
      rho <- plogis(unname(coef["eta0"]))
      if (!fixA)  J["W0", "logA"] <- betaOriginal["W0"]
      J["W0", "eta0"] <- A * rho * (1 - rho)
    }
    
    if (dFree) J["d", "lambda"] <- betaOriginal["d"]
    
    J
  }
  
  
  return(list(
    mname = "Richards",
    mexpr = mexpr,
    mformula = mformula,
    tparinit = tparinit,
    tparlower = tparlower,
    tparupper = tparupper,
    backtransform = backtransform,
    jacobian = jacobian
  ))
}

