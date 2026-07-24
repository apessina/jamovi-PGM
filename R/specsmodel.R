
# build models equations. provide specs and complementary functions

## Richards - Tjørve and Tjørve, 2010 (Eq. 9)
spec_richards <- function(t, y) {
  
  ### Order observations
  ord <- order(t)
  t <- t[ord]
  y <- y[ord]
  
  
  ### Numerical growth rate
  dt <- diff(t)
  dy <- diff(y)
  
  validDiff <- (
    is.finite(dt) &
      is.finite(dy) &
      dt > 0
  )
  
  growthRate <- dy[validDiff] / dt[validDiff]
  growthTime <- (
    t[-1] + t[-length(t)]
  ) / 2
  growthTime <- growthTime[validDiff]
  
  
  ### Initial asymptote
  A0 <- max(y, na.rm = TRUE) * 1.10
  
  
  ### Initial inflection time
  if (
    length(growthRate) > 0 &&
    any(is.finite(growthRate))
  ) {
    
    Ti0 <- growthTime[
      which.max(growthRate)
    ]
    
  } else {
    
    Ti0 <- median(t, na.rm = TRUE)
  }
  
  
  ### Initial relative growth rate
  if (
    length(growthRate) > 0 &&
    any(is.finite(growthRate))
  ) {
    
    K0 <- max(
      growthRate,
      na.rm = TRUE
    ) / A0
    
  } else {
    
    K0 <- 1 / diff(range(t, na.rm = TRUE))
  }
  
  K0 <- max(
    K0,
    1e-4
  )
  
  
  ### Initial shape parameter
  d0 <- 1.5
  
  
  ### Model expression in natural parameterization
  mexpr <- expression(
    A * (1 + (d - 1) * exp(
      -K * (t - Ti) / d^(d / (1 - d))
    ))^(1 / (1 - d))
  )
  
  ### Initial values in natural parameterization
  parinit <- list(
    A  = A0,
    K  = K0,
    Ti = Ti0,
    d  = d0
  )
  
  ### Natural parameters restrictions
  parlower <- c(
    A = max(y, na.rm = TRUE),
    K = 1e-8,
    Ti = min(t, na.rm = TRUE),
    d = 1.05
  )
  
  parupper <- c(
    A = max(y, na.rm = TRUE) * 10,
    K = Inf,
    Ti = max(t, na.rm = TRUE),
    d = 20
  )
  
  
  ### Initial values for transformed fit
  tparinit <- list(
    logA = log(A0),
    logK = log(K0),
    Ti = Ti0,
    logdm105 = log(d0 - 1.05)
  )
  
  ### Transformed parameters restrictions
  tparlower <- c(
    logA = log(max(y, na.rm = TRUE)),
    logK = log(1e-8),
    Ti = min(t, na.rm = TRUE),
    logdm105 = -Inf
  )
  
  tparupper <- c(
    logA = log(max(y, na.rm = TRUE) * 10),
    logK = Inf,
    Ti = max(t, na.rm = TRUE),
    logdm105 = log(20 - 1.05)
  )
  
  ### Transformed formula
  mformula <- y ~
    exp(logA) *
    (
      1 +
        (0.05 + exp(logdm105)) *
        exp(
          -exp(logK) *
            (t - Ti) /
            (
              (1.05 + exp(logdm105)) ^
                (
                  (1.05 + exp(logdm105)) /
                    (1 - (1.05 + exp(logdm105)))
                )
            )
        )
    ) ^
    (
      1 /
        (1 - (1.05 + exp(logdm105)))
    )
  
  
  ### Back-transformation of transformed coef.s
  backtransform <- function(coef) {
    
    c(
      A = exp(unname(coef["logA"])),
      K = exp(unname(coef["logK"])),
      Ti = unname(coef["Ti"]),
      d = 1.05 + exp(
        unname(coef["logdm105"])
      )
    )
  }
  
  jacobian <- function(coef) {
    
    betaOriginal <- backtransform(coef)
    
    J <- matrix(
      0,
      nrow = 4,
      ncol = 4,
      dimnames = list(
        c("A", "K", "Ti", "d"),
        c(
          "logA",
          "logK",
          "Ti",
          "logdm105"
        )
      )
    )
    
    J["A", "logA"] <- betaOriginal["A"]
    J["K", "logK"] <- betaOriginal["K"]
    J["Ti", "Ti"] <- 1
    J["d", "logdm105"] <- betaOriginal["d"] - 1.05
    
    return(J)
  }
  
  
  mhtml <- paste0(
    '<i>A</i>(1 + (<i>d</i> − 1) exp(',
    '−<i>K</i>(<i>t</i> − <i>Ti</i>)/',
    '<i>d</i><sup><i>d</i>/(1 − <i>d</i>)</sup>))',
    '<sup>1/(1 − <i>d</i>)</sup>'
  )
  
  
  return(list(
    mname = "Richards",
    mhtml = mhtml,
    mexpr = mexpr,
    mformula = mformula,
    parinit = parinit,
    parlower = parlower,
    parupper = parupper,
    tparinit = tparinit,
    tparlower = tparlower,
    tparupper = tparupper,
    backtransform = backtransform,
    jacobian = jacobian
  ))
}

