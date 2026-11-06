
# additional functions and helpers for models specifications in specsmodel.R

## switch options name by model spec function
get_model_spec <- function(
    model,
    t,
    y,
    parameterization = "Ti",
    double = FALSE,
    fixValue = NULL
) {
  
  if (!model %in% c(
    "richards",
    "gompertz",
    "logistic",
    "bertalanffy"
  ))
    stop("ERROR: Invalid model family")
  
  spec_urichards(
    t = t,
    y = y,
    family = model,
    parameterization = parameterization,
    phases = if (isTRUE(double)) 2L else 1L,
    fixValue = fixValue
  )
}


## 
ur_family <- function(family) {
  
  switch(
    family,
    
    richards = list(
      name = "Richards",
      lambda = NULL
    ),
    
    gompertz = list(
      name = "Gompertz",
      lambda = 0
    ),
    
    logistic = list(
      name = "Logistic",
      lambda = log(2)
    ),
    
    bertalanffy = list(
      name = "von Bertalanffy",
      lambda = log(2 / 3)
    ),
    
    stop("Invalid U-Richards family")
  )
}


# Regular Richards core
pgm_R <- function(lambda, x) {
  u <- x * expm1(lambda)
  L <- ifelse(u == 0, 1, log1p(u) / u)
  exp(-x * L)
}


##
ur_h <- function(lambda, rho) {
  ell <- -log(rho)
  u <- ell * expm1(lambda)
  q <- ifelse(u == 0, 1, expm1(u) / u)
  ell * q
}


##
ur_x <- function(
    t,
    K,
    lambda,
    parameterization,
    Ti = NULL,
    rho = NULL
) {
  
  a <- expm1(lambda)
  p <- pgm_R(lambda, 1)
  c <- K * exp(lambda) / p
  
  if (parameterization == "Ti") {
    x <- exp(-c * (t - Ti))
  } else {
    h0 <- ur_h(lambda, rho)
    x <- h0 * exp(-c * t)
  }
  
  list(x = x, a = a, c = c)
}


##
ur_value <- function(
    t,
    A,
    K,
    lambda,
    parameterization,
    Ti = NULL,
    rho = NULL
) {
  
  s <- ur_x(
    t = t,
    K = K,
    lambda = lambda,
    parameterization = parameterization,
    Ti = Ti,
    rho = rho
  )
  
  A * pgm_R(lambda, s$x)
}


##
ur_derivatives <- function(
    t,
    A,
    K,
    lambda,
    parameterization,
    Ti = NULL,
    rho = NULL
) {
  
  s <- ur_x(
    t = t,
    K = K,
    lambda = lambda,
    parameterization = parameterization,
    Ti = Ti,
    rho = rho
  )
  
  x <- s$x
  a <- s$a
  c <- s$c
  
  W <- A * pgm_R(lambda, x)
  
  v  <- c * x / (1 + a * x)
  
  v1 <- v * (a * v - c)
  v2 <- (2 * a * v - c) * v1
  v3 <- 2 * a * v1^2 +
    (2 * a * v - c) * v2
  v4 <- 6 * a * v1 * v2 +
    (2 * a * v - c) * v3
  
  W1 <- W * v
  
  W2 <- W * (
    v^2 + v1
  )
  
  W3 <- W * (
    v^3 +
      3 * v * v1 +
      v2
  )
  
  W4 <- W * (
    v^4 +
      6 * v^2 * v1 +
      3 * v1^2 +
      4 * v * v2 +
      v3
  )
  
  W5 <- W * (
    v^5 +
      10 * v^3 * v1 +
      15 * v * v1^2 +
      10 * v^2 * v2 +
      10 * v1 * v2 +
      5 * v * v3 +
      v4
  )
  
  list(
    W = W,
    W1 = W1,
    W2 = W2,
    W3 = W3,
    W4 = W4,
    W5 = W5
  )
}

ur_derivative_functions <- function(
    beta,
    family,
    parameterization,
    fixValue = NULL
) {
  
  if ("A" %in% names(beta)) {
    
    A <- unname(beta["A"])
    
  } else {
    
    A <- suppressWarnings(
      as.numeric(fixValue)
    )
  }
  
  K <- unname(
    beta["K"]
  )
  
  lambda <- if ("d" %in% names(beta))
    log(unname(beta["d"]))
  else
    ur_family(family)$lambda
  
  Ti <- if (parameterization == "Ti")
    unname(beta["Ti"])
  else
    NULL
  
  rho <- if (parameterization == "W0")
    unname(beta["W0"]) / A
  else
    NULL
  
  
  calc <- function(t) {
    
    ur_derivatives(
      t = t,
      A = A,
      K = K,
      lambda = lambda,
      parameterization = parameterization,
      Ti = Ti,
      rho = rho
    )
  }
  
  
  list(
    
    W = function(t)
      calc(t)$W,
    
    W1 = function(t)
      calc(t)$W1,
    
    W1s = function(t) {
      z <- calc(t)
      z$W1^2
    },
    
    W_W1s = function(t) {
      z <- calc(t)
      z$W * z$W1^2
    },
    
    W2 = function(t)
      calc(t)$W2,
    
    W4 = function(t)
      calc(t)$W4,
    
    OGF = function(t) {
      z <- calc(t)
      z$W1 * z$W2
    },
    
    OGF3 = function(t) {
      z <- calc(t)
      
      3 * z$W3^2 +
        4 * z$W2 * z$W4 +
        z$W1 * z$W5
    }
  )
}

