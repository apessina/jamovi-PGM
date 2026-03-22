
# build models equations and provide its specs

## Richards - Tjørve & Tjørve (Eq. 9) ----
spec_richards <- function(t, y) {
  html_eq <- paste0(
    '<i>A</i><sub>1</sub>(1 + (<i>d</i><sub>1</sub> − 1) exp(',
    '−<i>K</i><sub>1</sub>(<i>t</i> − <i>T</i><sub>i1</sub>)/',
    '<i>d</i><sub>1</sub><sup><i>d</i><sub>1</sub>/(1 − <i>d</i><sub>1</sub>)</sup>))',
    '<sup>1/(1 − <i>d</i><sub>1</sub>)</sup>'
  )
  
  a_sup <- max(y, na.rm = TRUE)
  a_inf <- min(y, na.rm = TRUE)
  gr <- diff(y) / diff(t)
  
  expr <- richards_component_expr(1)
  
  init <- list(
    A1  = a_sup,
    K1  = max(gr, na.rm = TRUE) / max(a_sup, 1e-6),
    Ti1 = t[which.max(gr)],
    d1  = 1.2
  )
  
  lower <- c(A1 = 0, K1 = 0, Ti1 = min(t), d1 = 1.05)
  upper <- c(A1 = Inf, K1 = Inf, Ti1 = max(t), d1 = Inf)
  
  return(list(
    html_eq = html_eq, expr = expr, init = init, 
    lower = lower, upper = upper, n_comp = 1
  ))
}

## Double Richards (2 components) ----
spec_drichards <- function(t, y) {

  a_sup <- max(y, na.rm = TRUE)
  a_inf <- min(y, na.rm = TRUE)
  gr <- diff(y) / diff(t)  
    
  n_comp <- 2
  comps <- lapply(seq_len(n_comp), richards_component_expr)
  expr <- sum_exprs(comps)
  
  # heuristic init.
  t_min <- min(t); t_max <- max(t)
  t_span <- t_max - t_min
  ti_guess <- seq(
    t_min + 0.25 * t_span,
    t_max - 0.25 * t_span,
    length.out = n_comp
  )
  
  A_guess <- rep(a_sup / n_comp, n_comp)
  K_guess <- rep(max(gr, na.rm = TRUE) / max(a_sup, 1e-6), n_comp)
  d_guess <- rep(1.2, n_comp)
  
  init <- list()
  lower <- c()
  upper <- c()
  
  for (j in seq_len(n_comp)) {
    init[[paste0("A", j)]]  <- A_guess[j]
    init[[paste0("K", j)]]  <- K_guess[j]
    init[[paste0("Ti", j)]] <- ti_guess[j]
    init[[paste0("d", j)]]  <- d_guess[j]
    
    lower[paste0("A", j)]  <- 0
    lower[paste0("K", j)]  <- 0
    lower[paste0("Ti", j)] <- t_min
    lower[paste0("d", j)]  <- 1.05
    
    upper[paste0("A", j)]  <- Inf
    upper[paste0("K", j)]  <- Inf
    upper[paste0("Ti", j)] <- t_max
    upper[paste0("d", j)]  <- Inf
  }
  
  html_terms <- vapply(seq_len(n_comp), function(j) {
    paste0(
      '<i>A</i><sub>', j, '</sub>(1 + (<i>d</i><sub>', j, '</sub> − 1) exp(',
      '−<i>K</i><sub>', j, '</sub>(<i>t</i> − <i>T</i><sub>i', j, '</sub>)/',
      '<i>d</i><sub>', j, '</sub><sup><i>d</i><sub>', j, '</sub>/(1 − <i>d</i><sub>', j, '</sub>)</sup>))',
      '<sup>1/(1 − <i>d</i><sub>', j, '</sub>)</sup>'
    )
  }, character(1))
  html_eq <- paste(html_terms, collapse = " + ")
  
  return(list(
    html_eq = html_eq, expr = expr, init = init, 
    lower = lower, upper = upper, n_comp = n_comp
  ))
}

## Triple Richards (3 components) ----
spec_trichards <- function(t, y) {

  a_sup <- max(y, na.rm = TRUE)
  a_inf <- min(y, na.rm = TRUE)
  gr <- diff(y) / diff(t)  
    
  n_comp <- 3
  comps <- lapply(seq_len(n_comp), richards_component_expr)
  expr <- sum_exprs(comps)
  
  # heuristic init.
  t_min <- min(t); t_max <- max(t)
  t_span <- t_max - t_min
  ti_guess <- seq(
    t_min + 0.25 * t_span,
    t_max - 0.25 * t_span,
    length.out = n_comp
  )
  
  A_guess <- rep(a_sup / n_comp, n_comp)
  K_guess <- rep(max(gr, na.rm = TRUE) / max(a_sup, 1e-6), n_comp)
  d_guess <- rep(1.2, n_comp)
  
  init <- list()
  lower <- c()
  upper <- c()
  
  for (j in seq_len(n_comp)) {
    init[[paste0("A", j)]]  <- A_guess[j]
    init[[paste0("K", j)]]  <- K_guess[j]
    init[[paste0("Ti", j)]] <- ti_guess[j]
    init[[paste0("d", j)]]  <- d_guess[j]
    
    lower[paste0("A", j)]  <- 0
    lower[paste0("K", j)]  <- 0
    lower[paste0("Ti", j)] <- t_min
    lower[paste0("d", j)]  <- 1.05
    
    upper[paste0("A", j)]  <- Inf
    upper[paste0("K", j)]  <- Inf
    upper[paste0("Ti", j)] <- t_max
    upper[paste0("d", j)]  <- Inf
  }
  
  html_terms <- vapply(seq_len(n_comp), function(j) {
    paste0(
      '<i>A</i><sub>', j, '</sub>(1 + (<i>d</i><sub>', j, '</sub> − 1) exp(',
      '−<i>K</i><sub>', j, '</sub>(<i>t</i> − <i>T</i><sub>i', j, '</sub>)/',
      '<i>d</i><sub>', j, '</sub><sup><i>d</i><sub>', j, '</sub>/(1 − <i>d</i><sub>', j, '</sub>)</sup>))',
      '<sup>1/(1 − <i>d</i><sub>', j, '</sub>)</sup>'
    )
  }, character(1))
  html_eq <- paste(html_terms, collapse = " + ")
  
  return(list(
    html_eq = html_eq, expr = expr, init = init, 
    lower = lower, upper = upper, n_comp = n_comp
  ))
}

## Generalized Logistic - Koya & Goshu ----
spec_glogistic <- function(t, y) {
  html_eq <- paste0(
    '<i>A</i><sub>L</sub> + (<i>A</i> − <i>A</i><sub>L</sub>)',
    '[1 − <i>B</i>e<sup>−<i>k</i>(<i>t</i> − <i>&mu;</i>)</sup>]<sup><i>m</i></sup>'
  )

  expr <- substitute(
    AL + (A - AL) * (1 - B * exp(-k * (t - mu)))^m
  )
  
  t_min <- min(t, na.rm = TRUE)
  t_max <- max(t, na.rm = TRUE)
  y_min <- min(y, na.rm = TRUE)
  y_max <- max(y, na.rm = TRUE)
  t_span <- max(t_max - t_min, 1e-6)
  
  k_guess <- 1 / t_span
  if (!is.finite(k_guess) || k_guess <= 0)
    k_guess <- 0.1
  
  init <- list(
    AL = y_min,
    A  = y_max,
    B  = 0.5,
    k  = k_guess,
    mu = t_min,
    m  = 1
  )
  
  lower <- c(
    AL = -Inf,
    A  = 0,
    B  = 1e-6,
    k  = 1e-6,
    mu = t_min - t_span,
    m  = 0.1
  )
  
  upper <- c(
    AL = y_max,
    A  = Inf,
    B  = 0.999,
    k  = Inf,
    mu = t_min,
    m  = 5
  )
  
  return(list(
    html_eq = html_eq, expr = expr, init = init,
    lower = lower, upper = upper, n_comp = 1
  ))
}

## Generalized Weibull - Koya & Goshu ----
spec_gweibull <- function(t, y) {
  html_eq <- paste0(
    '<i>A</i>[1 − <i>B</i>e<sup>−<i>k</i>(',
    '((<i>t</i> − <i>&mu;</i>)/<i>&delta;</i>)<sup><i>&nu;</i></sup>',
    ')</sup>]'
  )
  
  expr <- substitute(
    A * (1 - B * exp(-k * (((t - mu) / delta)^nu)))
  )
  
  t_min <- min(t, na.rm = TRUE)
  t_max <- max(t, na.rm = TRUE)
  y_max <- max(y, na.rm = TRUE)
  t_span <- max(t_max - t_min, 1e-6)
  
  k_guess <- 1 / t_span
  if (!is.finite(k_guess) || k_guess <= 0)
    k_guess <- 0.1
  
  init <- list(
    A     = y_max,
    B     = 0.5,
    k     = k_guess,
    mu    = t_min - 0.25 * t_span,
    delta = max(t_span / 2, 1e-3),
    nu    = 1
  )
  
  lower <- c(
    A     = 0,
    B     = 1e-6,
    k     = 1e-6,
    mu    = t_min - 2 * t_span,
    delta = 1e-4,
    nu    = 0.5
  )
  
  upper <- c(
    A     = Inf,
    B     = 0.999,
    k     = Inf,
    mu    = t_min - 1e-6,
    delta = Inf,
    nu    = 3
  )
  
  return(list(
    html_eq = html_eq, expr = expr, init = init,
    lower = lower, upper = upper, n_comp = 1
  ))
}

## Get selected model specs ----
get_model_spec = function(model, t, y) {
  switch(
    model,
    richards = spec_richards(t, y),
    glogistic = spec_glogistic(t, y),
    gweibull = spec_gweibull(t, y),
    drichards = spec_drichards(t, y),
    trichards = spec_trichards(t, y),
    stop("ERROR: Invalid model equation")
  )
}