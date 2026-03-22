
# helpers to build models equations

## Expression for each Richards component ----
richards_component_expr = function(idx) {
  A  <- as.name(paste0("A", idx))
  K  <- as.name(paste0("K", idx))
  Ti <- as.name(paste0("Ti", idx))
  d  <- as.name(paste0("d", idx))
  
  ### Tjørve & Tjørve (Eq. 9)
  substitute(
    A * (1 + (d - 1) * exp(-K * (t - Ti) / (d^(d / (1 - d)))))^(1 / (1 - d)),
    list(A = A, K = K, Ti = Ti, d = d)
  )
}

## Concatenate Richards components ----
sum_exprs = function(expr_list) {
  Reduce(function(a, b) call("+", a, b), expr_list)
}

## Build single equation with Richards components ----
build_richards_n = function(n_comp, t, y) {
  comps <- lapply(seq_len(n_comp), richards_component_expr)
  expr <- sum_exprs(comps)
  
  list(expr = expr)
}