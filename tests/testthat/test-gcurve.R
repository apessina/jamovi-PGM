# Helpers ---------------------------------------------------------------

make_setup <- function(var, model = "logistic") {
  list(
    setupId = "test_1",
    setupLabel = "M1",
    var = var,
    model = model,
    errorVarType = "none",
    errorCorType = "none"
  )
}


empty_parameter_terms <- function(parameterization = "Ti") {
  list(
    list(param = "A", terms = list()),
    list(param = "K", terms = list()),
    list(
      param = if (parameterization == "W0") "W0" else "Ti",
      terms = list()
    ),
    list(param = "d", terms = list())
  )
}


# Explicit U-Logistic in Ti parameterization:
logistic_ti <- function(t, A, K, Ti) {
  A / (1 + exp(-4 * K * (t - Ti)))
}


make_logistic_data <- function() {
  
  t <- seq(0, 20, by = 0.5)
  
  y <- logistic_ti(
    t = t,
    A = 100,
    K = 0.08,
    Ti = 10
  )
  
  # Small deterministic perturbation prevents an exact zero-residual fit
  y <- y + 0.10 * sin(1.7 * t)
  
  data.frame(
    time = t,
    weight = y
  )
}


# 1. Basic NLS ----------------------------------------------------------

test_that("basic Logistic Ti model recovers known parameters", {
  
  data <- make_logistic_data()
  
  results <- PGM::gcurve(
    data = data,
    deps = "weight",
    time = "time",
    modelParam = "Ti",
    modeling = list(
      make_setup("weight")
    ),
    fixedTerms = empty_parameter_terms("Ti"),
    randomGroups = empty_parameter_terms("Ti")
  )
  
  par <- results$parTable$asDF
  
  expect_equal(
    par$Parameter,
    c("A", "K", "Ti")
  )
  
  expect_true(
    all(is.finite(par$Estimate))
  )
  
  expect_equal(
    par$Estimate[par$Parameter == "A"],
    100,
    tolerance = 1
  )
  
  expect_equal(
    par$Estimate[par$Parameter == "K"],
    0.08,
    tolerance = 0.01
  )
  
  expect_equal(
    par$Estimate[par$Parameter == "Ti"],
    10,
    tolerance = 0.2
  )
})


# 2. W0 parameterization ------------------------------------------------

test_that("W0 parameterization returns the expected initial value", {
  
  data <- make_logistic_data()
  
  results <- PGM::gcurve(
    data = data,
    deps = "weight",
    time = "time",
    modelParam = "W0",
    modeling = list(
      make_setup("weight")
    ),
    fixedTerms = empty_parameter_terms("W0"),
    randomGroups = empty_parameter_terms("W0")
  )
  
  par <- results$parTable$asDF
  
  expect_equal(
    par$Parameter,
    c("A", "K", "W0")
  )
  
  expectedW0 <- logistic_ti(
    t = 0,
    A = 100,
    K = 0.08,
    Ti = 10
  )
  
  expect_equal(
    par$Estimate[par$Parameter == "A"],
    100,
    tolerance = 1
  )
  
  expect_equal(
    par$Estimate[par$Parameter == "K"],
    0.08,
    tolerance = 0.01
  )
  
  expect_equal(
    par$Estimate[par$Parameter == "W0"],
    expectedW0,
    tolerance = 0.3
  )
})


# 3. Fixed effects / GNLS -----------------------------------------------

test_that("parameter-specific fixed effects use the GNLS path", {
  
  t <- seq(0, 20, by = 0.5)
  
  dataA <- data.frame(
    time = t,
    weight = logistic_ti(t, A = 100, K = 0.08, Ti = 10) +
      0.10 * sin(1.7 * t),
    treatment = factor("A")
  )
  
  dataB <- data.frame(
    time = t,
    weight = logistic_ti(t, A = 120, K = 0.08, Ti = 10) +
      0.10 * cos(1.3 * t),
    treatment = factor("B")
  )
  
  data <- rbind(dataA, dataB)
  data$treatment <- factor(data$treatment)
  
  fixed <- empty_parameter_terms("Ti")
  
  # Fixed effect on A only
  fixed[[1]]$terms <- list(
    list("treatment")
  )
  
  results <- PGM::gcurve(
    data = data,
    deps = "weight",
    time = "time",
    estTables = TRUE,
    modelParam = "Ti",
    modeling = list(
      make_setup("weight")
    ),
    fixedTerms = fixed,
    randomGroups = empty_parameter_terms("Ti")
  )
  
  est <- results$estTable$asDF
  par <- results$parTable$asDF
  
  estimator <- est$method[
    est$category == "Estimator"
  ]
  
  expect_match(
    estimator,
    "Generalized NLS"
  )
  
  # A should have a reference estimate plus a treatment contrast/profile
  expect_gt(
    sum(startsWith(par$Parameter, "A")),
    1
  )
  
  expect_true(
    all(is.finite(par$Estimate))
  )
})


# 4. Random effects / NLME ----------------------------------------------

test_that("parameter-specific random effects use the NLME path", {
  
  t <- seq(0, 20, by = 1)
  
  subjectA <- c(
    88, 92, 96, 100, 104,
    108, 112, 116, 120, 124
  )
  
  pieces <- lapply(
    seq_along(subjectA),
    function(i) {
      
      data.frame(
        time = t,
        weight =
          logistic_ti(
            t,
            A = subjectA[i],
            K = 0.08,
            Ti = 10
          ) +
          0.15 * sin(t + i),
        id = factor(
          paste0("S", i)
        )
      )
    }
  )
  
  data <- do.call(rbind, pieces)
  data$id <- factor(data$id)
  
  random <- empty_parameter_terms("Ti")
  
  # Random intercept for transformed A (logA) by individual
  random[[1]]$terms <- list(
    list("id")
  )
  
  results <- PGM::gcurve(
    data = data,
    deps = "weight",
    time = "time",
    estTables = TRUE,
    modelParam = "Ti",
    modeling = list(
      make_setup("weight")
    ),
    fixedTerms = empty_parameter_terms("Ti"),
    randomGroups = random
  )
  
  est <- results$estTable$asDF
  ran <- results$randomTable$asDF
  
  estimator <- est$method[
    est$category == "Estimator"
  ]
  
  expect_match(
    estimator,
    "NL Mixed Effects"
  )
  
  expect_true(
    nrow(ran) >= 1
  )
  
  expect_true(
    any(ran$Parameter == "A")
  )
  
  expect_true(
    all(
      is.finite(ran$SD) &
        ran$SD >= 0
    )
  )
})