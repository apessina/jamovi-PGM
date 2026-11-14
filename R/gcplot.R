
# functions to populate plots states

## provide vectors for growthPlot state
populate_growthplot <- function(
    op,
    tCurves,
    xSamples,
    ySamples,
    functions,
    cp = NULL,
    dep = "", time = "",
    yLegend = NULL
) {
  
  growth <- make_plot_ypred(
    t = tCurves, 
    curve = "growthFunction",
    functions = functions,
    dep = dep, time = time
  )
  
  curves <- list(
    growth = list(
      x = tCurves,
      y = growth$ypred,
      axisLabel = growth$axisLabel,
      legend = if (!is.null(yLegend)) yLegend else growth$legend 
    )
  )
  
  if (op$growthRightAxis) {
    
    rightAxis <- make_plot_ypred(
      t = tCurves, 
      curve = op$growthRightAxisFunc,
      functions = functions,
      dep = dep, time = time
    )
    
    curves$rightAxis <- list(
      x = tCurves,
      y = rightAxis$ypred,
      axisLabel = rightAxis$axisLabel,
      legend = rightAxis$legend
    )
    
  }
  
  if (op$growthNoAxis) {
    
    noAxis <- make_plot_ypred(
      t = tCurves, 
      curve = op$growthNoAxisFunc,
      functions = functions
    )
    
    curves$noAxis <- list(
      x = tCurves,
      y = noAxis$ypred,
      legend = noAxis$legend
    )
    
  }
  
  if ("critPoints" %in% op$plotsElements && !is.null(cp)) {
    critPoints <- get_plot_cpoints(
      op = op,
      yfunc = functions$W,
      critPoints = cp
    )
  } else {
    critPoints <- NULL
  }
  
  list(
    timeLabel = time,
    curves = curves,
    cps = critPoints,
    samples = list(
      x = xSamples, 
      y = ySamples
    )
  )
}


## provide vectors for growthComparisonPlot global lists
prepare_comparison_data <- function(
    op,
    tCurves,
    functions,
    cp = NULL,
    dep = "",
    time = ""
) {
  
  ## Growth curve
  curves <- list(x = tCurves, growth = functions$W(tCurves))
  
  ## Right-axis function
  if (op$comparisonRightFunc != "none") {
    
    rightAxis <- make_plot_ypred(
      t = tCurves,
      curve = op$comparisonRightFunc,
      functions = functions,
      dep = dep,
      time = time
    )
    
    curves$rightAxis <- rightAxis$ypred
    curves$rightAxisLabel <- rightAxis$axisLabel
    curves$rightAxisLegend <- rightAxis$legend
    
  }
  
  ## Critical points
  if ("critPoints" %in% op$plotsElements && !is.null(cp)) {
    
    critPoints <- get_plot_cpoints(
      op = op,
      yfunc = functions$W,
      critPoints = cp
    )
    
  } else {
    critPoints <- NULL
  }
  
  list(
    curves = curves,
    cps = critPoints
  )
}


## Grouping unit for residual ACF
prepare_diagnostic_id <- function(
    op,
    groups,
    data,
    nlmegroups = NULL
) {
  diagnosticId <- NULL
  
  corrDiag <- unique(as.character(unlist(groups, use.names = FALSE)))
  
  ## residual-correlation grouping, when explicitly specified
  if (length(corrDiag) == 1L) {
    diagnosticId <- data[[corrDiag]]
    
  } else if (length(corrDiag) > 1L) {
    diagnosticId <- interaction(
      data[, corrDiag, drop = FALSE],
      drop = TRUE,
      sep = "/"
    )
    
    ## Otherwise use the innermost NLME grouping level
  } else if (!is.null(nlmegroups) && ncol(nlmegroups) > 0L) {
    diagnosticId <- nlmegroups[[ncol(nlmegroups)]]
  }
  
  diagnosticId
}
  
