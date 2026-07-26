
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
    tCurvesGlobal,
    tSamplesGlobal,
    ySamples,
    functions,
    cp = NULL
) {
  
  
  ySamples <- c(
    ySamples, rep(NA, length(tSamplesGlobal) - length(ySamples))
  )
  
  
  ## Fill length with NA to plot trimmed depSetupLabel lengths
  Wpred <- functions$W(tCurves)
  
  curves <- list(
    x = c(
      tCurves,
      rep(NA, length(tCurvesGlobal) - length(tCurves))
    ),
    growth = c(
      Wpred,
      rep(NA, length(tCurvesGlobal) - length(Wpred))
    )
  )
  
  
  if (op$comparisonRightFunc != "none") {
    
    raPred <- make_plot_ypred(
      t = tCurves,
      curve = op$comparisonRightFunc,
      functions = functions
    )$ypred
    
    curves$rightAxis <- c(
      raPred,
      rep(NA, length(tCurvesGlobal) - length(raPred))
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
    ySamples = ySamples,
    curves = curves,
    cps = critPoints
    )
}

