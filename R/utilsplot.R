
# additional functions and helpers for plots in gcplot.R

## provide vectors and labels for plots, by a selected function
make_plot_ypred <- function(
    t,
    curve,
    functions,
    dep="", time=""
) {
  
  
  if (curve == "growthFunction") {
    ypred <- functions$W(t)
    axisLabel <- paste0(dep)
    legend <- "Growth function"
    
  } else if (curve == "ogf") {
    ypred <- functions$OGF(t)
    axisLabel <- paste0("OGF (", dep, "² / ", time, "³)")
    legend <- "Ontogenetic growth force"
    
  } else if (curve == "growthRate") {
    ypred <- functions$W1(t)
    axisLabel <- paste0("Growth rate (", dep, " / ", time, ")")
    legend <- "Growth rate"
    
  } else if (curve == "acceleration") {
    ypred <- functions$W2(t)
    axisLabel <- paste0("Acceleration (", dep, " / ", time, "²)")
    legend <- "Acceleration"
    
  } else if (curve == "squaredRate") {
    ypred <- functions$W1s(t)
    axisLabel <- paste0("y'(t)² (", dep, "² / ", time, "⁴)")
    legend <- "Squared growth rate"
    
  } else if (curve == "gfSquaredRate") {
    ypred <- functions$W_W1s(t)
    axisLabel <- paste0("y(t) * y'(t)² (", dep, "³ / ", time, "⁵)")
    legend <- "y(t) * y'(t)²"
    
  } else {
    stop("ERROR: Invalid curvetion nmxion in make_plot_ypred()")
  }
  
  
  list(
    ypred = ypred,
    axisLabel = axisLabel,
    legend = legend
  )
}


## provide critical points (x, y) locations for plots
get_plot_cpoints <- function(
    op,
    yfunc,
    critPoints
) {
  
  cpLocs <- list()
  
  
  iPoint <- critPoints$iPoint
  
  if (op$iPoints) {
    cpLocs[["I"]] <- list(
      x = iPoint,
      y = yfunc(iPoint)
    )
  }
  
  lagPoints <- critPoints$lagPoints
  ogfPoints <- critPoints$ogfPoints
  accPoints <- critPoints$accPoints
  asyPoints <- critPoints$asyPoints
  
  cpSpecs <- list(
    list(group = "lagPoints", nmx = "ogf0", tloc = lagPoints$OGF0, label = "F0"),
    list(group = "lagPoints", nmx = "paa", tloc = lagPoints$PAA, label = "PAA"),
    list(group = "lagPoints", nmx = "tangent", tloc = lagPoints$tangent, label = "Tan."),
    list(group = "lagPoints", nmx = "threshold", tloc = lagPoints$threshold, label = "Thr."),
    list(group = "ogfPoints", nmx = "ogfMax", tloc = ogfPoints$F1, label = "F1"),
    list(group = "accPoints", nmx = "accMax", tloc = accPoints$P1, label = "P1"),
    list(group = "ogfPoints", nmx = "ogfMin", tloc = ogfPoints$F2, label = "F2"),
    list(group = "accPoints", nmx = "accMin", tloc = accPoints$P2, label = "P2"),
    list(group = "asyPoints", nmx = "ogf3", tloc = asyPoints$OGF3, label = "F3"),
    list(group = "asyPoints", nmx = "pda", tloc = asyPoints$PDA, label = "PDA")
  )
  
  for (cpSpec in cpSpecs) {
    
    if (!is.na(cpSpec$tloc) && cpSpec$nmx %in% op[[cpSpec$group]]) {
      
      cpLocs[[cpSpec$label]] <- list(
        x = cpSpec$tloc,
        y = yfunc(cpSpec$tloc)
      )
      
    }
    
  }
  
  
  cpLocs
}


## better limits to avoid "out-of-axis" points
better_lim <- function(
  vec, 
  n=5
) {
  
  ticks <- pretty(vec, n=n)
  
  c(min(ticks), max(ticks))
}


## draw a critical point inside plot area
draw_cpoint <- function(
  x,
  y,
  label = "",
  pcolor = "black",
  scolor = "slategrey"
) {
    
  segments(
    x, par("usr")[3], x, y,
    col = scolor, lwd = 2, lty = 4
  )
  
  points(x, y, pch = 19, col = pcolor)
    
  text(
    x, y, labels = label,
    pos = 3, offset = 1, col = pcolor, font = 2
  )
  
}


## draw a new curve inside plot area
draw_new_curve <- function(
  x,
  y,
  axisLabel = NULL,
  linetype = "dashed",
  linecolor = "black"
) {
  
  par(new = TRUE) # new plot, same area
  
  plot(x, y, type = "l", col = linecolor, axes = FALSE,
       xlab = "", ylab = "", lty = linetype, ylim = better_lim(y))
  
  if (!is.null(axisLabel))
    axis(side = 4) # add new y-axis on the right side of the plot
    mtext(axisLabel, side = 4, line = 2) # y-axis label
    
}

