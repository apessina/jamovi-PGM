
# functions related with growth curve geometry analysis


## Find first sign crossing
find_first_crossing <- function(values, start = 1) {
  
  values <- as.numeric(values)
  idx <- which(is.finite(values))
  
  if (length(idx) < 2)
    return(NA_integer_)
  
  crossing <- which(diff(sign(values)) != 0)
  crossing <- crossing[crossing >= start]
  
  if (length(crossing) == 0)
    return(NA_integer_)
  
  crossing[1]
}


## Critical points calculation
compute_critical_points <- function(
  op,
  t,
  Wderivs,
  thValue = NULL
) {
  
  ## Empty points structure
  note <- NULL
  iPoint <- NA_real_
  accPoints <- list(P1=NA_real_, P2=NA_real_)
  ogfPoints <- list(F1=NA_real_, F2=NA_real_)
  
  lagPoints <- list(
    OGF0=NA_real_, 
    PAA=NA_real_,
    tangent=NA_real_, 
    threshold=NA_real_
  )
  
  asyPoints <- list(
    OGF3=NA_real_, 
    PDA=NA_real_
  )
  
  
  ## Check if any critical point was requested
  anyCritPoint <- any(
    op$iPoints,
    length(op$accPoints) > 0,
    length(op$ogfPoints) > 0,
    length(op$lagPoints) > 0,
    length(op$asyPoints) > 0
  )
  
  if (anyCritPoint) {
    
  
    ## Check if the data has an inflection point
    W2_pred <- Wderivs$W2(t)
    i_Pi <- find_first_crossing(W2_pred)
    
    
    ## Calculate critical points
    if (is.finite(i_Pi)) {
      
      ### Inflection Point
      if (op$iPoints) {
        iPoint <- t[i_Pi]
      }
      
      
      ## Length of x-axis
      len <- length(t)
      
      ## Calculate OGF 3rd derivative vector
      if (any("ogf0" %in% op$lagPoints, "ogf3" %in% op$asyPoints)) {
        OGF3_pred <- Wderivs$OGF3(t)
      }
      
      ## Calculate growth function 4th derivative vector
      if (any("paa" %in% op$lagPoints, "pda" %in% op$asyPoints)) {
        W4_pred <- Wderivs$W4(t)
      }
      
      
      ## by Second Derivative
      
      ### P1
      needP1 <- any(
        "ogf0" %in% op$lagPoints,
        "paa" %in% op$lagPoints,
        "accMax" %in% op$accPoints
      )
      if (needP1) {
        
        i_P1 <- which.max(W2_pred[1:i_Pi])
        
        accPoints[["P1"]] <- t[i_P1]
      }  
      
      ### P2
      needP2 <- any(
        "ogf3" %in% op$asyPoints,
        "pda" %in% op$asyPoints,
        "accMin" %in% op$accPoints
      )
      if (needP2) {
        
        i_P2 <- i_Pi - 1 + which.min(W2_pred[i_Pi:len])
        
        accPoints[["P2"]] <- t[i_P2]
      }
      
      
      ## Effective Growth Onset
      
      ### OGF0
      if ("ogf0" %in% op$lagPoints) {
        
        if (i_P1 >= 2) {
          
          i_F0 <- which.max(OGF3_pred[1:i_P1])
          F0 <- t[i_F0]
          
          if (isTRUE(all.equal(F0, t[i_P1])) || isTRUE(all.equal(F0, min(t))))
            F0 <- NA_real_
        
        } else {
          F0 <- NA_real_
        }
        
        lagPoints[["OGF0"]] <- F0
      }
      
      ### PAA
      if ("paa" %in% op$lagPoints) {
        
        if (i_P1 >= 2) {
        
          i_PAA_rel <- find_first_crossing(W4_pred[1:i_P1])
          
          if (is.finite(i_PAA_rel)) {
            
            i_PAA <- i_PAA_rel
            PAA <- t[i_PAA]
            
          } else {
            PAA <- NA_real_
          }
          
          if (isTRUE(all.equal(PAA, t[i_P1])) || isTRUE(all.equal(PAA, min(t))))
            PAA <- NA_real_
          
        } else {
          PAA <- NA_real_
        }
        
        lagPoints[["PAA"]] <- PAA
      }
      
      ### Tangent method
      if ("tangent" %in% op$lagPoints) {
        
        W1_pred <- Wderivs$W1(t)
        
        i_star <- which.max(W1_pred)
        t_star <- t[i_star]
        
        slope <- W1_pred[i_star]
        y_star <- Wderivs$W(t_star)
        
        y0 <- Wderivs$W(min(t))
        
        if (is.finite(slope) && abs(slope) > .Machine$double.eps) {
          
          t_lag_tangent <- t_star - (y_star-y0) / slope
          
          if (isTRUE(all.equal(t_lag_tangent, min(t))))
            t_lag_tangent <- NA_real_
        
        } else {
          t_lag_tangent <- NA_real_
        }
        
        lagPoints[["tangent"]] <- t_lag_tangent
      }
      
      ### Threshold
      if ("threshold" %in% op$lagPoints) {
        
        if (!is.null(thValue) && is.finite(thValue)) {
          
          y_vals <- sapply(t, Wderivs$W)
          
          y_thresh <- Wderivs$W(max(t)) * thValue # custom
          i_thresh <- which(y_vals >= y_thresh)
          
          if (length(i_thresh) > 0) {
            t_lag_threshold <- t[min(i_thresh)]
            
            if (isTRUE(all.equal(t_lag_threshold, min(t))))
              t_lag_threshold <- NA_real_
            
          } else {
            t_lag_threshold <- NA_real_
          }
          
        } else {
          t_lag_threshold <- NA_real_
        }
        
        lagPoints[["threshold"]] <- t_lag_threshold
      }
      
      
      ## Close to Asymptote
      
      ### F3
      if ("ogf3" %in% op$asyPoints) {
        
        i_F3_rel <- find_first_crossing(OGF3_pred[i_P2:len])
        
        if (is.finite(i_F3_rel)) {
          
          i_F3 <- i_P2 - 1 + i_F3_rel
          F3 <- t[i_F3]
          
        } else {
          F3 <- NA_real_
        }
        
        asyPoints[["OGF3"]] <- F3
      }
            
      ### PDA
      if ("pda" %in% op$asyPoints) {
        
        i_PDA_rel <- find_first_crossing(W4_pred[i_P2:len])
        
        if (is.finite(i_PDA_rel)) {
          
          i_PDA <- i_P2 - 1 + i_PDA_rel
          PDA <- t[i_PDA]
          
        } else {
          PDA <- NA_real_
        }
        
        asyPoints[["PDA"]] <- PDA
      }
      
      
      ## Ontogenetic Growth Force
      
      if (any(c("ogfMax", "ogfMin") %in% op$ogfPoints)) { 
        OGF_pred <- Wderivs$OGF(t)
      }
      
      ### Fi - inflection point
      if ("ogfI" %in% op$ogfPoints) {
        
        ogfPoints[["Fi"]] <- t[i_Pi]
      }
      
      ### F1
      if ("ogfMax" %in% op$ogfPoints) {
        
        i_F1 <- which.max(OGF_pred[1:i_Pi])
        
        ogfPoints[["F1"]] <- t[i_F1]
      }
      
      ### F2
      if ("ogfMin" %in% op$ogfPoints) {
        
        i_F2 <- i_Pi - 1 + which.min(OGF_pred[i_Pi:len])
        
        ogfPoints[["F2"]] <- t[i_F2]
      }
      
  
    } else {
      note <- "No inflection point found"
    }
    
  
  }
  
  
  list(
    note = note,
    iPoint = iPoint,
    ogfPoints = ogfPoints,
    accPoints = accPoints,
    lagPoints = lagPoints,
    asyPoints = asyPoints
  )
}


## Area Under the Curve calculation
compute_auc <- function(
  op,
  t, 
  Wderivs
) {
  
  ## Geometric integral descriptors
  aucGf <- NA_real_
  aucSgr <- NA_real_
  aucGfSgr <- NA_real_
  
    
  # Integration interval
  if (op$aucInterval) {
    aucLower <- op$aucLower
    aucUpper <- op$aucUpper
  } else {
    aucLower <- min(t, na.rm=TRUE)
    aucUpper <- max(t, na.rm=TRUE)
  }
  
  
  ## Growth function
  if ("growthFunction" %in% op$aucFunction) {
    aucGf <- integrate(
      Wderivs$W, lower = aucLower, upper = aucUpper
    )$value
  }
  
  ## Square of the growth rate
  if ("squaredRate" %in% op$aucFunction) {
    aucSgr <- integrate(
      Wderivs$W1s, lower = aucLower, upper = aucUpper
    )$value
  }
  
  ## Weighted square of the growth rate
  if ("gfSquaredRate" %in% op$aucFunction) {
    aucGfSgr <- integrate(
      Wderivs$W_W1s, lower = aucLower, upper = aucUpper
    )$value
  }

  
  list(
    interval = list(
      aucLower = aucLower,
      aucUpper = aucUpper
    ),
    metrics = list(
      aucGf    = aucGf,
      aucSgr   = aucSgr,
      aucGfSgr = aucGfSgr
    )
  )
}
