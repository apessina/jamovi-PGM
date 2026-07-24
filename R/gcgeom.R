
# functions related with growth curve geometry analysis

## creates a derivative tree structure to output derivatives functions  
## from the fitted growth curve expression
get_W_derivatives <- function(op, W_expr, t) {
  
  ## Conditions tree
  anyCritPoint <- any(
    op$iPoints,
    length(op$ogfPoints) > 0,
    length(op$accPoints) > 0,
    length(op$lagPoints) > 0,
    length(op$asyPoints) > 0
  )
  
  ogf3Conditions <- any(
    "ogf0" %in% op$lagPoints, 
    "ogf3" %in% op$asyPoints
  )
  
  ogfConditions <- any(
    ogf3Conditions,
    "ogfMax" %in% op$ogfPoints, 
    "ogfMin" %in% op$ogfPoints,
    "ogfI" %in% op$ogfPoints,
    op$comparisonRightFunc=="ogf",
    op$growthRightAxisFunc=="ogf", 
    op$growthNoAxisFunc=="ogf"
  )
  
  w4Conditions <- any(
    "paa" %in% op$lagPoints,
    "pda" %in% op$asyPoints
  )
  
  w2Conditions <- any(
    anyCritPoint,
    ogf3Conditions,
    ogfConditions,
    w4Conditions,
    op$comparisonRightFunc=="acceleration",
    op$growthRightAxisFunc=="acceleration", 
    op$growthNoAxisFunc=="acceleration"
  )
  
  w_w1sConditions <- any(
    "gfSquaredRate" %in% op$aucFunction, 
    op$comparisonRightFunc=="gfSquaredRate",
    op$growthRightAxisFunc=="gfSquaredRate", 
    op$growthNoAxisFunc=="gfSquaredRate"
  )
  
  w1sConditions <- any(
    w_w1sConditions,
    "squaredRate" %in% op$aucFunction, 
    op$comparisonRightFunc=="squaredRate",
    op$growthRightAxisFunc=="squaredRate", 
    op$growthNoAxisFunc=="squaredRate"
  )
  
  w1Conditions <- any(
    ogf3Conditions,
    ogfConditions,
    w4Conditions,
    w2Conditions,
    w1sConditions,
    op$comparisonRightFunc=="growthRate",
    op$growthRightAxisFunc=="growthRate", 
    op$growthNoAxisFunc=="growthRate"
  )
  
  
  ## Growth function
  W <- function(t) eval(W_expr, envir = list(t = t)) # W(t)
  
  ## Derivatives of Growth Function
  if (w1Conditions) {
    
    ### GF First derivative - Growth Rate
    W1_expr <- D(W_expr, "t") # 1st derivative
    W1 <- function(t) eval(W1_expr, envir = list(t = t)) # W'(t)
    
    ### Additional geometric functions based on W and W1
    if (w1sConditions) {
      W1s_expr <- call("^", W1_expr, 2) # W'(t)²
      W1s <- function(t) eval(W1s_expr, envir = list(t = t))
      
      if (w_w1sConditions) {
        W_W1s_expr <- call("*", W_expr, W1s_expr) # W(t) * W'(t)²
        W_W1s <- function(t) eval(W_W1s_expr, envir = list(t = t))
      }
      
    }
    
    ### GF Second derivative - Acceleration
    if (w2Conditions) {
      W2_expr <- D(W1_expr, "t") # 2nd derivative
      W2 <- function(t) eval(W2_expr, envir = list(t = t)) # W''(t)
      
      ## GF 4th derivative for PDA
      if (w4Conditions) {
        W3_expr <- D(W2_expr, "t") # 3rd derivative
        W4_expr <- D(W3_expr, "t") # 4th derivative
        W4 <- function(t) eval(W4_expr, envir = list(t = t)) # W''''(t) (PDA)
      }
      
      ## Ontogenetic Growth Force
      if (ogfConditions) {
        OGF_expr <- call("*", W1_expr, W2_expr) # W'(t) * W''(t)
        OGF <- function(t) eval(OGF_expr, envir = list(t = t)) # OGF(t)
        
        ## OGF 3rd derivative for F0 and F3
        if (ogf3Conditions) {
          OGF1_expr <- D(OGF_expr, "t") # 1st OGF derivative
          OGF2_expr <- D(OGF1_expr, "t") # 2nd OGF derivative 
          OGF3_expr <- D(OGF2_expr, "t") # 3rd OGF derivative 
          OGF3 <- function(t) eval(OGF3_expr, envir = list(t = t)) # OGF'''(t)
          
        }
        
      }
      
    }
    
  }
  
  
  rep_NA <- function(x) rep(NA_real_, length(x))
  
  list(
    W = W,
    W1 = if (exists("W1")) W1 else rep_NA, 
    W1s = if (exists("W1s")) W1s else rep_NA, 
    W_W1s = if (exists("W_W1s")) W_W1s else rep_NA, 
    W2 = if (exists("W2")) W2 else rep_NA, 
    W4 = if (exists("W4")) W4 else rep_NA, 
    OGF = if (exists("OGF")) OGF else rep_NA, 
    OGF3 = if (exists("OGF3")) OGF3 else rep_NA
  )
}


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
