# =========================================================
# functions related with critical points calculation
# =========================================================

compute_critical_points <- function(
    t_new,
    W1_pred,
    W2_pred,
    OGF_pred,
    OGF3_pred,
    W, W2, W4, 
    OGF, OGF3,
    thVal = NULL
) {
    
  ## Length of x-axis
  len <- length(t_new)
  
  ## Check if the data have an inflection point
  zero_acc <- which(diff(sign(W2_pred))!=0)
  if (!length(zero_acc)==0) {
    inflection = TRUE
    
    ## by Ontogenetic Growth Force
    ### Fi - inflection point
    i_Fi <- which(diff(sign(OGF_pred))!=0)[1]
    Fi <- t_new[i_Fi]
    ### F1
    i_F1 <- which.max(OGF(t_new[1:i_Fi]))
    F1 <- t_new[i_F1]
    ### F2
    i_F2 <- i_Fi - 1 + which.min(OGF(t_new[i_Fi:len]))
    F2 <- t_new[i_F2]
    ### List of calculated F-Points
    f_points <- list(F1=F1, Fi=Fi, F2=F2)
    
    ## by Growth Rate and Acceleration
    ### Pi
    i_Pi <- zero_acc[1]
    Pi <- t_new[i_Pi]
    ### P1
    i_P1 <- which.max(W2(t_new[1:i_Pi]))
    P1 <- t_new[i_P1]
    ### P2
    i_P2 <- i_Pi - 1 + which.min(W2(t_new[i_Pi:len]))
    P2 <- t_new[i_P2]
    ### List of calculated P-Points
    p_points <- list(P1=P1, Pi=Pi, P2=P2)
    
    ## End of Lag Phase
    ### OGF0 - end of lag phase
    i_F0 <- which.max(OGF3(t_new[1:i_F1]))
    F0 <- t_new[i_F0]
    if (F0==F1 | F0==0)
      F0 <- NA
    ### tangent method
    t_star <- t_new[which.max(W1_pred)]
    slope <- max(W1_pred)
    y_star <- W(t_star)
    y0 <- W(min(t_new)) # min(W_fun(t_new))
    t_lag_tangent <- t_star - (y_star-y0) / slope
    if (t_lag_tangent==0)
      t_lag_tangent <- NA
    ### threshold method
    y_vals <- sapply(t_new, W)
    y_thresh <- W(max(t_new)) * thVal ## custom by user
    t_lag_threshold <- t_new[min(which(y_vals >= y_thresh))]
    if (t_lag_threshold==0)
      t_lag_threshold <- NA
    ### List of calculated points
    l_points <- list(OGF0=F0, tang=t_lag_tangent, thres=t_lag_threshold)
    
    ## Close to Asymptote
    ### F3
    i_F3 <- i_F2 + which(diff(sign(OGF3(t_new[i_F2:len])))!=0)[1]
    F3 <- t_new[i_F3]
    ### PDA
    i_PDA <- i_P2 + which(diff(sign(W4(t_new[i_P2:len])))!=0)[1]
    PDA <- t_new[i_PDA]
    ### List of calculated A-Points
    a_points <- list(OGF3=F3, PDA=PDA)

  } else {
    inflection = FALSE
    f_points <- list(F1=NA, Fi=NA, F2=NA)
    p_points <- list(P1=NA, Pi=NA, P2=NA)
    l_points <- list(OGF0=NA, tang=NA, thres=NA)
    a_points <- list(OGF3=NA, PDA=NA)
  }
  
  list(
    f_points = f_points,
    p_points = p_points,
    l_points = l_points,
    a_points = a_points,
    inflection = inflection
  )
}