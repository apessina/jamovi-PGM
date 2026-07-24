
# functions to populate results tables

## Note helper
add_note <- function(key, ref, note, table) {
  
  if (!is.null(note)) {
    
    table$setNote(
      key  = key, 
      note = paste(ref, "-", note),
      init = FALSE
    )
    
  }
  
}


## Estimation info
populate_etable <- function(
  table,
  setup,
  label,
  est,
  opt,
  randomEff = NULL,
  weq,
  ceq,
  note = NULL
) {
  
  Li <- 0
  
  Li <- Li + 1
  table$addRow(
    rowKey  = paste0(setup, "L", Li),
    
    values  = list(
      var   = label,
      category = "Estimator",
      method = est
    )
    
  )
  
  Li <- Li + 1
  table$addRow(
    rowKey  = paste0(setup, "L", Li),
    
    values  = list(
      category = "Optimization",
      method = opt
    )
    
  )
  
  if (!is.null(randomEff)) {
    Li <- Li + 1
    table$addRow(
      rowKey  = paste0(setup, "L", Li),
      
      values  = list(
        category = "Random Effect",
        method = randomEff
      )
      
    )
  }
  
  Li <- Li + 1
  table$addRow(
    rowKey  = paste0(setup, "L", Li),
    
    values  = list(
      category = "Weight Function",
      method = weq
    )
    
  )
  
  Li <- Li + 1
  table$addRow(
    rowKey  = paste0(setup, "L", Li),
    
    values  = list(
      category = "Residual Correlation",
      method = ceq
    )
    
  )
  
  add_note(key = setup, ref = label, note, table)
}


## Parameters info
populate_ptable <- function(
  table,
  setup,
  label,
  pardf,
  note = NULL
) {
  
  if (!is.data.frame(pardf)) {
    stop("ERROR: pardf must be a dataframe in populate_ptable()")
  }
  
  for (i in seq_len(nrow(pardf))) {
    
    table$addRow(
      
      rowKey = paste0(setup, "L", i),
      
      values = list(
        var        = if (i==1) label,
        Parameter  = rownames(pardf)[i],
        Estimate   = pardf$parvalue[i],
        SE         = pardf$parerror[i],
        Lower      = pardf$parlower[i],
        Upper      = pardf$parupper[i]
      )
      
    )
    
  }
  
  add_note(key = setup, ref = label, note, table)
}


## Goodness of Fit metrics
populate_goftable <- function(
  table,
  setup,
  label,
  metrics,
  note = NULL
) {

  if (!is.list(metrics)) {
    stop("ERROR: metrics must be a list in populate_goftable()")
  }  
  
  table$addRow(
  
    rowKey = setup,
    
    values = list(
      var   = label,
      AIC   = metrics$AIC,
      AICc  = metrics$AICc,
      BIC   = metrics$BIC,
      r2    = metrics$R2,
      r2adj = metrics$R2adj
    )
  
  )
  
  add_note(key = setup, ref = label, note, table)
  
  table$setVisible(visible = TRUE)
}

## Error metrics
populate_emtable <- function(
  table,
  setup,
  label,
  metrics,
  note = NULL
) {
  
  if (!is.list(metrics)) {
    stop("ERROR: metrics must be a list in populate_emtable()")
  }
  
  table$addRow(
    
    rowKey = setup, 
    
    values = list(
      var   = label,
      RMSE  = metrics$RMSE,
      MAE   = metrics$MAE,
      MedAE = metrics$MedAE,
      sMAPE = metrics$sMAPE,
      RRMSE = metrics$RRMSE*100
    )
  
  )

  add_note(key = setup, ref = label, note, table)
  
  table$setVisible(visible = TRUE)
}

## Area Under the Curve metrics
populate_auctable <- function(
  table,
  setup,
  label,
  interval,
  metrics, 
  note = NULL
) {

  if (!is.list(metrics)) {
    stop("ERROR: metrics must be a list in populate_auctable()")
  }
    
  table$addRow(
      
    rowKey = setup, 
    
    values = list(
      var   = label,
      Lower = interval$aucLower,
      Upper = interval$aucUpper,
      gf    = metrics$aucGf,
      sgr   = metrics$aucSgr,
      gfSgr = metrics$aucGfSgr
    )
    
  )
  
  add_note(key = setup, ref = label, note, table)
  
  table$setVisible(visible = TRUE)
}

## Critical points location
populate_cptable <- function(
  table,
  setup,
  label,
  taxis, # name of time variable
  daxis, # name of dep variable
  critPoints, # list of critical points
  dfunc, # function of dep variable
  note = NULL
) {
  
  if (!is.list(critPoints)) {
    stop("ERROR: critPoints must be a list in populate_cptable()")
  }
  
  ## Location on time axis
  table$addRow(
    
    rowKey = paste0(setup, "_time"), 
    
    values = list(
      var       = label,
      Axis      = taxis, # op$time
      F0        = critPoints$lagPoints$OGF0,
      PAA       = critPoints$lagPoints$PAA,
      Tangent   = critPoints$lagPoints$tangent,
      Threshold = critPoints$lagPoints$threshold,
      P1        = critPoints$accPoints$P1,
      F1        = critPoints$ogfPoints$F1,
      iPoint    = critPoints$iPoint,
      F2        = critPoints$ogfPoints$F2,
      P2        = critPoints$accPoints$P2,
      F3        = critPoints$asyPoints$OGF3,
      PDA       = critPoints$asyPoints$PDA
    )
    
  )
  
  ## Location on dependent axis
  table$addRow(
    
    rowKey = paste0(setup, "_dep"), 
    
    values = list(
      Axis      = daxis,
      F0        = dfunc(critPoints$lagPoints$OGF0),
      PAA       = dfunc(critPoints$lagPoints$PAA),
      Tangent   = dfunc(critPoints$lagPoints$tangent),
      Threshold = dfunc(critPoints$lagPoints$threshold),
      P1        = dfunc(critPoints$accPoints$P1),
      F1        = dfunc(critPoints$ogfPoints$F1),
      iPoint    = dfunc(critPoints$iPoint),
      P2        = dfunc(critPoints$accPoints$P2),
      F2        = dfunc(critPoints$ogfPoints$F2),
      F3        = dfunc(critPoints$asyPoints$OGF3),
      PDA       = dfunc(critPoints$asyPoints$PDA)
    )
  
  )

  add_note(key = setup, ref = label, note, table)
  
  table$setVisible(visible = TRUE)
}


## Populate the Rank column of aucTable
rank_auctable <- function(
  table,
  aucs
) {
  
  aucDataFrame <- as.data.frame(
    do.call(
      rbind,
      lapply(aucs, function(x) {
        unlist(
          x$metrics[c(
            "aucGf",
            "aucSgr",
            "aucGfSgr"
          )]
        )
      })
    )
  )
  
  maxValues <- apply(
    aucDataFrame,
    2,
    max,
    na.rm = TRUE
  )
  maxValues[maxValues == 0] <- NA_real_
  
  rankDataFrame <- sweep(
    aucDataFrame,
    MARGIN = 2,
    STATS = maxValues,
    FUN = "/"
  ) * 100
  
  for (setupId in rownames(rankDataFrame)) {
    table$setRow(
      
      rowKey = setupId,
      
      values = list(
        gfRank    = rankDataFrame[setupId, "aucGf", drop = TRUE],
        sgrRank   = rankDataFrame[setupId, "aucSgr", drop = TRUE],
        gfSgrRank = rankDataFrame[setupId, "aucGfSgr", drop = TRUE]
      )
      
    )
  }
  
}