
# additional functions and helpers for models specifications in specsmodel.R

## switch options name by model spec function
get_model_spec <- function(model, t, y) {
  
  switch(
    model,
    richards = spec_richards(t, y),
    stop("ERROR: Invalid model equation")
  )
  
}


## populate HTML block in results with model name and equation 
populate_eqhtml <- function(results, equations) {
  
  if (!is.list(equations)) {
    stop("ERROR: equations must be a list in populate_eqhtml()")
  }
  
  
  blocks <- character()
  
  for (mname in names(equations)) {
    
    blocks <- c(blocks, 
      paste0(
        '<p style="font-family:Segoe UI, sans-serif; ',
        'font-size:14;"><ul><li>',mname,': ',
        equations[[mname]],'</li></ul>',
        '</p></div>'
      )
    )
    
  }
  
  equationsHtml <- paste0(blocks, collapse = "\n")
  
  
  results$eqHtml$setContent(
    paste0(
      '<div style="margin-top:6px;">',
      
      '<div style="font-family:Segoe UI, sans-serif; ',
      'font-size:12px; margin-bottom:3px;">',
      'Models Equations</div>',
      
      equationsHtml
    )
  )
  
}

