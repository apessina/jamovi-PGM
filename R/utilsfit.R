
# 


## 
pgm_fixed_start <- function(
    formulas,
    data,
    baseStart
) {
  
  start <- unlist(
    lapply(names(formulas), function(p) {
      
      X <- model.matrix(
        delete.response(terms(formulas[[p]])), data = data
      )
      
      values <- setNames(numeric(ncol(X)), colnames(X))
      
      if ("(Intercept)" %in% names(values))
        values["(Intercept)"] <- baseStart[[p]]
      
      if (length(values) == 1L && identical(names(values), "(Intercept)")) {
        names(values) <- p
      } else {
        names(values) <- paste(p, names(values), sep = ".")
      }
      
      values
    }),
    use.names = TRUE
  )
  
  if (anyNA(start) || any(!is.finite(start))) {
    stop(
      "Non-finite initial values were generated for the fixed effects.",
      call. = FALSE
    )
  }
  
  start
}


## recognize data nesting hierarchy
pgm_isnested <- function(data, inner, outer) {
  
  groups <- split(data[[outer]], data[[inner]], drop = TRUE)
  
  all(vapply(groups, function(x) length(unique(x)) == 1L, logical(1)))
}


## 
pgm_group_hierarchy <- function(
    data,
    randomGroups,
    corrGroups = character(0),
    groupComponents = list()
) {
  
  randomGroups <- lapply(randomGroups, unique)
  
  corrGroups <- unique(corrGroups)
  
  randomVars <- unique(unlist(randomGroups, use.names = FALSE))
  
  allVars <- unique(c(randomVars, corrGroups))
  
  if (!length(allVars))
    return(list(random = character(0), correlation = character(0)))
  
  if (!all(allVars %in% names(data))) {
    
    missingVars <- setdiff(allVars, names(data))
    
    stop(
      "Grouping variables not found in fitting data: ",
      paste(missingVars, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  
  
  ## Components represented by each grouping variable
  components_of <- function(g) {
    if (length(groupComponents) && g %in% names(groupComponents))
      return(groupComponents[[g]])
    g
  }
  
  if (length(allVars) == 1L) {
    
    return(list(
      random = randomVars,
      correlation =
        if (length(corrGroups))
          allVars
      else
        character(0)
    ))
  }
  
  ## Pairwise ordering:
  ## row must precede column
  before <- matrix(
    FALSE,
    nrow = length(allVars),
    ncol = length(allVars),
    dimnames = list(allVars, allVars)
  )
  
  
  ## Explicit correlation hierarchy supplied by user
  add_order <- function(x) {
    
    x <- unique(x)
    
    if (length(x) < 2L)
      return()
    
    for (i in seq_len(length(x) - 1L)) {
      for (j in seq.int(i + 1L, length(x))) {
        before[x[[i]], x[[j]]] <<- TRUE
      }
    }
  }
  
  add_order(corrGroups)
  
  ## Infer additional nesting from data
  pairs <- utils::combn(allVars, 2L)
  
  for (i in seq_len(ncol(pairs))) {
    
    a <- pairs[1L, i]
    b <- pairs[2L, i]
    
    if (before[a, b] || before[b, a])
      next
    
    aInB <- pgm_isnested(data, a, b)
    bInA <- pgm_isnested(data, b, a)
    
    if (aInB && !bInA) {
      before[b, a] <- TRUE
    } else if (bInA && !aInB) {
      before[a, b] <- TRUE
    }
  }
  
  ## Unique outer -> inner ordering
  remaining <- allVars
  hierarchy <- character(0)
  
  while (length(remaining)) {
    
    roots <- remaining[
      vapply(
        remaining,
        function(x) !any(before[remaining, x]),
        logical(1)
      )
    ]
    
    if (!length(roots))
      stop("Conflicting grouping hierarchy.", call. = FALSE)
    
    if (length(roots) > 1L) {
      stop(
        "The selected grouping variables do not form a unique nested hierarchy. ",
        "The following variables appear to be crossed or otherwise unordered: ",
        paste(
          roots,
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }
    
    hierarchy <- c(hierarchy, roots)
    
    remaining <- setdiff(remaining, roots)
  }
  
  randomOrder <- hierarchy[hierarchy %in% randomVars]
  
  if (!length(corrGroups))
    return(list(random = randomOrder, correlation = character(0)))
  
  ## Complete hierarchy down to deepest requested
  ## correlation grouping variable
  deepestCor <- max(match(corrGroups, hierarchy))
  
  baseCorOrder <- hierarchy[seq_len(deepestCor)]
  
  ## No random effects:
  ## ordinary correlation hierarchy
  if (!length(randomOrder))
    return(list(random = character(0), correlation = baseCorOrder))
  
  deepestRandom <- max(match(randomOrder, hierarchy))
  
  ## Atomic variables represented by the random hierarchy.
  ## A composite group such as Species:Replicate covers
  ## both Species and Replicate.
  randomComponents <- unique(
    unlist(lapply(randomOrder, components_of), use.names = FALSE)
  )
  
  ## Correlation levels occurring above or within the
  ## deepest random level must either be actual random
  ## groups or components represented by a composite group.
  beforeRandom <- baseCorOrder[match(baseCorOrder, hierarchy) <= deepestRandom]
  
  unresolved <- setdiff(setdiff(beforeRandom, randomVars), randomComponents)
  
  if (length(unresolved)) {
    stop(
      "The correlation and random-effects grouping ",
      "hierarchies are incompatible. Correlation grouping ",
      "variables not represented by the random-effects ",
      "hierarchy: ",
      paste(unresolved, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  
  ## Correlation must use the random grouping hierarchy
  ## at least as deeply as the random effects.
  corOrder <- randomOrder
  
  ## If correlation extends deeper than the random
  ## hierarchy, append those additional correlation levels.
  if (deepestCor > deepestRandom) {
    
    deeper <- hierarchy[seq.int(deepestRandom + 1L, deepestCor)]
    
    deeper <- deeper[!deeper %in% randomVars]
    
    corOrder <- unique(c(corOrder, deeper))
  }
  
  
  list(
    random = randomOrder,
    correlation = corOrder
  )
}


## Normalize jamovi Terms while preserving interactions
pgm_normalize_terms <- function(terms) {
  
  if (is.null(terms) || !length(terms))
    return(list())
  
  ## Backward compatibility with flat main-effect vectors
  if (is.character(terms))
    return(lapply(terms, as.character))
  
  lapply(
    terms,
    function(term) as.character(unlist(term, use.names = FALSE))
  )
}


## Convert jamovi Terms to R formula terms
pgm_compose_terms <- function(terms) {
  
  terms <- pgm_normalize_terms(terms)
  
  if (!length(terms))
    return(character(0))
  
  jmvcore::composeTerms(terms)
}


## Variables represented by a jamovi Terms object
pgm_term_variables <- function(terms) {
  
  unique(unlist(pgm_normalize_terms(terms), use.names = FALSE))
}


## Fixed-effects formula for one structural parameter
pgm_fixed_formula <- function(
    parameter,
    terms
) {
  
  stats::reformulate(pgm_compose_terms(terms), response = parameter)
}


##
pgm_prepare_random_groups <- function(
    data,
    randomTerms
) {
  
  allTerms <- do.call(c, unname(randomTerms))
  
  if (is.null(allTerms) || !length(allTerms)) {
    return(list(
      data = data,
      groups = randomTerms,
      labels = character(0),
      components = list()
    ))
  }
  
  term_key <- function(term) paste(term, collapse = "\034")
  
  ## Remove duplicated grouping terms across parameters
  keys <- vapply(allTerms, term_key, character(1))
  
  keep <- !duplicated(keys)
  
  allTerms <- allTerms[keep]
  keys <- keys[keep]
  
  internal <- character(length(allTerms))
  
  labels <- character(length(allTerms))
  
  
  for (i in seq_along(allTerms)) {
    
    term <- allTerms[[i]]
    
    labels[[i]] <- paste(term, collapse = ":")
    
    ## Simple grouping variable
    if (length(term) == 1L) {
      
      internal[[i]] <- term[[1]]
      
    } else {
      
      ## Internal variable for an interaction grouping factor
      candidate <- paste0(".pgmGroup", i)
      
      j <- i
      
      while (candidate %in% names(data)) {
        j <- j + 1L
        candidate <- paste0(".pgmGroup", j)
      }
      
      internal[[i]] <- candidate
      
      data[[candidate]] <- do.call(
        interaction,
        c(
          unname(data[term]),
          list(
            drop = TRUE,
            lex.order = TRUE,
            sep = ":"
          )
        )
      )
    }
  }
  
  keyToInternal <- setNames(internal, keys)
  
  groups <- lapply(
    randomTerms,
    function(terms) {
      
      if (!length(terms))
        return(character(0))
      
      vapply(
        terms,
        function(term) keyToInternal[[term_key(term)]],
        character(1)
      )
    }
  )
  
  
  names(groups) <- names(randomTerms)
  
  
  list(
    data = data,
    groups = groups,
    labels = setNames(labels, internal),
    components = setNames(allTerms, internal)
  )
}


##
pgm_group_label <- function(
    group,
    labels
) {
  
  if (length(labels) && group %in% names(labels))
    return(unname(labels[[group]]))
  
  group
}

pgm_validate_fixed_terms <- function(
    data,
    fixedTerms
) {
  
  problems <- character(0)
  
  for (p in names(fixedTerms)) {
    terms <- fixedTerms[[p]]
    if (!length(terms))
      next
    vars <- unique(unlist( terms, use.names = FALSE)
    )
    
    
    for (v in vars) {
      
      if (!v %in% names(data)) {
        stop(
          "Fixed-effect variable not found in fitting data: ",
          v,
          ".",
          call. = FALSE
        )
      }
      
      x <- data[[v]]
      
      ## Factors require at least two observed levels
      if (is.factor(x) || is.character(x) || is.logical(x)) {
        observed <- unique(x[!is.na(x)])
        if (length(observed) < 2L)
          problems <- c(problems, paste0(v, " (parameter ", p, ")"))
      }
    }
  }
  
  
  if (length(problems)) {
    
    stop(
      paste0(
        "Fixed-effect factor(s) have fewer than two observed levels ",
        "in the data used for this fit: ",
        paste(unique(problems), collapse = ", "),
        ". The corresponding main effect or interaction cannot be estimated."
      ),
      call. = FALSE
    )
  }
  
  
  invisible(TRUE)
}


##
pgm_validate_fixed_design <- function(
    formulas,
    data
) {
  
  for (p in names(formulas)) {
    
    X <- stats::model.matrix(
      stats::delete.response(stats::terms(formulas[[p]])), data = data
    )
    
    if (anyNA(X) || any(!is.finite(X))) {
      stop(
        "The fixed-effects design matrix for parameter ",
        p,
        " contains missing or non-finite values.",
        call. = FALSE
      )
    }
    
    qrx <- qr(X)
    
    if (qrx$rank < ncol(X)) {
      
      dependent <- colnames(X)[qrx$pivot[seq.int(qrx$rank + 1L, ncol(X))]]
      
      stop(
        paste0(
          "The fixed-effects design for parameter ",
          p,
          " is rank deficient. Non-estimable coefficient(s): ",
          paste(dependent, collapse = ", "),
          ". Check for empty combinations of factor levels, ",
          "redundant terms, or perfectly collinear predictors."
        ),
        call. = FALSE
      )
    }
  }
  
  invisible(TRUE)
}

pgm_parameter_terms <- function(x) {
  
  lapply(x, function(item) {
    
    if (
      is.list(item) &&
      !is.null(names(item)) &&
      "terms" %in% names(item)
    ) {
      return(item$terms)
    }
    
    item
  })
}