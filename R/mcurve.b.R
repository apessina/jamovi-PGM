
# This file is a generated template, your changes will not be overwritten

mcurveClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "mcurveClass",
    inherit = mcurveBase,
    private = list(
      
      ##### Main Function #####
      
      .run = function() {
        
        ##### Get User Data #####
        
        ## Option values into shorter variable names
        deps  <- self$options$deps
        time <- self$options$time
        
        ## Check if variables have any data
        if (is.null(deps) || is.null(time))
          return()
        
        ## Get the data
        data <- self$data
        
        ## Convert time to appropriate data type
        data[[time]] <- jmvcore::toNumeric(data[[time]])
        
        ##### Models Functions #####
        
        ## Richards Model
        richards <- function(t, y) {
          
          ### Define model expression
          s_expr <- expression(
            A * ( 1 + (d - 1) * exp(-K * (t - Ti) / (d / (1 - d))) ) ^ (1 / (1 - d))
          )
          
          ### Initial parameters values
          gr <- diff(y)/diff(t) # growth rate vector (1st derivative)
          A_init <- max(y, na.rm=TRUE) # asymptotic - maximum y value
          K_init <- max(gr) / A_init # relative max. growth rate
          Ti_init <- t[which.max(gr)] # inflection point
          d_inits <- seq(0.5, 3, by=0.01) # fixed "d" values for grid search
          
          ### Parameters restriction
          if (self$options$pConstraint=="strict") {
            d_inits <- d_inits[d_inits<0.95|d_inits>1.05] # remove d~1
          } else {
            d_inits <- d_inits[d_inits!=1] # remove d=1
          }
          
          ### Store best iteration results
          best_model <- NULL
          best_SSe <- Inf
          best_d <- NA
          
          ### Iteration to estimate best parameters for each fixed "d" value
          for (d_init in d_inits) {
            try({
              fit <- nlsLM(
                y ~ A * (1 + (d_init - 1) * exp(-K * (t - Ti) / (d_init / (1 - d_init)))) ^ (1 / (1 - d_init)),
                start = list(A=A_init, K=K_init, Ti=Ti_init),
                lower = c(A=0, K=0, Ti=min(t)), upper = c(Inf, Inf, max(t)),
                control = nls.lm.control(maxiter=1000, ftol=1e-10)
              )
              SSe <- sum(residuals(fit)^2)
              if (SSe < best_SSe) {
                best_model <- fit
                best_SSe <- SSe
                best_d <- d_init
              }
            }, silent=TRUE)
          }
          
          ### Get final parameters in the expression
          params <- as.list(coef(best_model)) # estimated A, K, Ti
          params$d <- best_d # best d value on grid search
          f_expr <- eval(substitute(expression(
            A * (1 + (d - 1) * exp(-K * (x - Ti) / (d / (1 - d)))) ^ (1 / (1 - d))
          ), params))
          
          return(list(s_expr=s_expr, params=params, f_expr=f_expr))
        }
        
        ##### Curve "resolution" #####
        
        ## Giving time a new "resolution" for smooth lines
        t_raw <- data[[time]]
        t_raw_res <- length(t_raw) * self$options$res # res is the "resolution factor"
        t_raw_new <- seq(min(t_raw), max(t_raw), length.out=t_raw_res) # same range, more points
        
        ##### Iteration p/ dep #####
        t_df <- data.frame(t=t_raw)
        t_new_df <- data.frame(t=t_raw_new)
        
        for (dep in deps) {
        
          ## Convert dep to appropriate data type
          data[[dep]] <- jmvcore::toNumeric(data[[dep]])
          
          ## Remove NAs
          data_cln <- na.omit(data[, c(time, dep)])
          
          ## Define time and dep variables
          x_raw <- data_cln[[time]]
          y_raw <- data_cln[[dep]]
          
          ## Check if both variables have at least two unique values
          if (length(unique(x_raw)) < 2 || length(unique(y_raw)) < 2)
            stop("At least two unique (x, y) pairs are required for growth analysis")
          
          ## Aggregate multiple y-values with same x-value
          trim_perc <- ifelse(self$options$trim, self$options$tPerc/100, 0)
          if (self$options$agg == "median") {
            agg_vals <- aggregate(y_raw, by=list(x_raw), FUN=median)
          } else {
            agg_vals <- aggregate(y_raw, by=list(x_raw), FUN=function(x) 
              mean(x, trim=trim_perc))
          }
          t <- agg_vals[[1]]
          y <- agg_vals[[2]]
          
          ##### Data Modeling #####
          
          ## Adjust selected model to user data
          if (self$options$model == "richards") {
            model <- richards(t, y)
          }
          
          ## Get model information
          s_expr <- model$s_expr # raw expression
          params <- model$params # values of parameters
          n_params <- length(params) # number of parameters
          W_expr <- model$f_expr # final/fit expression
          
          ## First 2 symbolic derivatives as expression
          W1_expr <- D(W_expr, "x") # 1st derivative
          W2_expr <- D(W1_expr, "x") # 2nd derivative
          
          ## Ontogenetic Growth Force expression
          OGF_expr <- call("*", W1_expr, W2_expr) # W'(t) * W''(t)
          
          ### First 3 symbolic derivatives of OGF expressions
          OGF1_expr <- D(OGF_expr, "x") # 1st OGF derivative
          OGF2_expr <- D(OGF1_expr, "x") # 2nd OGF derivative 
          OGF3_expr <- D(OGF2_expr, "x") # 3rd OGF derivative 
          
          ## Expressions to functions
          W  <- function(x) eval(W_expr) # W(t)
          W1 <- function(x) eval(W1_expr) # W'(t)
          W2 <- function(x) eval(W2_expr) # W''(t)
          OGF <- function(x) eval(OGF_expr) # OGF(t)
          OGF3 <- function(x) eval(OGF3_expr) # OGF'''(t)
          
          ##### Goodness-Of-Fit #####
          
          ## Number of samples
          n = length(y)
          
          ## Degrees of freedom
          df1 <- n_params - 1
          df2 <- n - n_params
          
          ## Sum of squares
          SSt <- sum((y - mean(y)) ^ 2) # total
          SSm <- sum((W(t) - mean(y)) ^ 2) # model
          SSe <- sum((y - W(t)) ^ 2) # error
          
          ## Goodness-of-fit metrics
          R2 <- 1 - SSe / SSt # R²
          AIC <- n * log(SSe/n) + 2 * n_params
          AICc <- AIC + (2 * n_params * (n_params + 1)) / (n-n_params - 1)
          BIC <- n * log(SSe/n) + log(n) * n_params
          
          ### Global F-test
          F_s <- ((SSt - SSe) / df1) / (SSe / df2) # F statistics
          F_p <- pf(F_s, df1, df2, lower.tail=FALSE) # F p-value
          
          ##### Curve "resolution" #####
          
          ## Trim t_raw_new unit the max value of t.
          ## This allows graphical comparison in the same scale,
          ## even if measures stopped at a certain age for one variable
          t_new <- t_raw_new[t_raw_new <= max(t)]
          
          ## Apply Data Modeling functions to t_new
          W_pred <- W(t_new) # vector of weight over time
          
          ##### Model Information #####
          
          ## Show model equation as text
          str_expr <- paste(deparse(s_expr[[1]]), collapse="")
          eq_str <- paste0("\nModel equation:\n", str_expr)
          self$results$text$setContent(eq_str)
          
          ## Estimated parameters table
          pTable <- self$results$pTable
          pTable$setRow(rowKey=dep, values=list(
            A=format(round(params$A, 2), nsmall=2),
            d=format(round(params$d, 2), nsmall=2),
            K=format(round(params$K, 2), nsmall=2),
            Ti=format(round(params$Ti, 2), nsmall=2)
          ))
          
          ## Goodness-Of-Fit table
          tableFit <- self$results$fitq
          tableFit$setRow(rowKey=dep, values=list(
            AIC=format(round(AIC, 2), nsmall=2),
            AICc=format(round(AICc, 2), nsmall=2),
            BIC=format(round(BIC, 2), nsmall=2),
            R2=format(round(R2, 3), nsmall=3),
            fdf1=round(df1),
            fdf2=round(df2),
            f=format(round(F_s, 3), nsmall=3),
            fp=F_p
          ))

          ## Fill df length with NA to plot different dep lengths...
          ## ...and save values for plots
          t_df[[dep]] <- c(y_raw, rep(NA, nrow(t_df)-length(y_raw)))
          t_new_df[[dep]] <- c(W_pred, rep(NA, nrow(t_new_df)-length(W_pred)))
          # ge_df[[dep]] <- ge_pred
          
        }          
        
        ##### Plots Data #####
        
        ## Data for next functions
        private$.prep_mplot(t_df, t_new_df)
        
      }, # close .run
      
      ##### Plots Functions #####
      
      # Prepare data for Model Plot
      .prep_mplot = function(t_df, t_new_df) {
        ## Set plot with model data as dataframe
        image <- self$results$mplot
        image$setState(list(
          data_p=t_df,
          data_m=t_new_df
        )
        )
      }, # close .prep_mplot
      
      # Model Plot Function
      .mplot = function(image, ...) {
        
        ## Check if there is any data to plot
        if (is.null(image$state))
          return(FALSE)
        
        data_p <- image$state$data_p
        data_m <- image$state$data_m
        
        ## Time columns and y columns names
        t_p <- data_p$t
        t_m <- data_m$t
        y_cols <- setdiff(names(data_p), "t")
        y_pred_cols <- setdiff(names(data_m), "t")
        
        ## Generate line colors in the standard jamovi palette
        cores_pred <- jmvcore::colorPalette(n=length(y_pred_cols))
        
        ## Better limits to avoid points without axis labels
        better_lim <- function(vec, n=5) {
          ticks <- pretty(vec, n=n)
          c(min(ticks), max(ticks))
        }
        
        ## Join all y data to calculate limits
        ylim_total <- range(unlist(data_p[y_cols]), unlist(data_m[y_pred_cols]), na.rm = TRUE)
        ylim_total <- better_lim(ylim_total)
        
        ## Deactivate the default box around the plotting area
        par(bty='L')
        
        ## Create and set-up the plot
        plot(t_p, data_p[[y_cols[1]]], col="grey", pch=16, xlab="", ylab="",
             ylim = ylim_total)
        
        ## Add the other variables points
        if (length(y_cols) > 1) {
          for (col in y_cols[-1]) {
            points(t_p, data_p[[col]], col="grey", pch=16)
          }
        }
        
        ## Add predicted curves
        for (i in seq_along(y_pred_cols)) {
          col_name <- y_pred_cols[i]
          lines(t_m, data_m[[col_name]], col=cores_pred[i], lwd=2)
        }
        
        ## Make legend
        legend("bottom", inset=-0.2, legend=y_pred_cols, col=cores_pred, 
               lty=1, lwd=2, horiz=TRUE, bty="n", xpd=TRUE)
        
        ## Notify the rendering system that we have plotted something
        TRUE
        
      } # close .mplot
      
    ) # close list
) # close R6Class

