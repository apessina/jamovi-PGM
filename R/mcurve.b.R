
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
        Ke_df <- data.frame(t=t_raw_new)
        action_val <- list()
        
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
          
          ## Expressions to functions
          W  <- function(x) eval(W_expr) # W(t)
          W1 <- function(x) eval(W1_expr) # W'(t)
          
          ## Define "Kinetic Energy" function
          Ke <- function(x) 0.5 * W(x) * W1(x)^2
          
          ##### Evaluation Metrics #####
          
          ## Number of samples
          n = length(y)
          
          ## Degrees of freedom
          df1 <- n_params - 1 # model
          df2 <- n - n_params # error
          
          ## Sum of squares
          SSt <- sum((y - mean(y)) ^ 2) # total
          SSm <- sum((W(t) - mean(y)) ^ 2) # model
          SSe <- sum((y - W(t)) ^ 2) # error
          
          ## Goodness-of-fit metrics
          R2 <- 1 - SSe / SSt # R²
          R2_adj <- 1 - (SSe/df2) / (SSt/(n - 1)) # Adjusted R²
          AIC <- n * log(SSe/n) + 2 * n_params
          AICc <- AIC + (2 * n_params * (n_params + 1)) / (n-n_params - 1)
          BIC <- n * log(SSe/n) + log(n) * n_params
          
          ### Global F-test
          F_s <- ((SSt - SSe) / df1) / (SSe / df2) # F statistics
          F_p <- pf(F_s, df1, df2, lower.tail=FALSE) # F p-value
          
          ## Residuals
          res <- y - W(t)
          MSE  <- mean(res^2)
          
          ## Error metrics
          RMSE <- sqrt(MSE)
          MAE <- mean(abs(res))
          MedAE <- median(abs(res))
          sMAPE <- mean(2 * abs(res) / (abs(y) + abs(W(t)))) * 100
          RRMSE <- RMSE / mean(y)
          
          ##### Curve "resolution" #####
          
          ## Trim t_raw_new unit the max value of t.
          ## This allows graphical comparison in the same scale,
          ## even if measures stopped at a certain age for one variable
          t_new <- t_raw_new[t_raw_new <= max(t)]
          
          ## Apply Data Modeling functions to t_new
          W_pred <- W(t_new) # vector of weight over time
          Ke_pred <- Ke(t_new) # vector of Ke over time
          
          ## Calculate the "Action"
          if (self$options$intInterval) {
            intKe_lower <- self$options$intL
            intKe_upper <- self$options$intU
          } else {
            intKe_lower <- min(t_new, na.rm=TRUE)
            intKe_upper <- max(t_new, na.rm=TRUE)
          }
          intKe_pred <- integrate(Ke, lower=intKe_lower, upper=intKe_upper)$value
          
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
            R2_adj=format(round(R2_adj, 3), nsmall=3),
            RMSE=format(round(RMSE, 3), nsmall=3),
            MAE=format(round(MAE, 3), nsmall=3),
            MedAE=format(round(MedAE, 3), nsmall=3),
            sMAPE=sprintf("%.3f%%", round(sMAPE, 3)),
            RRMSE=sprintf("%.3f%%", round(RRMSE*100, 3)),
            fdf1=round(df1),
            fdf2=round(df2),
            f=format(round(F_s, 3), nsmall=3),
            fp=F_p
          ))
          
          ## Fill the "Action" table
          tableKe <- self$results$action
          tableKe$setRow(rowKey=dep, values=list(
            Lower=format(round(intKe_lower, 2), nsmall=2),
            Upper=format(round(intKe_upper, 2), nsmall=2),
            Action=intKe_pred
          ))
          action_val[dep] <- intKe_pred

          ## Fill df length with NA to plot trimmed dep lengths
          ## Save values for plots
          t_df[[dep]] <- c(y_raw, rep(NA, nrow(t_df)-length(y_raw)))
          t_new_df[[dep]] <- c(W_pred, rep(NA, nrow(t_new_df)-length(W_pred)))
          Ke_df[[dep]] <- c(Ke_pred, rep(NA, nrow(t_new_df)-length(Ke_pred)))
          
        } # close loop
        
        ## Rank Action table
        action_num <- unlist(action_val)
        rank_val <- (action_num / max(action_num)) * 100
        for (dep in deps) {
          tableKe$setRow(rowKey=dep, values=list(
            Rank=sprintf("%.2f%%", rank_val[dep])
          ))
        }
        
        ##### Plots Data #####
        
        ## Data for next functions
        private$.prep_mplot(t_df, t_new_df)
        private$.prep_aplot(Ke_df)
        
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
        colors_pred <- jmvcore::colorPalette(n=length(y_pred_cols))
        
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
          lines(t_m, data_m[[col_name]], col=colors_pred[i], lwd=2)
        }
        
        ## Make legend
        legend("bottom", inset=-0.2, legend=y_pred_cols, col=colors_pred, 
               lty=1, lwd=2, horiz=TRUE, bty="n", xpd=TRUE)
        
        ## Notify the rendering system that we have plotted something
        TRUE
        
      }, # close .mplot
      
      # Prepare the data for "Action" plot
      .prep_aplot = function(Ke_df) {
        ## Set-up the plot with the model data in a dataframe
        image <- self$results$aplot
        image$setState(list(
          data_Ke=Ke_df
        )
        )
      }, # close .prep_aplot
      
      # Create the derivative plot function
      .aplot = function(image, ...) {
        
        ## Check if there is data for the plot
        if (is.null(image$state))
          return(FALSE)
        
        data_Ke <- image$state$data_Ke
        
        ## Time columns and y columns names
        t <- data_Ke$t
        Ke_cols <- setdiff(names(data_Ke), "t") # all columns, except t
        
        ## Generate line colors in the standard jamovi palette
        colors_pred <- jmvcore::colorPalette(n=length(Ke_cols))
        
        ## Deactivate the default box around the plotting area
        par(bty='L')
        
        ## Better limits to avoid points without axis labels
        better_lim <- function(vec, n=5) {
          ticks <- pretty(vec, n = n)
          c(min(ticks), max(ticks))
        }
        
        ## Join all y data to calculate limits
        ylim_total <- range(unlist(data_Ke[Ke_cols]), na.rm=TRUE)
        ylim_total <- better_lim(ylim_total)
        
        ## Create and set-up the plot with 1st curve
        plot(t, data_Ke[[Ke_cols[1]]], type="l", col=colors_pred[1],
             lty=1, xlab="", ylab="", ylim=ylim_total)
        
        ## Add other predicted curves
        if (length(Ke_cols) > 1) {
          for (i in 2:length(Ke_cols)) {
            col_name <- Ke_cols[i]
            lines(t, data_Ke[[col_name]], col=colors_pred[i], lty=1)
          }
        }
        
        ## Make legend
        legend("bottom", inset=-0.2, legend=Ke_cols, col=colors_pred, lty=1,
               horiz=TRUE, bty="n", lwd=2, xpd=TRUE)
        
        ## Notify the rendering system that we have plotted something
        TRUE
        
      } # close .aplot
      
    ) # close list
) # close R6Class

