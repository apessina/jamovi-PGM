
# This file is a generated template, your changes will not be overwritten

scurveClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "scurveClass",
    inherit = scurveBase,
    private = list(
      
      ##### Main Function #####
      
      .run = function() {
        
        ##### Data -----
        
        ## Option values into shorter variable names
        dep  <- self$options$dep
        time <- self$options$time
        
        ## Check if variables have any data
        if (is.null(dep) || is.null(time))
          return()
        
        ## Get the data
        data <- self$data[, c(dep, time), drop=FALSE]
        
        ## Convert to appropriate data types
        data[[dep]] <- jmvcore::toNumeric(data[[dep]])
        data[[time]] <- jmvcore::toNumeric(data[[time]])
        
        ## Remove NAs
        data <- na.omit(data)
        
        ## Convert to appropriate data types
        y_raw <- data[[dep]]
        x_raw <- data[[time]]
        
        ## Check if both variables have at least two unique values
        if (length(unique(x_raw)) < 2 || length(unique(y_raw)) < 2)
          stop("At least two unique (x, y) pairs are required for growth analysis")
        
        ##### Aggregation -----
        
        ## Aggregate multiple y-values with same x-value
        agg_vals <- aggregate(y_raw, by=list(x_raw), FUN=mean)
        t <- agg_vals[[1]]
        y <- agg_vals[[2]]
        
        ##### Models -----
        
        ## Some definitions for parameters domains
        gr <- diff(y)/diff(t) # growth rate vector (1st derivative)
        a_sup <- max(y, na.rm=TRUE)
        a_inf <- min(y, na.rm=TRUE)
        
        if (self$options$model == "richards") {
          expr <- expression(
            A * (1 + (d - 1) * exp(-K * (t - Ti) / (d^(d / (1 - d)))))^(1 / (1 - d))
          )
          init <- list(A=a_sup, K=max(gr)/a_sup, Ti=t[which.max(gr)], d=1.2)
          lower <- c(A=0, K=0, Ti=min(t), d=1.05)
          upper <- c(Inf, Inf, max(t), Inf)
        }
      
        
        ## Show model equation as text
        eq_str <- paste0("Model:\n", paste(expr[[1]], collapse=""))
        self$results$text$setContent(eq_str)

        ##### Modeling -----
        
        ## Iteration to estimate best parameters for each fixed "d" value
        fit <- nlsLM(
          formula=as.formula(call("~", quote(y), expr[[1]]), env=.GlobalEnv), 
          data=data.frame(y=y ,t=t), start=init, lower=lower, upper=upper,
          control = nls.lm.control(maxiter=1000, ftol=1e-10)
        )
        
        ## Get final parameters in the expression
        params <- coef(fit) # estimated parameters
        
        ## Model equation with fit parameters as expression
        W_expr <- as.expression(
          do.call('substitute', list(expr[[1]], as.list(params)))
        )
        
        ## First 2 symbolic derivatives as expression
        W1_expr <- D(W_expr, "t") # 1st derivative
        W2_expr <- D(W1_expr, "t") # 2nd derivative
        
        ## 3rd and 4th derivatives for PDA
        W3_expr <- D(W2_expr, "t") # 3rd derivative
        W4_expr <- D(W3_expr, "t") # 4th derivative
        
        ## Ontogenetic Growth Force expression
        OGF_expr <- call("*", W1_expr, W2_expr) # W'(t) * W''(t)
        
        ### First 3 symbolic derivatives of OGF expressions
        OGF1_expr <- D(OGF_expr, "t") # 1st OGF derivative
        OGF2_expr <- D(OGF1_expr, "t") # 2nd OGF derivative 
        OGF3_expr <- D(OGF2_expr, "t") # 3rd OGF derivative 
        
        ## Expressions to functions
        W  <- function(t) eval(W_expr, envir=list(t=t)) # W(t) 
        W1 <- function(t) eval(W1_expr, envir=list(t=t)) # W'(t)
        W2 <- function(t) eval(W2_expr, envir=list(t=t)) # W''(t)
        W4 <- function(t) eval(W4_expr, envir=list(t=t)) # W''''(t) (PDA)
        OGF <- function(t) eval(OGF_expr, envir=list(t=t)) # OGF(t)
        OGF3 <- function(t) eval(OGF3_expr, envir=list(t=t)) # OGF'''(t)
        
        ##### Parameters -----
        
        ## Number of samples
        n = length(y)
        
        ## Number of parameters
        n_params <- length(params)

        ## Degrees of freedom
        df1 <- n_params - 1 # model
        df2 <- n - n_params # error
        
        ## Covariance Matrix
        covm <- vcov(fit)
        
        ## Standard Errors 
        se_vals <- sqrt(diag(covm))
        
        ## t Values - Wald t-test
        t_vals <- params / se_vals
        
        ## p-values - Wald t-test
        p_vals <- 2 * pt(abs(t_vals), df=df2, lower.tail=FALSE)
        
        ## Confidence Intervals - Wald t-test
        alpha <- 0.05 # fixed in 95% for now
        t_crit <- qt(1 - alpha / 2, df=df2)
        lower_ci <- params - t_crit * se_vals
        upper_ci <- params + t_crit * se_vals
        
        ## Parameters DF
        params_df <- data.frame(
          param=names(params),
          estim=params,
          lower=lower_ci,
          upper=upper_ci,
          se=se_vals,
          t=t_vals,
          p=p_vals
        )
        
        ##### Evaluation -----
        
        ## Residuals / Error
        res <- residuals(fit)
        
        ## Sum of squares
        SSt <- sum((y - mean(y)) ^ 2) # total
        SSe <- sum(res ^ 2) # error
        
        ## Goodness-of-fit metrics
        R2 <- 1 - SSe / SSt # R²
        R2_adj <- 1 - (SSe/df2) / (SSt/(n - 1)) # Adjusted R²
        AIC <- n * log(SSe/n) + 2 * n_params
        AICc <- AIC + (2 * n_params * (n_params + 1)) / (n-n_params - 1)
        BIC <- n * log(SSe/n) + log(n) * n_params
        
        ## Error metrics
        MSE  <- mean(res ^ 2)
        RMSE <- sqrt(MSE)
        MAE <- mean(abs(res))
        MedAE <- median(abs(res))
        sMAPE <- mean(2 * abs(res) / (abs(y) + abs(fitted(fit)))) * 100
        RRMSE <- RMSE / mean(y)
        
        ##### Curve Resolution -----
        
        ## Giving t a new "resolution" for smooth lines
        t_res <- length(t) * self$options$res # res is the "resolution factor"
        t_new <- seq(min(t), max(t), length.out=t_res) # same range, more points
        
        ## Apply Data Modeling functions to t_new
        W_pred <- W(t_new) # vector of weight over time
        W1_pred <- W1(t_new) # vector of growth rate over time
        W2_pred <- W2(t_new) # vector of acceleration over time
        W4_pred <- W4(t_new) # vector of 4th derivative over time (PDA)
        OGF_pred <- OGF(t_new) # vector of OGF over time
        OGF3_pred <- OGF3(t_new) # vector of the OGF 3rd Derivative over time
        
        ##### Critical Points -----
        
        ## Length of x-axis
        len <- length(t_new)
        
        ## Check if the data have an inflection point
        zero_acc <- which(diff(sign(W2_pred))!=0)
        if (!length(zero_acc)==0) {
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
          y_thresh <- W(max(t_new)) * self$options$thVal ## custom by user
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
          f_points <- list(F1=NA, Fi=NA, F2=NA)
          p_points <- list(P1=NA, Pi=NA, P2=NA)
          l_points <- list(OGF0=NA, tang=NA, thres=NA)
          a_points <- list(OGF3=NA, PDA=NA)
        }
        
        ##### Results -----
        
        ## Add dependent variable name in the section title
        self$results$text$setTitle(paste("Results for", self$options$dep))

        ## Estimated Parameters and Model Evaluation tables
        pTable <- self$results$pTable
        tableFit <- self$results$fitq
        
        if (fit$convInfo$isConv) {
          for (i in seq_len(n_params)) { # one new row per parameter
            pTable$addRow(rowKey=i, values=list(
              var = if (i==1) dep, # variable name
              Parameter=params_df$param[i],
              Estimate=format(round(params_df$estim[i], 3), nsmall=2),
              Lower=format(round(params_df$lower[i], 3), nsmall=2),
              Upper=format(round(params_df$upper[i], 3), nsmall=2),
              SE=format(round(params_df$se[i], 3), nsmall=2),
              Statistics=format(round(params_df$t[i], 3), nsmall=2),
              pvalue=params_df$p[i]
            ))
          }
          pTable$setNote("conv", "This version fits the curve to means at each time point. Uncertainty metrics are approximate and may be optimistically biased under heteroscedasticity.", init=FALSE)
          
          tableFit$setRow(rowNo=1, values=list(
            AIC=format(round(AIC, 2), nsmall=2),
            AICc=format(round(AICc, 2), nsmall=2),
            BIC=format(round(BIC, 2), nsmall=2),
            R2=format(round(R2, 3), nsmall=3),
            R2_adj=format(round(R2_adj, 3), nsmall=3),
            RMSE=format(round(RMSE, 3), nsmall=3),
            MAE=format(round(MAE, 3), nsmall=3),
            MedAE=format(round(MedAE, 3), nsmall=3),
            sMAPE=sprintf("%.3f%%", round(sMAPE, 3)),
            RRMSE=sprintf("%.3f%%", round(RRMSE*100, 3))
          ))
          tableFit$setNote("conv", "This version fits the curve to means at each time point. Uncertainty metrics are approximate and may be optimistically biased under heteroscedasticity.", init=FALSE)
        } else {
          pTable$setNote("conv", "Model didn't converge.", init=FALSE)
          tableFit$setNote("conv", "Model didn't converge.", init=FALSE)
        }
        
        ## Key Growth Points table
        tableFp <- self$results$fpoints
        if (!length(zero_acc)==0) {
          tableFp$addRow(rowKey=self$options$time, values=list(
            OGF0=format(round(l_points$OGF0, 2), nsmall=2),
            Tangent=format(round(l_points$tang, 2), nsmall=2),
            Threshold=format(round(l_points$thres, 2), nsmall=2),
            F1=format(round(f_points$F1, 2), nsmall=2),
            Fi=format(round(f_points$Fi, 2), nsmall=2),
            F2=format(round(f_points$F2, 2), nsmall=2),
            P1=format(round(p_points$P1, 2), nsmall=2),
            Pi=format(round(p_points$Pi, 2), nsmall=2),
            P2=format(round(p_points$P2, 2), nsmall=2),
            OGF3=format(round(a_points$OGF3, 2), nsmall=2),
            PDA=format(round(a_points$PDA, 2), nsmall=2)
          ))
          tableFp$addRow(rowKey=self$options$dep, values=list(
            OGF0=format(round(W(l_points$OGF0), 2), nsmall=2),
            Tangent=format(round(W(l_points$tang), 2), nsmall=2),
            Threshold=format(round(W(l_points$thres), 2), nsmall=2),
            F1=format(round(W(f_points$F1), 2), nsmall=2),
            Fi=format(round(W(f_points$Fi), 2), nsmall=2),
            F2=format(round(W(f_points$F2), 2), nsmall=2),
            P1=format(round(W(p_points$P1), 2), nsmall=2),
            Pi=format(round(W(p_points$Pi), 2), nsmall=2),
            P2=format(round(W(p_points$P2), 2), nsmall=2),
            OGF3=format(round(W(a_points$OGF3), 2), nsmall=2),
            PDA=format(round(W(a_points$PDA), 2), nsmall=2)
          ))
        } else {
          tableFp$setNote("sig", "No inflection point found. Input data might not follow a sigmoidal trend.", init=FALSE)
        }
        
        ##### Plots Data #####
        
        ## Data for next functions
        private$.prep_mplot(x_raw, y_raw, t_new, W_pred, OGF_pred, 
                            OGF3_pred, f_points, p_points, l_points, a_points)
        private$.prep_dplot(t_new, W1_pred, W2_pred)
        
      }, # close .run
      
      ##### Plots Functions #####
      
      # Prepare data for Model Plot
      .prep_mplot = function(t, y, t_new, W, OGF, OGF3, Fp, Pp, Lp, Ap) {
        ## Set plot with model data as dataframe
        image <- self$results$mplot
        image$setState(list(
          data_p = data.frame(t=t, y=y),
          data_m = data.frame(t_new=t_new, W=W, OGF=OGF, OGF3=OGF3), 
          f_points = Fp, 
          p_points = Pp,
          l_points = Lp,
          a_points = Ap
          )
        )
      }, # close .prep_mplot
      
      # Model Plot Function
      .mplot = function(image, ...) {
        
        ## Check if there is any data to plot
        if (is.null(image$state))
          return(FALSE)
        
        ## Define data variables
        data_p <- image$state$data_p
        data_m <- image$state$data_m
        f_points <- image$state$f_points
        p_points <- image$state$p_points
        l_points <- image$state$l_points
        a_points <- image$state$a_points
        
        ## Deactivate the default box around plotting area
        par(bty='L')
        
        ## Better limits to avoid "out-of-axis" points
        better_lim <- function(vec, n=5) {
          ticks <- pretty(vec, n=n)
          c(min(ticks), max(ticks))
        }
        
        ## Redefine plot margins to fit axis labels
        par(mar = c(5.1, 3.1, 4.1, 3.1))  # c(5.1, 4.1, 4.1, 2.1)
        par(mgp = c(2, 0.5, 0)) # c(3, 1, 0)
        
        ## Create plot with samples points
        plot(data_p$t, data_p$y, col="grey", pch=16, 
             xlab="", ylab=self$options$dep,
             ylim=better_lim(data_p$y))
        ### Add W(t_new) curve
        lines(data_m$t_new, data_m$W, col="red", lwd=2)
        legends = c("Fitted curve") # legend label
        l_cols = c("red") # legend color
        l_ltys = c("solid") # legend line type
        
        ## Add Key Growth Points
        if (self$options$keyGrowth) {
          ## Function to retrieve aprox. y at gowth curve from an x value
          y_ <- function(x_) {
            return(approx(data_m$t_new, data_m$W, xout=x_, rule=2)$y)
          }
          ## Add A-Points vertical lines
          if ("ogf0" %in% self$options$lagEnd) {
            x_point <- l_points$OGF0
            y_point <- y_(x_point)
            segments(x_point, par("usr")[3], x_point, y_point, col="slategrey", lwd=2, lty=4)
            points(x_point, y_point, pch=19, col="black")
            text(x_point, y_point, labels="F0", pos=3, offset=1, col="black", font=2)
          }
          if ("tang" %in% self$options$lagEnd) {
            x_point <- l_points$tang
            y_point <- y_(x_point)
            segments(x_point, par("usr")[3], x_point, y_point, col="slategrey", lwd=2, lty=4)
            points(x_point, y_point, pch=19, col="black")
            text(x_point, y_point, labels="Tan.", pos=3, offset=1, col="black", font=2)
          }
          if ("thres" %in% self$options$lagEnd) {
            x_point <- l_points$thres
            y_point <- y_(x_point)
            segments(x_point, par("usr")[3], x_point, y_point, col="slategrey", lwd=2, lty=4)
            points(x_point, y_point, pch=19, col="black")
            text(x_point, y_point, labels="Thr.", pos=3, offset=1, col="black", font=2)
          }
          ## Add F-Points vertical lines
          if ("ogfMax" %in% self$options$fPoints) {
            x_point <- f_points$F1
            y_point <- y_(x_point)
            segments(x_point, par("usr")[3], x_point, y_point, col="slategrey", lwd=2, lty=4)
            points(x_point, y_point, pch=19, col="black")
            text(x_point, y_point, labels="F1", pos=3, offset=1, col="black", font=2)
          }
          if ("ogfI" %in% self$options$fPoints) {
            x_point <- f_points$Fi
            y_point <- y_(x_point)
            segments(x_point, par("usr")[3], x_point, y_point, col="slategrey", lwd=2, lty=4)
            points(x_point, y_point, pch=19, col="black")
            text(x_point, y_point, labels="Fi", pos=3, offset=1, col="black", font=2)
          }
          if ("ogfMin" %in% self$options$fPoints) {
            x_point <- f_points$F2
            y_point <- y_(x_point)
            segments(x_point, par("usr")[3], x_point, y_point, col="slategrey", lwd=2, lty=4)
            points(x_point, y_point, pch=19, col="black")
            text(x_point, y_point, labels="F2", pos=3, offset=1, col="black", font=2)
          }
          ## Add P-Points vertical lines
          if ("accMax" %in% self$options$pPoints) {
            x_point <- p_points$P1
            y_point <- y_(x_point)
            segments(x_point, par("usr")[3], x_point, y_point, col="slategrey", lwd=2, lty=4)
            points(x_point, y_point, pch=19, col="black")
            text(x_point, y_point, labels="P1", pos=3, offset=1, col="black", font=2)
          }
          if ("accI" %in% self$options$pPoints) {
            x_point <- p_points$Pi
            y_point <- y_(x_point)
            segments(x_point, par("usr")[3], x_point, y_point, col="slategrey", lwd=2, lty=4)
            points(x_point, y_point, pch=19, col="black")
            text(x_point, y_point, labels="Pi", pos=3, offset=1, col="black", font=2)
          }
          if ("accMin" %in% self$options$pPoints) {
            x_point <- p_points$P2
            y_point <- y_(x_point)
            segments(x_point, par("usr")[3], x_point, y_point, col="slategrey", lwd=2, lty=4)
            points(x_point, y_point, pch=19, col="black")
            text(x_point, y_point, labels="P2", pos=3, offset=1, col="black", font=2)
          }
          ## Add Asymptote vertical lines
          if ("ogf3" %in% self$options$asymptote) {
            x_point <- a_points$OGF3
            y_point <- y_(x_point)
            segments(x_point, par("usr")[3], x_point, y_point, col="slategrey", lwd=2, lty=4)
            points(x_point, y_point, pch=19, col="black")
            text(x_point, y_point, labels="F3", pos=3, offset=1, col="black", font=2)
          }
          if ("pda" %in% self$options$asymptote) {
            x_point <- a_points$PDA
            y_point <- y_(x_point)
            segments(x_point, par("usr")[3], x_point, y_point, col="slategrey", lwd=2, lty=4)
            points(x_point, y_point, pch=19, col="black")
            text(x_point, y_point, labels="PDA", pos=3, offset=1, col="black", font=2)
          }
        }
        
        ## Add OGF(t_new) curve
        if (self$options$ogf) {
          par(new = TRUE) # new plot, same area
          plot(data_m$t_new, data_m$OGF, type="l", col="darkgreen",axes=FALSE,
               xlab="", ylab="", ylim=better_lim(data_m$OGF))
          axis(side=4) # add OGF y-axis on right side of the plot
          mtext(paste0("OGF (", self$options$dep, "² / ", self$options$time, "³)"), 
                side=4, line=2) # y-axis label
          legends = c(legends, "Ontogenetic growth force") # legend label
          l_cols = c(l_cols, "darkgreen") # legend color
          l_ltys = c(l_ltys, "solid") # legend line type
        }
        
        ## Add OGF'''(t_new) curve
        if (self$options$ogf_s) {
          par(new = TRUE) # new plot, same area
          plot(data_m$t_new, data_m$OGF3, type="l", col="lightgreen",axes=FALSE,
               xlab="", ylab="", lty="dashed", ylim=better_lim(data_m$OGF3))
          legends = c(legends, "OGF 3rd Derivative") # legend label
          l_cols = c(l_cols, "lightgreen") # legend color
          l_ltys = c(l_ltys, "dashed") # legend line type
        }
        
        ## Add lines legend
        legend("bottom", inset=-0.2, legend=legends, col=l_cols, lty=l_ltys, 
               horiz=TRUE, bty="n", lwd=2, xpd=TRUE)
        
        ## Notify the rendering system that we have plotted something
        TRUE
        
      }, # close .mplot
      
      # Prepare data for Derivative Plot
      .prep_dplot = function(t_new, W1, W2) {
        ## Set plot with model data as dataframe
        image <- self$results$dplot
        image$setState(data.frame(t_new=t_new, W1=W1, W2=W2))
      }, # close .prep_dplot
      
      # Create the derivative plot function
      .dplot = function(image, ...) {
        
        ## Check if there is any data to plot
        if (is.null(image$state))
          return(FALSE)
        
        ## Deactivate the default box around plotting area
        par(bty='L')
        
        ## Better limits to avoid "out-of-axis" points
        better_lim <- function(vec, n=5) {
          ticks <- pretty(vec, n = n)
          c(min(ticks), max(ticks))
        }
        
        ## Redefine plot margins to fit axis labels
        par(mar = c(5.1, 3.1, 4.1, 3.1))  # c(5.1, 4.1, 4.1, 2.1)
        par(mgp = c(2, 0.5, 0)) # c(3, 1, 0)
        
        ## Create plot with W'(t_new) curve
        plot(image$state$t_new, image$state$W1, type="l", col = "yellow3",
             xlab="", ylab=paste0("Growth rate (", self$options$dep, " / ", self$options$time, ")"), 
             ylim=better_lim(image$state$W1))
        legends = c("Growth rate") # legend label
        l_cols = c("yellow3") # legend color
        
        ## Add W''(t_new) curve
        par(new = TRUE) # new plot, same area
        plot(image$state$t_new, image$state$W2, type="l", col="navy", axes=FALSE, 
             xlab="", ylab="", ylim=better_lim(image$state$W2))
        axis(side = 4) # add acceleration y-axis on right side of the plot
        mtext(paste0("Acceleration (", self$options$dep, " / ", self$options$time, "²)"), 
              side=4, line=2) # y-axis label
        legends = c(legends, "Acceleration") # legend label
        l_cols = c(l_cols, "navy") # legend color
        
        ## Add lines legend
        legend("bottom", inset=-0.2, legend=legends, col=l_cols, lty=1, 
               horiz=TRUE, bty="n", lwd=2, xpd=TRUE)
        
        ## Notify the rendering system that we have plotted something
        TRUE
        
      } # close .dplot
      
    ) # close list
) # close R6Class
