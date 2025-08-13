
# This file is a generated template, your changes will not be overwritten

scurveClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "scurveClass",
    inherit = scurveBase,
    private = list(
      
      ##### Main Function #####
      
      .run = function() {
        
        ##### Get User Data #####
        
        ## Option values into shorter variable names
        dep  <- self$options$dep
        time <- self$options$time
        
        ## Check if variables have any data
        if (is.null(dep) || is.null(time))
          return()
        
        ## Get the data
        data <- self$data
        
        ## Convert to appropriate data types
        data[[dep]] <- jmvcore::toNumeric(data[[dep]])
        data[[time]] <- jmvcore::toNumeric(data[[time]])
        
        ## Remove NAs
        data <- na.omit(data)
        
        ## Define data variables
        y_raw <- data[[dep]]
        x_raw <- data[[time]]
        
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
        
        ## 3rd and 4th derivatives for PDA
        W3_expr <- D(W2_expr, "x") # 3rd derivative
        W4_expr <- D(W3_expr, "x") # 4th derivative
        
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
        W4 <- function(x) eval(W4_expr) # W''''(t) (PDA)
        OGF <- function(x) eval(OGF_expr) # OGF(t)
        OGF3 <- function(x) eval(OGF3_expr) # OGF'''(t)
        
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
        
        ##### Key Growth Points #####
        
        ## Length of x-axis
        len <- length(t_new)
        
        ## Check if the data have an inflection point
        zero_acc <- which(diff(sign(W2_pred))!=0)
        if (length(zero_acc)==0) {
          self$results$fpoints$setNote("sig", "No inflection point found. Input data might not follow a sigmoidal trend.", init=FALSE)
          f_points <- list(F1=NA, Fi=NA, F2=NA)
          p_points <- list(P1=NA, Pi=NA, P2=NA)
          l_points <- list(OGF0=NA, tang=NA, thres=NA)
          a_points <- list(OGF3=NA, PDA=NA)
        } else {
          
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
          ### threshold method
          y_vals <- sapply(t_new, W)
          y_thresh <- W(max(t_new)) * self$options$thVal ## custom by user
          t_lag_threshold <- t_new[min(which(y_vals >= y_thresh))]
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
        
        }
        
        ##### Model Information #####
        
        ## Add dependent variable name in the section title
        self$results$text$setTitle(paste("Results for", self$options$dep))
        
        ## Show model equation as text
        str_expr <- paste(deparse(s_expr[[1]]), collapse="")
        eq_str <- paste0("\nModel equation:\n", str_expr)
        self$results$text$setContent(eq_str)

        ## Estimated parameters table
        pTable <- self$results$pTable
        pTable$setRow(rowNo=1, values=list(
          A=format(round(params$A, 2), nsmall=2),
          d=format(round(params$d, 2), nsmall=2),
          K=format(round(params$K, 2), nsmall=2),
          Ti=format(round(params$Ti, 2), nsmall=2)
        ))
        
        ## Evaluation Metrics table
        tableFit <- self$results$fitq
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
          RRMSE=sprintf("%.3f%%", round(RRMSE*100, 3)),
          fdf1=round(df1),
          fdf2=round(df2),
          f=format(round(F_s, 3), nsmall=3),
          fp=F_p
        ))
        
        ## Key Growth Points table
        tableFp <- self$results$fpoints
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
        legends = c("Adjusted function") # legend label
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
