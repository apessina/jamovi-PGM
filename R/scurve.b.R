scurveClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "scurveClass",
    inherit = scurveBase,
    private = list(
      
      .run = function() {
        
        ##### Data -----
        
        ## Option values into shorter variable names
        deps  <- self$options$deps
        time <- self$options$time
        ids <- NULL # self$options$ids
        
        ## Check if variables have any data and get the data
        if (!is.null(deps) & !is.null(time)) {
          if (!is.null(ids)) {
            data <- self$data[, c(deps, time, ids), drop=FALSE]
          } else {
            data <- self$data[, c(deps, time), drop=FALSE]
          }
        } else {
          return()
        }
        
        ## Convert time to appropriate data type
        data[[time]] <- jmvcore::toNumeric(data[[time]])
        
        ##### Curve Resolution -----
        
        ## Giving time a new "resolution" for smooth lines
        t_raw <- data[[time]]
        n_points <- min(3000, max(500, round((max(t_raw) - min(t_raw)) * self$options$res)))
        t_global <- seq(min(t_raw), max(t_raw), length.out=n_points)
        
        ##### deps Iteration -----
        raw_dt <- data.frame(t = t_raw)
        deps_dt <- list(t = t_global)
        c_points <- list()
        auc_vals <- list()
        
        for (dep in deps) {
          
          ## Convert dep to appropriate data type
          data[[dep]] <- jmvcore::toNumeric(data[[dep]])  
          
          ## Remove NAs
          data_cln <- na.omit(data[, c(time, dep)])
          
          ## Convert to appropriate data types
          x_raw <- data_cln[[time]]
          y_raw <- data_cln[[dep]]
          
          ## Check if both variables have at least two unique values
          if (length(unique(x_raw)) < 2 || length(unique(y_raw)) < 2)
            stop("At least two unique (x, y) pairs are required for growth analysis")
  
          ##### Model sup data -----
          
          agg_vals <- aggregate(y_raw, by = list(x_raw), FUN = mean)
          t <- agg_vals[[1]]
          y <- agg_vals[[2]]
          
          ##### Model specs -----
          
          spec <- get_model_spec(self$options$model, t, y)
          
          self$results$equation$setContent(paste0(
            '<div style="margin-top:10px;">',
            '<div style="font-family:Segoe UI, sans-serif; font-size:12px; margin-bottom:3px;">Model Equation</div>',
            '<p style="text-align:center; font-family:Cambria, serif; font-size:16px;">',
            spec$html_eq,
            '</p>',
            '</div>'
          ))
          
          expr  <- spec$expr
          init  <- spec$init
          lower <- spec$lower
          upper <- spec$upper
          
          ##### Modeling -----
          
          ## Initial data structure for modeling
          data_fit <- data.frame(y=as.numeric(y_raw), t=as.numeric(x_raw))
          wfunc <- "All weights = 1"
          w <- rep(1, length(data_fit$y))
          
          ## Check if initial values return finite output
          start_env <- c(list(t = data_fit$t), init)
          start_val <- eval(expr, envir = start_env)
          if (any(!is.finite(start_val))) {
            stop("ERROR: Initial values generate non-finite model predictions")
          }
          
          # Check: (t - mu)/delta must stay positive for G. Weibull
          if (self$options$model == "genweibull") {
            z0 <- (data_fit$t - init$mu) / init$delta
            if (any(z0 <= 0)) {
              stop("ERROR: Invalid initial domain for G. Weibull")
            }
          }
          
          ## Fit selected model
          fit <- gsl_nls(
            fn = as.formula(call("~", quote(y), expr), env = .GlobalEnv),
            data = data_fit,
            start = init,
            lower = lower,
            upper = upper,
            algorithm = "lm"
          )
          
          ## Apply error weights if requested
          if (self$options$wtype != "none") {
            
            ### FGLS / IRLS
            eps   <- 1e-6
            max_iter <- 20
            for (k in 1:max_iter) {
              mu <- pmax(fitted(fit), eps)
              r  <- resid(fit)
              r2 <- pmax(r^2, eps)
              
              if (self$options$wtype == "power") {
                m <- lm(log(r2) ~ log(mu))
                delta <- coef(m)[2] / 2
                w <- mu^(-2 * delta)
                wfunc <- paste0("Var(ε) = σ² * μ^(2δ), δ=",round(delta, 2),
                                ", ",k," of ",max_iter," IRLS Iterations")
                
              } else if (self$options$wtype == "exp") {
                mu_s <- (mu - mean(mu)) / sd(mu)
                m <- lm(log(r2) ~ mu_s)
                delta_s <- coef(m)[2] / 2
                w <- exp(-2 * delta_s * mu_s)
                wfunc <- paste0("Var(ε) = σ² * exp(2δ·z(μ)), δ=",round(delta_s, 2),
                                ", ",k," of ",max_iter," IRLS Iterations")
                
              } else if (self$options$wtype == "tpoly") { 
                t0 <- pmax(data_fit$t, eps)
                m <- lm(log(r2) ~ poly(t0, 3, raw = TRUE))
                vhat <- exp(predict(m, newdata = data.frame(t0 = t0)))
                vhat <- pmax(vhat, quantile(vhat, 0.02, na.rm = TRUE))
                w <- 1 / vhat
                wfunc <- paste0("Var(ε) = σ² * f(t), ",
                                k," of ",max_iter," IRLS Iterations")
             
              } else {
                stop("ERROR: invalid weight type")
              }
              
              w <- w / median(w, na.rm = TRUE)
              cap <- quantile(w, 0.99, na.rm = TRUE)
              w[w > cap] <- cap
              floor_w <- quantile(w, 0.01, na.rm = TRUE)
              w[w < floor_w] <- floor_w
              data_fit$w <- w
              
              beta_old <- coef(fit)
              fit <- gsl_nls(
                fn = as.formula(call("~", quote(y), expr), env = .GlobalEnv),
                data = data_fit,
                lower = lower,
                upper = upper,
                algorithm = "lm",
                start = as.list(beta_old),
                weights = w
              )
              
              if (max(abs((coef(fit)-beta_old) / pmax(abs(beta_old), eps))) < 1e-3) break
            }
            
          }
          
          ## Model intra-individual correlation if ids
          icorr <- !is.null(ids) && self$options$cor_type != "none"
          if (icorr) {
            # to-do
          }
          
          ## Get estimated parameters
          beta_hat <- coef(fit)
          
          ## Model equation with fit parameters as expression
          W_expr <- do.call("substitute", list(expr, as.list(beta_hat)))
          
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
          W  <- function(t) eval(W_expr, envir = list(t = t)) # W(t) 
          W1 <- function(t) eval(W1_expr, envir = list(t = t)) # W'(t)
          W2 <- function(t) eval(W2_expr, envir = list(t = t)) # W''(t)
          W4 <- function(t) eval(W4_expr, envir = list(t = t)) # W''''(t) (PDA)
          OGF <- function(t) eval(OGF_expr, envir = list(t = t)) # OGF(t)
          OGF3 <- function(t) eval(OGF3_expr, envir = list(t = t)) # OGF'''(t)
          
          ##### Parameters -----
          
          ## Number of samples
          n <- nrow(data_fit)
          
          ## Number of parameters
          n_params <- length(beta_hat)
          
          ## Degrees of freedom
          df2 <- n - n_params
  
          ## Covariance Matrix
          covm <- try(vcov(fit), silent = TRUE)
  
          ## Standard Errors 
          if (inherits(covm, "try-error")) {
            se_vals <- rep(NA_real_, length(beta_hat))
          } else {
            se_vals <- sqrt(diag(covm))
          }
          
          ## Inference - Wald test
          stat_vals <- beta_hat / se_vals
          
          if (self$options$wtype == "none") {
            p_vals <- 2 * pt(abs(stat_vals), df = df2, lower.tail = FALSE)
            crit <- qt(1 - 0.05/2, df = df2)
            infer_note <- "Wald t (iid)"
          } else {
            p_vals <- 2 * pnorm(abs(stat_vals), lower.tail = FALSE)
            crit <- qnorm(1 - 0.05/2)
            infer_note <- "Wald z (FGLS approx)"
          }
          
          ## Confidence Intervals
          lower_ci <- beta_hat - crit * se_vals
          upper_ci <- beta_hat + crit * se_vals
          
          ## Parameters DF
          params_df <- data.frame(
            param=names(beta_hat),
            estim=beta_hat,
            lower=lower_ci,
            upper=upper_ci,
            se=se_vals,
            stat=stat_vals,
            p=p_vals
          )
          
          ##### Evaluation -----
          
          ## Residuals / Error
          res <- resid(fit)
          
          ## Fitted values
          yhat <- fitted(fit)
          
          ## Sum of squares
          SSt <- sum((y_raw - mean(y_raw))^2)
          SSe <- sum(res^2)
          
          ## Goodness-of-fit metrics
          R2 <- 1 - SSe / SSt
          R2_adj <- 1 - (SSe/df2) / (SSt/(n - 1))
          
          if (self$options$wtype == "none") {
            AIC <- n * log(SSe/n) + 2 * n_params
            AICc <- AIC + (2 * n_params * (n_params + 1)) / (n - n_params - 1)
            BIC <- n * log(SSe/n) + log(n) * n_params
          } else { # normal with Var(e_i) = sigma2 / w_i (w = 1/Var relative)
            w_i <- data_fit$w
            sigma2_hat <- sum(data_fit$w * res^2) / n
            logLik <- -0.5 * sum(
              log(2*pi*sigma2_hat / w_i) + (w_i * res^2) / sigma2_hat
            )
            AIC <- -2 * logLik + 2 * n_params
            BIC <- -2 * logLik + log(n) * n_params
            AICc <- AIC + (2 * n_params * (n_params + 1)) / (n - n_params - 1)
          }
          
          ## Error metrics
          MSE  <- sum(w * res^2) / sum(w)
          RMSE <- sqrt(MSE)
          MAE <- mean(abs(res))
          MedAE <- median(abs(res))
          sMAPE <- mean(2 * abs(res) / (abs(y_raw) + abs(yhat))) * 100
          RRMSE <- RMSE / weighted.mean(y_raw, w)
          
          ##### dep Resolution -----
          
          ## Trim t_global unit the max value of t.
          ### This allows graphical comparison in the same scale,
          ### even if measures stopped at a certain age for one variable
          t_new <- t_global[t_global >= min(t) & t_global <= max(t)]
          
          ## Apply Data Modeling functions to t_new
          W_pred <- W(t_new) # vector of weight over time
          W1_pred <- W1(t_new) # vector of growth rate over time
          W2_pred <- W2(t_new) # vector of acceleration over time
          W4_pred <- W4(t_new) # vector of 4th derivative over time (PDA)
          OGF_pred <- OGF(t_new) # vector of OGF over time
          OGF3_pred <- OGF3(t_new) # vector of the OGF 3rd Derivative over time
          
          ## Calculate AUC
          if (self$options$intInterval) {
            auc_lower <- self$options$intL
            auc_upper <- self$options$intU
          } else {
            auc_lower <- min(t_new, na.rm=TRUE)
            auc_upper <- max(t_new, na.rm=TRUE)
          }
          auc_vals[[dep]] <- integrate(W, lower=auc_lower, upper=auc_upper)$value
          
          ##### Critical Points -----
          cp <- compute_critical_points(
            t_new = t_new,
            W1_pred = W1_pred,
            W2_pred = W2_pred,
            OGF_pred = OGF_pred,
            OGF3_pred = OGF3_pred,
            W = W, W2 = W2, W4 = W4,
            OGF = OGF, OGF3 = OGF3,
            thVal = self$options$thVal
          )
          c_points[[dep]] <- cp
          
          ##### Results -----
          
          ## Estimation Info
          eTable <- self$results$eTable
          eTable$addRow(
            rowKey=paste0(dep,"L1"), 
            values=list(
              var=dep,
              infoh="Estimator",
              infor="Nonlinear Least Squares"
            )
          )
          eTable$addRow(
            rowKey=paste0(dep,"L2"), 
            values=list(
              infoh="Optimization",
              infor="Levenberg–Marquardt"
            )
          )
          eTable$addRow(
            rowKey=paste0(dep,"L3"), 
            values=list(
              infoh="Weight Function",
              infor=wfunc
            )
          )
  
          ## Estimated Parameters and Model Evaluation
          pTable <- self$results$pTable
          gofTable <- self$results$meval$gof
          emTable <- self$results$meval$em
          
          if (fit$convInfo$isConv) {
            
            for (i in seq_len(n_params)) {
              pTable$addRow(
                rowKey = paste0(dep,"L",i),
                values = list(
                  var = if (i==1) dep,
                  Parameter=params_df$param[i],
                  Estimate=format(round(params_df$estim[i], 3), nsmall=2),
                  Lower=format(round(params_df$lower[i], 3), nsmall=2),
                  Upper=format(round(params_df$upper[i], 3), nsmall=2),
                  SE=format(round(params_df$se[i], 3), nsmall=2),
                  Statistics=format(round(params_df$stat[i], 3), nsmall=2),
                  pvalue=params_df$p[i]
                )
              )
            }
            pTable$setNote("conv", infer_note, init=FALSE)
            
            gofTable$setRow(rowKey=dep, values=list(
              AIC=format(round(AIC, 2), nsmall=2),
              AICc=format(round(AICc, 2), nsmall=2),
              BIC=format(round(BIC, 2), nsmall=2),
              R2=format(round(R2, 3), nsmall=3),
              R2_adj=format(round(R2_adj, 3), nsmall=3)
            ))
            
            emTable$setRow(rowKey=dep, values=list(
              RMSE=format(round(RMSE, 3), nsmall=3),
              MAE=format(round(MAE, 3), nsmall=3),
              MedAE=format(round(MedAE, 3), nsmall=3),
              sMAPE=sprintf("%.3f%%", round(sMAPE, 3)),
              RRMSE=sprintf("%.3f%%", round(RRMSE*100, 3))
            ))
            
          } else {
            pTable$addRow(rowKey=1)
            pTable$setNote("conv", "Model didn't converge.", init=FALSE)
            gofTable$setNote("conv", "Model didn't converge.", init=FALSE)
            emTable$setNote("conv", "Model didn't converge.", init=FALSE)
          }
          
          ## AUC table
          aucTable <- self$results$aucTable
          aucTable$setRow(rowKey=dep, values=list(
            Lower=format(round(auc_lower, 2), nsmall=2),
            Upper=format(round(auc_upper, 2), nsmall=2),
            AUC=auc_vals[[dep]]
          ))
          
          ## Key Growth Points table
          tableFp <- self$results$fpoints
          if (cp$inflection) {
            tableFp$addRow(rowKey=self$options$time, values=list(
              F0=format(round(cp$l_points$OGF0, 2), nsmall=2),
              Tangent=format(round(cp$l_points$tang, 2), nsmall=2),
              Threshold=format(round(cp$l_points$thres, 2), nsmall=2),
              F1=format(round(cp$f_points$F1, 2), nsmall=2),
              Fi=format(round(cp$f_points$Fi, 2), nsmall=2),
              F2=format(round(cp$f_points$F2, 2), nsmall=2),
              P1=format(round(cp$p_points$P1, 2), nsmall=2),
              Pi=format(round(cp$p_points$Pi, 2), nsmall=2),
              P2=format(round(cp$p_points$P2, 2), nsmall=2),
              F3=format(round(cp$a_points$OGF3, 2), nsmall=2),
              PDA=format(round(cp$a_points$PDA, 2), nsmall=2)
            ))
            tableFp$addRow(rowKey=dep, values=list(
              F0=format(round(W(cp$l_points$OGF0), 2), nsmall=2),
              Tangent=format(round(W(cp$l_points$tang), 2), nsmall=2),
              Threshold=format(round(W(cp$l_points$thres), 2), nsmall=2),
              F1=format(round(W(cp$f_points$F1), 2), nsmall=2),
              Fi=format(round(W(cp$f_points$Fi), 2), nsmall=2),
              F2=format(round(W(cp$f_points$F2), 2), nsmall=2),
              P1=format(round(W(cp$p_points$P1), 2), nsmall=2),
              Pi=format(round(W(cp$p_points$Pi), 2), nsmall=2),
              P2=format(round(W(cp$p_points$P2), 2), nsmall=2),
              F3=format(round(W(cp$a_points$OGF3), 2), nsmall=2),
              PDA=format(round(W(cp$a_points$PDA), 2), nsmall=2)
            ))
          } else {
            tableFp$setNote("sig", "No inflection point found. Input data might not follow a sigmoidal trend.", init=FALSE)
          }
          
          ##### Fitted Mean CI -----
          ### Delta method (Wald approx.)
          if (!inherits(covm, "try-error")) {
            grad_exprs <- lapply(names(beta_hat), function(p) D(expr, p))
            grad_fun <- function(t0) {
              sapply(grad_exprs, function(g)
                eval(g, envir = c(list(t = t0), as.list(beta_hat))))
            }
            G <- t(sapply(t_new, grad_fun))
            var_fit <- rowSums((G %*% covm) * G)
            se_fit <- sqrt(pmax(var_fit, 0))
            W_lower <- W_pred - crit * se_fit
            W_upper <- W_pred + crit * se_fit
          } else {
            W_lower <- rep(NA_real_, length(t_new))
            W_upper <- rep(NA_real_, length(t_new))
          }
        
          ##### Global Lists ----
          ## Fill length with NA to plot trimmed dep lengths
          raw_dt[[dep]] <- c(
            y_raw, rep(NA, length(t_raw)-length(y_raw))
          )
          deps_dt[[dep]] <- data.frame(
            W_pred = c(
              W_pred, rep(NA, length(t_global)-length(W_pred))
            ),
            W_lower = c(
              W_lower, rep(NA, length(t_global)-length(W_lower))
            ),
            W_upper = c(
              W_upper, rep(NA, length(t_global)-length(W_upper))
            )
          )
          
          ##### Render Data -----
          
          ## Single model plot
          mplot <- self$results$mplot$addItem(key=dep)
          mplot$setTitle(paste("Growth Curve for", dep))
          mplot$setState(list(
            dep = dep,
            data_p = data.frame(t = x_raw, y = y_raw),
            data_m = data.frame(t_new = t_new, W = W_pred, 
                                W_lower = W_lower, W_upper = W_upper, 
                                OGF = OGF_pred, OGF3 = OGF3_pred),
            f_points = cp$f_points,
            p_points = cp$p_points,
            l_points = cp$l_points,
            a_points = cp$a_points
          ))
          
          ## Derivatives plot
          if (self$options$sndPlot) {
            dplot <- self$results$dplot$addItem(key=dep)
            dplot$setTitle(paste("Derivatives from", dep))
            dplot$setState(list(
              dep = dep,
              t_new = t_new,
              W1 = W1_pred,
              W2 = W2_pred
            ))
          }
          
          ## Residuals plot
          if (self$options$resPlot) {
            resplot <- self$results$resplot$addItem(key=dep)
            resplot$setTitle(paste("Residuals from", dep))
            resplot$setState(list(
              dep = dep,
              t = x_raw,
              w = w,
              res = res,
              fitted = fitted(fit),
              dfres = df.residual(fit)
            ))
          }
          
        } # close loop
        
        ##### Global Info -----
        
        ## Rank AUC table
        auc_num <- unlist(auc_vals)
        rank_val <- (auc_num / max(auc_num)) * 100
        for (dep in deps) {
          aucTable$setRow(
            rowKey=dep, 
            values=list(
              Rank=sprintf("%.2f%%", rank_val[[dep]])
            )
          )
        }
        
        ## Multiple models plot
        if (self$options$multiPlot) {
          multiplot <- self$results$multiplot
          multiplot$setState(list(
            data_p = raw_dt,
            data_m = deps_dt,
            c_points = c_points
          ))
        }
        
      }, # close .run
      
      ##### Render Functions -----
      
      # Residuals plot
      .resplot = function(image, ...) {
        
        ## Check if there is any data to plot
        if (is.null(image$state))
          return(FALSE)
        
        t <- image$state$t
        w <- image$state$w
        r <- image$state$res
        f <- image$state$fitted
        dfres = image$state$dfres
        
        ## Standardized residuals (based on model sigma)
        eps <- 1e-12
        w2 <- pmax(w, eps)
        sse_w <- sum(w2 * r^2, na.rm = TRUE)
        sigma2_hat <- sse_w / max(dfres, 1)
        rs <- as.numeric(r * sqrt(w2) / sqrt(sigma2_hat))
        
        ## Start plot 2x2
        op <- par(mfrow = c(2, 2)); on.exit(par(op), add = TRUE)
        
        ### Residuals vs Fitted
        plot(f, rs, xlab = "Fitted", ylab = "Std. residuals",
             main = "Residuals vs Fitted"); abline(h = 0, lty = 2)
        
        ### QQ-plot
        qqnorm(rs, main = "Normal Q–Q (std. residuals)"); qqline(rs)
        
        ### Scale–Location
        plot(f, sqrt(abs(rs)), xlab = "Fitted", ylab = "√|std. residual|",
             main = "Scale–Location")
        
        ### Residuals vs time
        plot(t, rs, xlab = self$options$time, ylab = "Std. residuals",
             main = paste0("Residuals vs ", self$options$time))
        abline(h = 0, lty = 2)
        
        ## Notify the rendering system that we have plotted something
        TRUE
        
      }, # close .resplot
      
      # Model plot
      .mplot = function(image, ...) {
        
        ## Check if there is any data to plot
        if (is.null(image$state))
          return(FALSE)
        
        ## Define data variables
        dep <- image$state$dep
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
        plot(data_p$t, data_p$y, xlab="", ylab=dep,
             col=if (self$options$dataPoints) "grey" else NA, pch=16,
             ylim=better_lim(data_p$y))
        
        ### Fitted function CI
        if (self$options$ciPolygon) {
          polygon(
            x = c(data_m$t_new, rev(data_m$t_new)),
            y = c(data_m$W_lower, rev(data_m$W_upper)),
            border = NA,
            col = adjustcolor("blue", alpha.f = 0.20)
          )
        }
        
        ### Add W(t_new) curve
        lines(data_m$t_new, data_m$W, col="red", lwd=2)
        legends = c("Fitted curve") # legend label
        l_cols = c("red") # legend color
        l_ltys = c("solid") # legend line type
        
        ## Add Key Growth Points
        if (self$options$keyGrowth) {
          
          y_ <- function(x_) {
            approx(data_m$t_new, data_m$W, xout = x_, rule = 2)$y
          }
          
          draw_growth_point <- function(x_point, label) {
            y_point <- y_(x_point)
            segments(x_point, par("usr")[3], x_point, y_point,
                     col = "slategrey", lwd = 2, lty = 4)
            points(x_point, y_point, pch = 19, col = "black")
            text(x_point, y_point, labels = label,
                 pos = 3, offset = 1, col = "black", font = 2)
          }
          
          point_defs <- list(
            list(option_group = "lagEnd",     opt = "ogf0",   point = l_points$OGF0, label = "F0"),
            list(option_group = "lagEnd",     opt = "tang",   point = l_points$tang, label = "Tan."),
            list(option_group = "lagEnd",     opt = "thres",  point = l_points$thres, label = "Thr."),
            list(option_group = "fPoints",    opt = "ogfMax", point = f_points$F1,    label = "F1"),
            list(option_group = "fPoints",    opt = "ogfI",   point = f_points$Fi,    label = "Fi"),
            list(option_group = "fPoints",    opt = "ogfMin", point = f_points$F2,    label = "F2"),
            list(option_group = "pPoints",    opt = "accMax", point = p_points$P1,    label = "P1"),
            list(option_group = "pPoints",    opt = "accI",   point = p_points$Pi,    label = "Pi"),
            list(option_group = "pPoints",    opt = "accMin", point = p_points$P2,    label = "P2"),
            list(option_group = "asymptote",  opt = "ogf3",   point = a_points$OGF3,  label = "F3"),
            list(option_group = "asymptote",  opt = "pda",    point = a_points$PDA,   label = "PDA")
          )
          
          for (d in point_defs) {
            if (!is.null(d$point) && !is.na(d$point) && d$opt %in% self$options[[d$option_group]]) {
              draw_growth_point(d$point, d$label)
            }
          }
        }
        
        ## Add OGF(t_new) curve
        if (self$options$ogf) {
          par(new = TRUE) # new plot, same area
          plot(data_m$t_new, data_m$OGF, type="l", col="darkgreen",axes=FALSE,
               xlab="", ylab="", ylim=better_lim(data_m$OGF))
          axis(side=4) # add OGF y-axis on right side of the plot
          mtext(paste0("OGF (", dep, "² / ", self$options$time, "³)"), 
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
      
      # Multiple plot
      .multiplot = function(image, ...) {
        
        ## Check if there is any data to plot
        if (is.null(image$state))
          return(FALSE)
        
        data_p <- image$state$data_p
        data_m <- image$state$data_m
        c_points <- image$state$c_points
        
        ## Time columns and y columns names
        t_p <- data_p$t
        t_m <- data_m$t
        y_deps <- setdiff(names(data_p), "t")
        y_pred_deps <- setdiff(names(data_m), "t")
        
        ## Generate colors in the standard jamovi palette
        colors_pred <- jmvcore::colorPalette(n=length(y_deps))
        
        ## Better limits to avoid points without axis labels
        better_lim <- function(vec, n=5) {
          ticks <- pretty(vec, n=n)
          c(min(ticks), max(ticks))
        }
        
        ## Join all y data to calculate limits
        ylim_total <- range(unlist(data_p[y_deps]), unlist(data_m[y_pred_deps]), na.rm = TRUE)
        ylim_total <- better_lim(ylim_total)
        
        ## Deactivate the default box around the plotting area
        par(bty='L')
        
        ## Create and set-up the plot
        plot(t_p, data_p[[y_deps[1]]], pch=16, 
             xlab="", ylab="", ylim = ylim_total,
             col = adjustcolor(colors_pred[1], alpha.f = 0.3)
        )
        
        ## Add the other variables points
        if (length(y_deps) > 1) {
          for (i in seq_along(y_deps)[-1]) {
            dep_name <- y_deps[i]
            dep_p <- data_p[[dep_name]]
            points(t_p, dep_p, pch = 16,
                   col = adjustcolor(colors_pred[i], alpha.f = 0.3))
          }
        }
        
        ## Add predicted curves
        for (i in seq_along(y_pred_deps)) {
          dep_name <- y_pred_deps[i]
          dep_m <- data_m[[dep_name]]
          
          ### CI Interval
          if (self$options$ciPolygon) {
            polygon(
              x = c(t_m, rev(t_m)),
              y = c(dep_m$W_lower, rev(dep_m$W_upper)),
              border = NA,
              col = adjustcolor("blue", alpha.f = 0.20)
            )
          }
          
          ### W(t) curve
          lines(t_m, dep_m$W_pred, col = colors_pred[i], lwd=2)
          
          ## Add Key Growth Points
          if (self$options$keyGrowth) {
            
            l_points <- c_points[[dep_name]]$l_points
            f_points <- c_points[[dep_name]]$f_points
            p_points <- c_points[[dep_name]]$p_points
            a_points <- c_points[[dep_name]]$a_points
            
            y_ <- function(x_) {
              approx(t_m, dep_m$W_pred, xout = x_, rule = 2)$y
            }
            
            draw_growth_point <- function(x_point, label, col) {
              y_point <- y_(x_point)
              segments(x_point, par("usr")[3], x_point, y_point,
                       col = col, lwd = 2, lty = 4)
              points(x_point, y_point, pch = 19, col = col)
              text(x_point, y_point, labels = label,
                   pos = 3, offset = 1, col = col, font = 2)
            }
            
            point_defs <- list(
              list(option_group = "lagEnd",     opt = "ogf0",   point = l_points$OGF0, label = "F0"),
              list(option_group = "lagEnd",     opt = "tang",   point = l_points$tang, label = "Tan."),
              list(option_group = "lagEnd",     opt = "thres",  point = l_points$thres, label = "Thr."),
              list(option_group = "fPoints",    opt = "ogfMax", point = f_points$F1,    label = "F1"),
              list(option_group = "fPoints",    opt = "ogfI",   point = f_points$Fi,    label = "Fi"),
              list(option_group = "fPoints",    opt = "ogfMin", point = f_points$F2,    label = "F2"),
              list(option_group = "pPoints",    opt = "accMax", point = p_points$P1,    label = "P1"),
              list(option_group = "pPoints",    opt = "accI",   point = p_points$Pi,    label = "Pi"),
              list(option_group = "pPoints",    opt = "accMin", point = p_points$P2,    label = "P2"),
              list(option_group = "asymptote",  opt = "ogf3",   point = a_points$OGF3,  label = "F3"),
              list(option_group = "asymptote",  opt = "pda",    point = a_points$PDA,   label = "PDA")
            )
            
            for (d in point_defs) {
              if (!is.null(d$point) && !is.na(d$point) && d$opt %in% self$options[[d$option_group]]) {
                draw_growth_point(d$point, d$label, colors_pred[i])
              }
            }
          }
        }
        
        ## Make legend
        legend("bottom", inset=-0.2, legend=y_pred_deps, col=colors_pred, 
               lty=1, lwd=2, horiz=TRUE, bty="n", xpd=TRUE)
        
        ## Notify the rendering system that we have plotted something
        TRUE
        
      }, # close .multiplot
      
      # Derivatives plot
      .dplot = function(image, ...) {
        
        ## Check if there is any data to plot
        if (is.null(image$state))
          return(FALSE)
        
        dep <- image$state$dep
        
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
             xlab="", ylab=paste0("Growth rate (", dep, " / ", self$options$time, ")"), 
             ylim=better_lim(image$state$W1))
        legends = c("Growth rate") # legend label
        l_cols = c("yellow3") # legend color
        
        ## Add W''(t_new) curve
        par(new = TRUE) # new plot, same area
        plot(image$state$t_new, image$state$W2, type="l", col="navy", axes=FALSE, 
             xlab="", ylab="", ylim=better_lim(image$state$W2))
        axis(side = 4) # add acceleration y-axis on right side of the plot
        mtext(paste0("Acceleration (", dep, " / ", self$options$time, "²)"), 
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
