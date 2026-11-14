gcurveClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "gcurveClass",
    inherit = gcurveBase,
    private = list(
      
      .init = function() {
        
        op <- self$options
        results <- self$results
        
        hasDeps <- length(op$deps) > 0L
        hasTime <- length(op$time) > 0L
        
        setups <- pgm_resolve_modeling(
          modeling = op$modeling,
          deps = op$deps,
          validate = FALSE
        )
        
        results$documentation$setVisible(!(hasDeps && hasTime))
        
        for (setup in setups) {
          
          setupId <- setup$setupId
          label <- paste0(setup$var, ": ", setup$setupLabel)
          
          results$modelEval$gofTable$addRow(
            rowKey = setupId,
            values = list(var = label)
          )
          
          results$modelEval$emTable$addRow(
            rowKey = setupId,
            values = list(var = label)
          )
          
          results$aucTable$addRow(
            rowKey = setupId,
            values = list(var = label)
          )
          
          results$critPointsTable$addRow(
            rowKey = paste0(setupId, "_time"),
            values = list(var = label, Axis = if (hasTime) op$time else "")
          )
          
          results$critPointsTable$addRow(
            rowKey = paste0(setupId, "_dep"),
            values = list(Axis = setup$var)
          )
        }
      },
      
      .run = function() {
        
        ##### Data -----
        
        ## Option values into shorter variable names
        op <- self$options
        results <- self$results
        deps <- op$deps
        time <- op$time
        corrGroups <- op$corrGroups
        
        ## Get required data
        if (length(deps) == 0L || length(time) == 0L)
          return()
        
        if (isTRUE(op$aucInterval) && op$aucUpper <= op$aucLower) {
          jmvcore::reject(
            "The upper AUC limit must be greater than the lower limit"
          )
        }
        
        fixedTermBlocks <- pgm_parameter_terms(op$fixedTerms)
        randomGroupBlocks <- pgm_parameter_terms(op$randomGroups)
        
        fixedTerms <- unique(unlist(fixedTermBlocks, use.names = FALSE))
        randomGroups <- unique(unlist(randomGroupBlocks, use.names = FALSE))
        vars <- unique(c(deps, time, corrGroups, fixedTerms, randomGroups))
        data <- self$data[, vars, drop = FALSE]
        
        ## Convert time to appropriate data type
        data[[time]] <- jmvcore::toNumeric(data[[time]])
        
        
        ##### Curve Resolution -----
        
        ## Giving time a new "resolution" for smooth lines
        if (length(unique(na.omit(data[[time]]))) < 2)
          jmvcore::reject("At least two different time values are required")

        trangeGlobal <- range(data[[time]], na.rm = TRUE)
        tnewSize <- min(3000, 
          max(500, round(diff(trangeGlobal) * op$timeAxisFactor))
        )
        tnewGlobal <- seq(trangeGlobal[1], trangeGlobal[2], length.out = tnewSize)
        
        
        ##### Model Setups Iteration -----
        
        ## Global variables
        ySamplesList <- list()
        curvesList <- list(tnewGlobal = tnewGlobal)
        cpList <- list()
        aucs <- list()
        equations <- list()
        
        ## Modeling Setup Array
        setups <- pgm_resolve_modeling(modeling = op$modeling, deps = deps)
        
        for (setup in setups) {
          
          setupId <- setup$setupId
          setupLabel <- setup$setupLabel
          dep <- setup$var
          model <- setup$model
          evartype <- setup$errorVarType
          ecortype <- setup$errorCorType
          
          ## Build composite dep: Setup label for results output
          depSetupLabel <- paste0(dep,": ",setupLabel)
          
          ## Convert dep to appropriate data type
          data[[dep]] <- jmvcore::toNumeric(data[[dep]])  
          
          ## Remove NAs within variables required for this fit
          varsFit <- unique(c(time, dep, corrGroups, fixedTerms, randomGroups))
          datacln <- na.omit(data[, varsFit, drop = FALSE])
          
          ## Clean data vectors
          xraw <- datacln[[time]]
          yraw <- datacln[[dep]]
          
          ## Check if both variables have at least two unique values
          if (length(unique(xraw)) < 2 || length(unique(yraw)) < 2) {
            jmvcore::reject(
              "At least two unique (x, y) pairs are required for growth analysis"
            )
          }
          
          ##### Mean points -----
          samplesMean <- aggregate(
            yraw,
            by = list(xraw),
            FUN = mean,
            na.rm = TRUE
          )
          
          t <- samplesMean[[1]]
          y <- samplesMean[[2]]
          
          
          ##### Modeling -----
          
          ## Get model expression, parameters limits
          ## initial and transformed initial values
          mspec <- get_model_spec(
            model = model,
            t = t,
            y = y,
            parameterization = op$modelParam,
            double = op$modelPhase == "double",
            fixValue = setup$fixValue
          )
          mexpr  <- mspec$mexpr
          
          ## Fit data -----
          fitData <- data.frame(
            t = as.numeric(xraw),
            y = as.numeric(yraw)
          )
          
          ### Grouping variables
          groupVars <- unique(c(corrGroups, randomGroups))
          for (g in groupVars)
            fitData[[g]] <- as.factor(datacln[[g]])
          
          ### Fixed-effect predictors
          for (v in setdiff(fixedTerms, groupVars))
            fitData[[v]] <- datacln[[v]]
          
          ## Fit -----
          fittedmodel <- pgm_fit(
            op       = op,
            fitData  = fitData,
            mexpr    = mexpr,
            mformula = mspec$mformula,
            tparinit = mspec$tparinit,
            evartype = evartype,
            ecortype = ecortype
          )

          
          ##### Parameters and evaluation -----
          
          fittedpar <- fitted_parameters(
            fit = fittedmodel$fit,
            mspec = mspec,
            fixedTerms = fittedmodel$fixedTerms,
            data = fittedmodel$fitData
          )
          betaHat <- fittedpar$beta
          
          fittedrandom <- NULL
          if (inherits(fittedmodel$fit, "nlme")) {
            fittedrandom <- fitted_random_effects(
              fit = fittedmodel$fit,
              randomGroups = fittedmodel$randomGroups,
              randomGroupLabels = fittedmodel$randomGroupLabels,
              mspec = mspec,
              betaTrans = fittedpar$betaTrans
            )
          }
          
          fittedeval <- fittedmodel_evaluation(
            op = op,
            fit = fittedmodel$fit,
            npar = fittedpar$nfixed
          )
          
          ##### dep Resolution -----
          
          ## Trim tnewGlobal unit the max value of t.
          ### This allows graphical comparison in the same plot,
          ### even if measures stopped at a certain age for one variable
          tnew <- tnewGlobal[tnewGlobal >= min(t) & tnewGlobal <= max(t)]
          
          
          ##### Growth Function -----
          
          ## Get growth function and its derivatives
          Wderivs <- ur_derivative_functions(
            beta = betaHat,
            family = model,
            parameterization = op$modelParam,
            fixValue = setup$fixValue
          )

          
          ##### Critical Points by depSetupLabel -----
          cp <- compute_critical_points(
            op = op,
            t = tnew, 
            Wderivs = Wderivs,
            thValue = op$thValue
          )
          
          
          ##### AUC by setupId -----
          aucs[[setupId]] <- compute_auc(
            op      = op,
            t       = tnew, 
            Wderivs = Wderivs
          )
            
          
          ##### Results -----
          
          ## Estimation info
          if (op$estTables) {
            populate_etable(
              table = results$estTable,
              setup = setupId,
              label = depSetupLabel,
              est = fittedmodel$estInfo,
              opt = fittedmodel$optInfo,
              fixedInfo = fittedmodel$fixedInfo,
              asymptote = setup$fixValue,
              randomInfo = fittedmodel$randomInfo,
              weq = fittedmodel$evarEq,
              ceq = fittedmodel$ecorEq
            )
          }
          
          # Parameters info
          if (op$parTables) {
            note <- NULL
            if (any(lengths(fittedmodel$fixedTerms) > 0L))
              note <- paste0(
                "p-values are obtained on model estimation scale; ",
                "factor levels test contrast against reference level."
              )
            
            populate_ptable(
              table = results$parTable,
              setup = setupId,
              label = depSetupLabel,
              pardf = fittedpar$pardf,
              note = note
            )
          }
          
          if (op$parTables && !is.null(fittedrandom) && nrow(fittedrandom) > 0L) {
            results$randomTable$setVisible(TRUE)
            populate_rtable(
              table = results$randomTable,
              setup = setupId,
              label = depSetupLabel,
              randomdf = fittedrandom
            )
          }
          
          # Goodness of Fit metrics
          if (any(!is.na(unlist(fittedeval$gof)))) {
            populate_goftable(
              table   = results$modelEval$gofTable,
              setup   = setupId,
              label   = depSetupLabel,
              metrics = fittedeval$gof
            )
          }
          
          # Error metrics
          if (any(!is.na(unlist(fittedeval$em)))) {
            populate_emtable(
              table = results$modelEval$emTable,
              setup   = setupId,
              label   = depSetupLabel,
              metrics = fittedeval$em
            )
          }
          
          # Area Under the Curve metrics
          if (length(op$aucFunction) > 0L) {
            populate_auctable(
              table    = results$aucTable,
              setup    = setupId,
              label    = depSetupLabel,
              interval = aucs[[setupId]]$interval,
              metrics  = aucs[[setupId]]$metrics,
              note     = aucs[[setupId]]$note
            )
          }
          
          # Critical points location
          if (any(!is.na(unlist(cp)))) {
            populate_cptable(
              table     = results$critPointsTable,
              setup     = setupId,
              label     = depSetupLabel,
              taxis     = op$time,
              daxis     = dep,
              critPoints = cp, 
              dfunc     = Wderivs$W, 
              note      = if (!is.null(cp$note)) cp$note
            )
          }
          
          
          ##### Render Data -----
          
          ## Growth Curve plot
          if (op$growthPlots) {
            growthPlot <- results$growthPlot$addItem(key = setupId)
            
            growthPlot$setTitle(paste("Growth Curve for", depSetupLabel))
            
            growthPlot$setState(
              populate_growthplot(
                op = op,
                tCurves = tnew,
                xSamples = xraw,
                ySamples = yraw,
                functions = Wderivs,
                cp = cp,
                dep = dep, time = time,
                yLegend = depSetupLabel
              )
            )
          }
          
          ## Residuals plot
          if (op$residualPlots) {
            
            residualPlot <- results$residualPlot$addItem(
              key = setupId
            )
            
            residualPlot$setTitle(
              paste("Residual Diagnostics for", depSetupLabel)
            )
            
            residualPlot$setState(
              list(
                timeLabel = time,
                diagnostics = residual_diagnostics(
                  fit = fittedmodel$fit,
                  time = xraw,
                  id = prepare_diagnostic_id(
                    op = op,
                    groups = corrGroups,
                    data = datacln,
                    nlmegroups = if (inherits(fittedmodel$fit, "nlme")) {
                      fittedmodel$fit$groups
                    } else {NULL}
                  )
                )
              )
            )
          }
          
          ## Growth comparison plot
          comparisonData <- prepare_comparison_data(
            op = op,
            tCurves = tnew,
            functions = Wderivs,
            cp = cp,
            dep = dep,
            time = time
          )
          ySamplesList[[depSetupLabel]] <- list(x = xraw, y = yraw)
          curvesList[[depSetupLabel]] <- comparisonData$curves
          cpList[[depSetupLabel]] <- comparisonData$cps
          
        } # close loop
        
        
        ##### Global Info -----
        
        ## Rank AUC table
        if (isTRUE(results$aucTable$visible)) {
          rank_auctable(
            table = results$aucTable,
            aucs = aucs
          )
        }
        
        ## Multiple models plot
        if (op$growthComparisonPlots) {
          growthComparisonPlot <- results$growthComparisonPlot
          
          depNames <- unique(vapply(
            setups, function(setup) setup$var, character(1)
          ))
          depSetupLabels <- setdiff(names(curvesList), "tnewGlobal")
          growthAxisLabel <- if (length(depNames) == 1L)
            depNames[[1L]]
          else
            "Growth"
          rightAxisLabel <- NULL
          
          if (op$comparisonRightFunc != "none" && length(depSetupLabels) > 0L) {
            firstCurve <- curvesList[[depSetupLabels[[1L]]]]
            rightAxisLabel <- if (length(depNames) == 1L)
              firstCurve$rightAxisLabel
            else
              firstCurve$rightAxisLegend
          }
          
          growthComparisonPlot$setState(list(
            curvesList = curvesList,
            ySamplesList = ySamplesList,
            cpList = cpList,
            timeLabel = time,
            growthAxisLabel = growthAxisLabel,
            rightAxisLabel = rightAxisLabel
          ))
        }
        
      }, # close .run
      
      ##### Render Functions -----
      
      # Residuals plot
      .residualPlot = function(image, ...) {
        
        if (is.null(image$state))
          return(FALSE)
        
        diag <- image$state$diagnostics
        
        fittedValues <- diag$fittedValues
        residualNorm <- diag$residualNorm
        timeValues <- diag$time
        
        valid <- is.finite(fittedValues) &
          is.finite(residualNorm) &
          is.finite(timeValues)
        
        fittedValues <- fittedValues[valid]
        residualNorm <- residualNorm[valid]
        timeValues <- timeValues[valid]
        
        if (length(residualNorm) < 3L)
          return(FALSE)
        
        oldPlotPar <- par(
          mfrow = c(2, 2),
          mar = c(4.2, 4.2, 3.0, 1.0),
          mgp = c(2.4, 0.7, 0)
        )
        
        on.exit(par(oldPlotPar), add = TRUE)
        
        ## 1. Normalized residuals versus fitted
        plot(fittedValues, residualNorm, pch = 16, cex = 0.7,
          xlab = "Fitted values", ylab = "Normalized residuals",
          main = "Residuals vs fitted"
        )
        
        abline(h = 0, lty = 2)
        
        if (length(unique(fittedValues)) >= 4L) {
          smooth <- stats::lowess(fittedValues, residualNorm)
          lines(smooth, lwd = 2)
        }
        
        ## 2. Normalized residuals versus time
        plot(timeValues, residualNorm, pch = 16, cex = 0.7,
          xlab = image$state$timeLabel, ylab = "Normalized residuals",
          main = "Residuals vs time")
        
        abline(h = 0, lty = 2)
        
        if (length(unique(timeValues)) >= 4L) {
          smooth <- stats::lowess(timeValues, residualNorm)
          lines(smooth, lwd = 2)
        }
        
        ## 3. Normal Q-Q
        qqnorm(residualNorm, pch = 16, cex = 0.7,
          main = "Normal Q–Q", ylab = "Normalized residuals"
        )
        
        qqline(residualNorm, lty = 2)
        
        ## 4. ACF
        if ("id" %in% names(diag)) {
          
          acfData <- pooled_subject_acf(
            residuals = diag$residualNorm,
            time = diag$time,
            id = diag$id
          )
          
          if (!is.null(acfData)) {
            
            plot(
              acfData$lag, acfData$acf, type = "h", lwd = 2,
              xlab = "Within-group lag", ylab = "Mean residual ACF",
              main = "Within-group residual ACF", ylim = c(-1, 1)
            )
            abline(h = 0)
            
          } else {
            
            plot.new()
            text(0.5, 0.5,
              "Insufficient repeated observations\nfor within-group ACF"
            )
          }
          
        } else {
          
          stats::acf(residualNorm, main = "Residual ACF")
        }
        
        
        ## Notify the rendering system that we have plotted something
        TRUE
      }, # close .residualPlot
      
      
      # Model plot
      .growthPlot = function(image, ...) {
        
        ## Check if there is any data to plot
        if (is.null(image$state))
          return(FALSE)
        
        
        ## Options
        op <- self$options
        
        ## Define data variables
        samples <- image$state$samples
        curves <- image$state$curves
        cps <- image$state$cps
        
        
        ## Deactivate the default box around plotting area
        par(bty = 'L')
        
        ## Redefine plot margins to fit axis labels
        par(mar = c(5.1, 3.1, 4.1, 3.1))  # c(5.1, 4.1, 4.1, 2.1)
        par(mgp = c(2, 0.5, 0)) # c(3, 1, 0)
        
        ## Create plot with samples points
        samplesColor = if ("dataPoints" %in% op$plotsElements) "grey" else NA
        
        plot(
          samples$x, samples$y, 
          xlab = image$state$timeLabel, ylab = curves$growth$axisLabel, 
          col = samplesColor, pch = 16, ylim = better_lim(samples$y)
        )
        
        ### Add W(t) curve
        lines(curves$growth$x, curves$growth$y, 
              col = op$growthCurveColor, lwd = 2)
        
        legendLabels = c(curves$growth$legend) # legend label
        legendColors = c(op$growthCurveColor) # legend color
        legendTypes = c("solid") # legend line type
        
        ## Add Key Growth Points
        if (!is.null(cps)) {
          for (cpName in names(cps)) {
            cp <- cps[[cpName]]
            draw_cpoint(cp$x, cp$y, label = cpName)
          }
        }

        
        ## Add right y-axis curve
        if (op$growthRightAxis) {
          
          draw_new_curve(
            x = curves$rightAxis$x,
            y = curves$rightAxis$y,
            axisLabel = curves$rightAxis$axisLabel,
            linetype = op$growthRightAxisLine,
            linecolor = op$growthRightAxisColor
          )
          
          legendLabels = c(legendLabels, curves$rightAxis$legend) # legend label
          legendColors = c(legendColors, op$growthRightAxisColor) # legend color
          legendTypes = c(legendTypes, op$growthRightAxisLine) # legend line type
        }
        
        
        ## Add a third curve (just for the shape over time)
        if (op$growthNoAxis && !(
          op$growthRightAxis && op$growthNoAxisFunc == op$growthRightAxisFunc)
          ) 
        {
          
          draw_new_curve(
            x = curves$noAxis$x,
            y = curves$noAxis$y,
            linetype = op$growthNoAxisLine,
            linecolor = op$growthNoAxisColor
          )
          
          legendLabels = c(legendLabels, curves$noAxis$legend) # legend label
          legendColors = c(legendColors, op$growthNoAxisColor) # legend color
          legendTypes = c(legendTypes, op$growthNoAxisLine) # legend line type
        }
        
        
        ## Add lines legend
        legend(
          "bottom", inset = c(0, -0.2), legend = legendLabels, 
          col = legendColors, lty = legendTypes, horiz = TRUE, 
          bty = "n", lwd = 2, xpd = TRUE
        )
        
        
        ## Notify the rendering system that we have plotted something
        TRUE
        
      }, # close .growthPlot
      
      # Multiple plot
      .growthComparisonPlot = function(image, ...) {
        
        ## Check if there is any data to plot
        if (is.null(image$state))
          return(FALSE)
        
        ## Options
        op <- self$options
        
        ## Data
        curvesList <- image$state$curvesList
        ySamplesList <- image$state$ySamplesList
        
        ## Setup labels and global time grid
        tnewGlobal <- curvesList$tnewGlobal
        depSetupLabels <- setdiff(names(curvesList), "tnewGlobal")
        if (length(depSetupLabels) == 0L)
          return(FALSE)
        
        ## Generate colors in the standard jamovi palette
        jmvColors <- jmvcore::colorPalette(n = length(depSetupLabels))
        
        ## Global vertical limits
        ylimTotal <- range(
          unlist(lapply(ySamplesList[depSetupLabels], \(s) s$y)),
          na.rm = TRUE
        )
        ylimTotal <- better_lim(ylimTotal)
        
        ## Global horizontal limits
        xlimTotal <- range(tnewGlobal)
        
        ## Plot settings
        par(bty = 'L')
        
        ## Redefine plot margins to fit axis labels
        par(mar = c(5.1, 3.1, 4.1, 3.1))  # c(5.1, 4.1, 4.1, 2.1)
        par(mgp = c(2, 0.5, 0)) # c(3, 1, 0)
        
        ## First setup
        samples <- ySamplesList[[depSetupLabels[1]]]
        
        samplesColor <- NA
        if ("dataPoints" %in% op$plotsElements)
          samplesColor <- adjustcolor(jmvColors[1], alpha.f = 0.3)
        
        plot(
          samples$x, samples$y, pch = 16, xlim = xlimTotal, ylim = ylimTotal,
          xlab = image$state$timeLabel, ylab = image$state$growthAxisLabel,
          col = samplesColor
        )
        
        ## Remaining setups
        if ("dataPoints" %in% op$plotsElements &&
            length(depSetupLabels) > 1L) {
          
          for (i in seq_along(depSetupLabels)[-1]) {
            
            samples <- ySamplesList[[depSetupLabels[i]]]
            points(samples$x, samples$y, pch = 16,
                   col = adjustcolor(jmvColors[i], alpha.f = 0.3))
          }
        }
        
        ## Add growth curves
        for (i in seq_along(depSetupLabels)) {
          
          depSetupLabel <- depSetupLabels[i]
          x <- curvesList[[depSetupLabel]]$x
          ypred <- curvesList[[depSetupLabel]]$growth
          
          ### W(t) curve
          lines(x, ypred, col = jmvColors[i], lwd = 2)
          
          ## Add Key Growth Points
          cps <- image$state$cpList[[depSetupLabel]]
          
          if (!is.null(cps)) {
            for (cpName in names(cps)) {
              cp <- cps[[cpName]]
              draw_cpoint(cp$x, cp$y, label = cpName, 
                pcolor = jmvColors[i], scolor = jmvColors[i]
              )
            }
          }
          
        }
        
        ## Make legend
        legend(
          "bottom", inset = c(0, -0.2), legend = depSetupLabels, bty = "n",
          col = jmvColors, lty = 1, lwd = 2, horiz = TRUE, xpd = TRUE
        )
        
        ## Right Axis
        if (op$comparisonRightFunc != "none") {
          
          ylimRight <- better_lim(range(
            unlist(lapply(curvesList[depSetupLabels], \(item) item$rightAxis)),
            na.rm = TRUE
          ))
          
          par(new = TRUE)
          
          plot(tnewGlobal, rep(NA_real_, length(tnewGlobal)), type="n",
               axes=FALSE, xlab="", ylab="", xlim=xlimTotal, ylim=ylimRight)
          
          for (i in seq_along(depSetupLabels)) {
            depSetupLabel <- depSetupLabels[i]
            lines(curvesList[[depSetupLabel]]$x, 
                  curvesList[[depSetupLabel]]$rightAxis, 
                  col = jmvColors[i], lty = "dashed")
          }
          
          axis(4)
          mtext(image$state$rightAxisLabel, side = 4, line = 2)
        }
        
        ## Notify the rendering system that we have plotted something
        TRUE
        
      } # close .growthComparisonPlot
      
    ) # close list
) # close R6Class
