gcurveClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
    "gcurveClass",
    inherit = gcurveBase,
    private = list(
      
      .run = function() {
        
        ##### Data -----
        
        ## Option values into shorter variable names
        op <- self$options
        results <- self$results
        deps  <- op$deps
        time <- op$time
        iid <- op$iid
        
        ## Check if variables have any data and get the data
        if (!is.null(deps) & !is.null(time)) {
          if (!is.null(iid)) {
            data <- self$data[, c(deps, time, iid), drop=FALSE]
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
        if (length(unique(na.omit(data[[time]]))) < 2)
          stop("At least two different time values are required")
        trawGlobal <- data[[time]]
        tnewSize <- min(
          3000, 
          max(
            500, 
            round((max(trawGlobal) - min(trawGlobal)) * op$timeAxisFactor)
          )
        )
        tnewGlobal <- seq(min(trawGlobal), max(trawGlobal), length.out=tnewSize)
        
        
        ##### Model Setups Iteration -----
        
        ## Global variables
        ySamplesList <- list(trawGlobal = trawGlobal)
        curvesList <- list(tnewGlobal = tnewGlobal)
        cpList <- list()
        aucs <- list()
        equations <- list()
        
        ## Modeling Setup Array
        setups <- op$modeling
        
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
          
          ## Remove NAs
          varsFit <- c(time, dep)
          if (!is.null(iid))
            varsFit <- c(varsFit, iid)
          datacln <- na.omit(data[, varsFit, drop = FALSE])
          
          ## Clean data vectors
          xraw <- datacln[[time]]
          yraw <- datacln[[dep]]
          
          ## Check if both variables have at least two unique values
          if (length(unique(xraw)) < 2 || length(unique(yraw)) < 2)
            stop("At least two unique (x, y) pairs are required for growth analysis")
  
          
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
          mspec <- get_model_spec(model, t, y)
          mexpr  <- mspec$mexpr
          
          ## Store model name and equation
          if (!mspec$mname %in% equations){
            equations[[mspec$mname]] <- mspec$mhtml
          }
          
          ## Fit
          fittedmodel <- pgm_fit(
            op       = op,
            t        = xraw,
            y        = yraw,
            iid      = if (!is.null(iid)) datacln[[iid]] else NULL,
            mexpr    = mexpr,
            mformula = mspec$mformula,
            parinit  = mspec$parinit,
            tparinit  = mspec$tparinit,
            evartype = evartype,
            ecortype = ecortype
          )

          
          ##### Parameters and evaluation -----
          
          fittedpar <- fitted_parameters(
            fit  = fittedmodel$fit, 
            mspec = mspec
          )
          
          betaHat <- setNames(
            fittedpar$pardf[["parvalue"]], rownames(fittedpar$pardf)
          )
          
          fittedeval <- fittedmodel_evaluation(
            op   = op,
            y    = yraw,
            fit  = fittedmodel$fit,
            npar = nrow(fittedpar$pardf)
          )
          
          ##### dep Resolution -----
          
          ## Trim tnewGlobal unit the max value of t.
          ### This allows graphical comparison in the same plot,
          ### even if measures stopped at a certain age for one variable
          tnew <- tnewGlobal[tnewGlobal >= min(t) & tnewGlobal <= max(t)]
          
          
          ##### Growth Function -----
          
          ## Model equation with fitted parameters
          Wexpr <- do.call("substitute", list(
            mexpr[[1]], as.list(betaHat)
          ))
          
          ## Get growth function and its derivatives functions
          Wderivs <- get_W_derivatives(op, Wexpr, t)
          
          ## Apply growth function to tnew
          Wpred <- Wderivs$W(tnew)

          
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
              est   = fittedmodel$estInfo,
              opt   = fittedmodel$optInfo,
              randomEff = fittedmodel$randomEff,
              weq   = fittedmodel$evarEq,
              ceq   = fittedmodel$ecorEq
            )
          }
          
          # Parameters info
          if (op$parTables) {
            populate_ptable(
              table = results$parTable,
              setup = setupId,
              label = depSetupLabel,
              pardf = fittedpar$pardf
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
          if (any(!is.na(unlist(aucs[[setupId]]$metrics)))) {
            populate_auctable(
              table    = results$aucTable,
              setup    = setupId,
              label    = depSetupLabel,
              interval = aucs[[setupId]]$interval,
              metrics  = aucs[[setupId]]$metrics
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
              paste("Residual diagnostics for", depSetupLabel)
            )
            
            residualPlot$setState(
              list(
                timeLabel = time,
                diagnostics = residual_diagnostics(
                  fit = fittedmodel$fit,
                  time = xraw,
                  id = if (!is.null(iid)) datacln[[iid]] else NULL
                )
              )
            )
          }
          
          ## Growth comparison plot
          comparisonData <- prepare_comparison_data(
            op = op,
            tCurves = tnew,
            tCurvesGlobal = tnewGlobal,
            tSamplesGlobal = trawGlobal,
            ySamples = yraw,
            functions = Wderivs,
            cp = cp
          )
          ySamplesList[[depSetupLabel]]  <- comparisonData$ySamples
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
        
        ## Print used models equations
        if (op$eqHtmls) {
          populate_eqhtml(results, equations)
        }
        
        ## Multiple models plot
        if (op$growthComparisonPlots) {
          
          growthComparisonPlot <- results$growthComparisonPlot
          
          growthComparisonPlot$setState(list(
            curvesList = curvesList,
            ySamplesList = ySamplesList,
            cpList = cpList
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
            
            plot(acfData$lag, acfData$acf, type = "h", lwd = 2,
              xlab = "Within-subject lag", ylab = "Mean residual ACF",
              main = "Within-subject residual ACF", ylim = c(-1, 1)
            )
            
            abline(h = 0)
          }
          
        } else {
          
          stats::acf(
            residualNorm,
            main = "Residual ACF"
          )
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
        
        plot(samples$x, samples$y, xlab = "", ylab = curves$growth$axisLabel,
             col = samplesColor, pch=16, ylim = better_lim(samples$y))
        
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
        if (op$growthNoAxis & op$growthNoAxisFunc != op$growthRightAxisFunc) {
          
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
        legend("bottom", inset = -0.2, legend = legendLabels, col = legendColors,
               lty = legendTypes, horiz = TRUE, bty = "n", lwd = 2, xpd = TRUE)
        
        
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
        
        ## Time columns and y columns names
        tnewGlobal <- curvesList$tnewGlobal
        trawGlobal <- ySamplesList$trawGlobal
        depSetupLabels <- setdiff(names(curvesList), "tnewGlobal")

        ## Generate colors in the standard jamovi palette
        jmvColors <- jmvcore::colorPalette(n = length(depSetupLabels))
        
        ## Join all samples y to calculate limits
        ylimTotal <- range(unlist(ySamplesList[depSetupLabels]), na.rm = TRUE)
        ylimTotal <- better_lim(ylimTotal)
        
        ## Deactivate the default box around the plotting area
        par(bty='L')
        
        ## Create and set-up the plot
        samplesColor <- NA
        
        if ("dataPoints" %in% op$plotsElements)
          samplesColor <- adjustcolor(jmvColors[1], alpha.f = 0.3)
        
        plot(trawGlobal, ySamplesList[[depSetupLabels[1]]], pch=16, 
             xlab="", ylab="", ylim = ylimTotal, col = samplesColor)
        
        ## Add the other variables samples
        if (length(depSetupLabels) > 1) {
          
          for (i in seq_along(depSetupLabels)[-1]) {
            
            depSetupLabel <- depSetupLabels[i]
            samplesColor <- NA
            
            if ("dataPoints" %in% op$plotsElements){
              
              samplesColor <- adjustcolor(jmvColors[i], alpha.f = 0.3)
              depSetupLabel <- depSetupLabels[i]
              samples <- ySamplesList[[depSetupLabel]]
              
              points(trawGlobal, samples, pch = 16, col = samplesColor)
            }
            
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
        legend("bottom", inset = -0.2, legend = depSetupLabels, 
               col = jmvColors, lty = 1, lwd = 2, horiz = TRUE, 
               bty = "n", xpd = TRUE)
          
        
        ## Right Axis
        if (op$comparisonRightFunc != "none") {
          
          x <- curvesList[[depSetupLabel]]$x
          xlimRight <- c(min(trawGlobal), max(trawGlobal))
          
          ylimRight <- better_lim(range(unlist(
            lapply(curvesList[depSetupLabels], \(item) item$rightAxis)
          ), na.rm = TRUE))

          par(new = TRUE)
          
          plot(x, curvesList[[depSetupLabels[1]]]$rightAxis, type = "n", 
               axes = FALSE, xlab = "", ylab = "", ylim = ylimRight)
          
          for (i in seq_along(curvesList))
            lines(x, curvesList[[depSetupLabels[i]]]$rightAxis, 
                  col = jmvColors[i], lty = "dashed")
          
          axis(4)
        }
        
        
        ## Notify the rendering system that we have plotted something
        TRUE
        
      } # close .growthComparisonPlot
      
    ) # close list
) # close R6Class
