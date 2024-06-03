# plotCalibration2() will not run unless we assign this function
# from the riskRegression library
getLegendData <- function(object,
                          models,
                          times,
                          brier.in.legend=TRUE,
                          format.brier,
                          auc.in.legend=TRUE,
                          format.auc,
                          drop.null.model=TRUE,
                          scale=100,
                          digits=1,
                          ...){
  model=AUC=lower=upper=Brier=NULL
  if (missing(models)) {
    models <- names(object$models)
  }
  if(!is.null(object$null.model)){
    if (drop.null.model==TRUE) models <- models[models!=object$null.model]
  }
  maxlen <- max(nchar(as.character(models)))
  legend.text.models <- sprintf(paste0("%",maxlen,"s"),models)
  if (missing(format.auc))
    format.auc <- paste0("%1.",digits,"f [%1.",digits,"f;%1.",digits,"f]")
  if (missing(format.brier))
    format.brier <- paste0("%1.",digits,"f [%1.",digits,"f;%1.",digits,"f]")
  if (is.null(object$null.model) || drop.null.model[[1]]==TRUE){
    keep.null.model <- FALSE
  }else{
    if (brier.in.legend==TRUE){
      keep.null.model=TRUE
    }else{
      keep.null.model=match(object$null.model,models,nomatch=FALSE)
    }
  }
  if (auc.in.legend==TRUE){
    if (is.null(object$AUC)){
      warning("Cannot show AUC as it is not stored in object. Set metrics='auc' in the call of Score.")
      legend.text.auc <- NULL
    }else{
      auc.data <- object$AUC$score[(model%in%models)]
      if (object$response.type!="binary"){
        if (missing(times)){
          tp <- max(auc.data[["times"]])
          if (length(unique(auc.data$times))>1)
            warning("Time point not specified, use max of the available times: ",tp)
        } else{ ## can only do one time point
          tp <- times[[1]]
          if (!(tp%in%unique(auc.data$times)))
            stop(paste0("Requested time ",times[[1]]," is not in object"))
        }
        auc.data <- auc.data[times==tp]
      }else tp <- NULL
      if (keep.null.model==FALSE){
        if(!is.null(object$null.model)){
          auc.data <- auc.data[model!=object$null.model]
        }
      }
      ## user's order
      auc.data[,model:=factor(model,levels=models)]
      setkey(auc.data,model)
      legend.text.auc <- auc.data[,sprintf(fmt=format.auc,scale*AUC,scale*lower,scale*upper)]
    }
  }else{
    legend.text.auc <- NULL
  }
  if (brier.in.legend==TRUE){
    if (is.null(object$Brier)){
      warning("Cannot show Brier score as it is not stored in object. Set metrics='brier' in the call of Score.")
      legend.text.brier <- NULL
    }else{
      brier.data <- object$Brier$score[(model%in%models)]
      if (object$response.type!="binary"){
        if (missing(times)){
          tp <- max(brier.data[["times"]])
          if (length(unique(brier.data$times))>1)
            warning("Time point not specified, use max of the available times: ",tp)
        } else{ ## can only do one time point
          tp <- times[[1]]
          if (!(tp%in%unique(brier.data$times)))
            stop(paste0("Requested time ",times[[1]]," is not in object"))
        }
        brier.data <- brier.data[times==tp]
      }else tp <- NULL
      if (!is.null(object$null.model) && keep.null.model[[1]]==FALSE){
        brier.data <- brier.data[model!=object$null.model]
      }
      brier.data[,model:=factor(model,levels=models)]
      setkey(brier.data,model)
      legend.text.brier <- brier.data[,sprintf(fmt=format.brier,scale*Brier,scale*lower,scale*upper)]
    }
  }else{
    legend.text.brier <- NULL
  }
  out <- cbind(legend.text.models,
               AUC=legend.text.auc,
               Brier=legend.text.brier)
  out
}


plotCalibration2 <- function (x, models, times, method = "nne", cens.method, 
          round = TRUE, bandwidth = NULL, q = 10, bars = FALSE, hanging = FALSE, 
          names = "quantiles", pseudo = FALSE, rug, show.frequencies = FALSE, 
          plot = TRUE, add = FALSE, diag = !add, legend = !add, auc.in.legend, 
          brier.in.legend, axes = !add, xlim = c(0, 1), ylim = c(0, 
                                                                 1), xlab = ifelse(bars, "Risk groups", "Predicted risk"), 
          ylab, col, lwd, lty, pch, type, percent = TRUE, na.action = na.fail, 
          cex = 1, ...) 
{
  require(data.table)
  require(riskRegression)
  
  if (x$response.type != "binary" && missing(cens.method)) {
    cens.method <- "local"
    message("The default method for estimating calibration curves based on censored data has changed for riskRegression version 2019-9-8 or higher\nSet cens.method=\"jackknife\" to get the estimate using pseudo-values.\nHowever, note that the option \"jackknife\" is sensitive to violations of the assumption that the censoring is independent of both the event times and the covariates.\nSet cens.method=\"local\" to suppress this message.")
  }
  model = risk = event = status = NULL
  if (missing(ylab)) {
    if (bars == TRUE) {
      ylab = ""
    }
    else {
      if (x$cens.type == "rightCensored") {
        ylab = "Estimated actual risk"
      }
      else {
        ylab = "Observed frequency"
      }
    }
  }
  if (missing(auc.in.legend)) 
    auc.in.legend <- ("auc" %in% tolower(x$metrics))
  if (missing(brier.in.legend)) 
    brier.in.legend <- ("auc" %in% tolower(x$metrics))
  if (missing(pseudo) & missing(rug)) 
    if (x$cens.type == "rightCensored") {
      showPseudo <- FALSE
      showRug <- FALSE
    }
  else {
    showPseudo <- FALSE
    showRug <- TRUE
  }
  else {
    if (missing(pseudo) || pseudo[[1]] == FALSE) 
      showPseudo <- FALSE
    else showPseudo <- TRUE
    if (missing(rug) || rug[[1]] == FALSE) 
      showRug <- FALSE
    else showRug <- TRUE
  }
  pframe <- x$Calibration$plotframe
  if (is.null(pframe)) 
    stop("Object has no information for calibration plot.\nYou should call the function \"riskRegression::Score\" with plots=\"calibration\".")
  Rvar <- grep("^(ReSpOnSe|pseudovalue)$", names(pframe), 
               value = TRUE)
  if (!missing(models)) {
    fitted.models <- pframe[, unique(model)]
    if (all(models %in% fitted.models)) {
      pframe <- pframe[model %in% models]
    }
    else {
      if (all(is.numeric(models)) && (max(models) <= length(fitted.models))) {
        models <- fitted.models[models]
        pframe <- pframe[model %in% models]
      }
      else {
        stop(paste0("The requested models: ", models, 
                    "\ndo not all match the fitted models: ", 
                    paste0(fitted.models, collapse = ", ")))
      }
    }
  }
  data.table::setkey(pframe, model)
  if (x$response.type != "binary") {
    if (missing(times)) {
      tp <- max(pframe[["times"]])
      if (length(unique(pframe$times)) > 1) 
        warning("Time point not specified, use max of the available times: ", 
                tp)
    }
    else {
      tp <- times[[1]]
      if (!(tp %in% unique(pframe$times))) 
        stop(paste0("Requested time ", times[[1]], 
                    " is not in object"))
    }
    pframe <- pframe[times == tp]
  }
  else tp <- NULL
  NF <- pframe[, length(unique(model))]
  if (bars) {
    method = "quantile"
    if (!(NF == 1)) 
      stop(paste0("Barplots work only for one risk prediction model at a time. Provided are ", 
                  NF, "models."))
  }
  if (missing(lwd)) 
    lwd <- rep(3, NF)
  if (missing(col)) {
    if (bars) 
      col <- c("grey90", "grey30")
    else {
      cbbPalette <- c("#000000", "#E69F00", 
                      "#56B4E9", "#009E73", "#D55E00", 
                      "#0072B2", "#CC79A7", "#F0E442")
      if (NF > length(cbbPalette)) 
        col <- 1:NF
      else col <- cbbPalette[1:NF]
    }
  }
  if (missing(type)) {
    if (method == "quantile") {
      type <- rep("b", length.out = NF)
    }
    else {
      type <- rep("l", length.out = NF)
    }
  }
  if (missing(lty)) 
    lty <- rep(1, length.out = NF)
  if (missing(pch)) 
    pch <- rep(1, length.out = NF)
  if (length(lwd) < NF) 
    lwd <- rep(lwd, length.out = NF)
  if (length(type) < NF) 
    type <- rep(type, length.out = NF)
  if (length(lty) < NF) 
    lty <- rep(lty, length.out = NF)
  if (length(col) < NF) 
    col <- rep(col, length.out = NF)
  if (length(pch) < NF) 
    pch <- rep(pch, length.out = NF)
  modelnames <- pframe[, unique(model)]
  axis1.DefaultArgs <- list(side = 1, las = 1, at = seq(0, 
                                                        xlim[2], xlim[2]/4))
  axis2.DefaultArgs <- list(side = 2, las = 2, at = seq(0, 
                                                        ylim[2], ylim[2]/4), mgp = c(4, 1, 0))
  if (is.character(legend[[1]]) || legend[[1]] == TRUE || auc.in.legend == 
      TRUE || brier.in.legend == TRUE) {
    legend.data <- getLegendData(object = x, models = models, 
                                 times = tp, auc.in.legend = auc.in.legend, brier.in.legend = brier.in.legend, 
                                 drop.null.model = TRUE)
    if (is.character(legend)) 
      legend.text <- legend
    else {
      legend.text <- unlist(legend.data[, 1])
    }
    nrows.legend <- NROW(legend.data)
    if (nrows.legend == 1) {
      legend.lwd <- NA
    }
    else {
      legend.lwd <- lwd
    }
    legend.DefaultArgs <- list(legend = legend.text, lwd = legend.lwd, 
                               col = col, ncol = 1, lty = lty, cex = cex, bty = "n", 
                               y.intersp = 1, x = "topleft", title = "")
    if (NCOL(legend.data) > 1) {
      addtable2plot.DefaultArgs <- list(yjust = 1.18, cex = cex, 
                                        table = legend.data[, -1, drop = FALSE])
    }
    else {
      addtable2plot.DefaultArgs <- NULL
    }
  }
  else {
    legend.DefaultArgs <- NULL
    addtable2plot.DefaultArgs <- NULL
  }
  if (bars) {
    legend.DefaultArgs <- list(legend = modelnames, col = col, 
                               cex = cex, bty = "n", x = "topleft")
    names.DefaultArgs <- list(cex = 0.7 * par()$cex, y = c(-abs(diff(ylim))/15, 
                                                           -abs(diff(ylim))/25))
    frequencies.DefaultArgs <- list(cex = 0.7 * par()$cex, 
                                    percent = FALSE, offset = 0)
  }
  else {
    if (length(modelnames) <= 1 && sum(auc.in.legend + brier.in.legend) == 
        0) {
      legend = FALSE
    }
  }
  if (bars) {
    if (x$cens.type == "rightCensored") 
      legend.DefaultArgs$legend <- c("Predicted risk", 
                                     "Estimated actual risk")
    else legend.DefaultArgs$legend <- c("Predicted risk", 
                                        "Observed frequency")
  }
  lines.DefaultArgs <- list(pch = pch, type = type, cex = cex, 
                            lwd = lwd, col = col, lty = lty)
  abline.DefaultArgs <- list(lwd = 1, col = "red")
  if (missing(ylim)) {
    if (showPseudo[1] && !bars[1]) {
      ylim <- range(pframe[[Rvar]])
    }
    else ylim <- c(0, 1)
  }
  if (missing(xlim)) {
    xlim <- c(0, 1)
  }
  plot.DefaultArgs <- list(x = 0, y = 0, type = "n", 
                           ylim = ylim, xlim = xlim, ylab = ylab, xlab = xlab)
  rug.DefaultArgs <- list(col = col, lwd = lwd, side = 1, ticksize = 0.03)
  pseudo.DefaultArgs <- list(col = col, cex = cex, density = 55)
  if (bars) {
    barplot.DefaultArgs <- list(ylim = ylim, col = col, axes = FALSE, 
                                ylab = ylab, xlab = xlab, beside = TRUE, legend.text = NULL, 
                                cex.axis = cex, cex.lab = par()$cex.lab, cex.names = cex)
    control <- prodlim::SmartControl(call = list(...), keys = c("barplot", 
                                                                "legend", "addtable2plot", "axis2", 
                                                                "abline", "names", "frequencies"), 
                                     ignore = NULL, ignore.case = TRUE, defaults = list(barplot = barplot.DefaultArgs, 
                                                                                        addtable2plot = addtable2plot.DefaultArgs, abline = abline.DefaultArgs, 
                                                                                        legend = legend.DefaultArgs, names = names.DefaultArgs, 
                                                                                        frequencies = frequencies.DefaultArgs, axis2 = axis2.DefaultArgs), 
                                     forced = list(abline = list(h = 0)), verbose = TRUE)
  }
  else {
    control <- prodlim::SmartControl(call = list(...), keys = c("plot", 
                                                                "rug", "pseudo", "lines", "legend", 
                                                                "addtable2plot", "axis1", "axis2"), 
                                     ignore = NULL, ignore.case = TRUE, defaults = list(plot = plot.DefaultArgs, 
                                                                                        pseudo = pseudo.DefaultArgs, rug = rug.DefaultArgs, 
                                                                                        addtable2plot = addtable2plot.DefaultArgs, lines = lines.DefaultArgs, 
                                                                                        legend = legend.DefaultArgs, axis1 = axis1.DefaultArgs, 
                                                                                        axis2 = axis2.DefaultArgs), forced = list(plot = list(axes = FALSE), 
                                                                                                                                  axis1 = list(side = 1)), verbose = TRUE)
  }
  method <- match.arg(method, c("quantile", "nne"))
  getXY <- function(f) {
    risk = NULL
    p <- pframe[model == f, risk]
    jackF <- pframe[model == f][[Rvar]]
    switch(method, quantile = {
      if (length(q) == 1) groups <- quantile(p, seq(0, 
                                                    1, 1/q)) else {
                                                      groups <- q
                                                    }
      # ADJUSTMENT MADE HERE: We set the groups to be unique
      groups <- unique(groups)
      xgroups <- (groups[-(length(groups))] + groups[-1])/2
      # this is where the error occurs
      pcut <- cut(p, groups, include.lowest = TRUE)
      if (x$response.type == "binary") {
        plotFrame = data.frame(Pred = tapply(p, pcut, 
                                             mean), Obs = pmin(1, pmax(0, tapply(jackF, 
                                                                                 pcut, mean))))
      } else {
        if (x$response.type == "survival") {
          censcode <- pframe[status == 0, status[1]]
          qfit <- prodlim::prodlim(prodlim::Hist(time, 
                                                 status, cens.code = censcode) ~ pcut, data = pframe)
          plotFrame = data.frame(Pred = tapply(p, pcut, 
                                               mean), Obs = predict(qfit, newdata = data.frame(pcut = levels(pcut)), 
                                                                    cause = 1, mode = "matrix", times = tp, 
                                                                    type = "cuminc"))
        } else {
          censcode <- pframe[status == 0, event[1]]
          qfit <- prodlim::prodlim(prodlim::Hist(time, 
                                                 event, cens.code = censcode) ~ pcut, data = pframe)
          n.cause <- match(x$cause, x$states)
          plotFrame = data.frame(Pred = tapply(p, pcut, 
                                               mean), Obs = predict(qfit, newdata = data.frame(pcut = levels(pcut)), 
                                                                    cause = n.cause, mode = "matrix", times = tp, 
                                                                    type = "cuminc"))
        }
      }
      attr(plotFrame, "quantiles") <- groups
      plotFrame
    }, nne = {
      if (round == TRUE) {
        if (!is.null(bandwidth) && bandwidth[1] >= 1) {
        } else {
          p <- round(p, 2)
        }
      }
      p <- na.omit(p)
      if (no <- length(attr(p, "na.action"))) warning("calPlot: removed ", 
                                                      no, " missing values in risk prediction.", 
                                                      call. = FALSE, immediate. = TRUE)
      if (is.null(bandwidth)) {
        if (length(unique(p)) == 1) bw <- 1 else bw <- prodlim::neighborhood(p)$bandwidth
      } else {
        bw <- bandwidth
      }
      if (bw >= 1) {
        plotFrame <- data.frame(Pred = mean(p), Obs = mean(jackF))
      } else {
        if (x$response.type == "binary") {
          nbh <- prodlim::meanNeighbors(x = p, y = jackF, 
                                        bandwidth = bw)
          plotFrame <- data.frame(Pred = nbh$uniqueX, 
                                  Obs = nbh$averageY)
        } else {
          if (cens.method == "local") {
            if (x$response.type == "survival") {
              censcode <- pframe[status == 0, status[1]]
              pfit <- prodlim::prodlim(prodlim::Hist(time, 
                                                     status, cens.code = censcode) ~ p, data = pframe, 
                                       bandwidth = bandwidth)
              plotFrame = data.frame(Pred = sort(unique(p)), 
                                     Obs = predict(pfit, newdata = data.frame(p = sort(unique(p))), 
                                                   cause = 1, mode = "matrix", times = tp, 
                                                   type = "cuminc"))
            } else {
              censcode <- pframe[status == 0, event[1]]
              pfit <- prodlim::prodlim(prodlim::Hist(time, 
                                                     event, cens.code = censcode) ~ p, data = pframe, 
                                       bandwidth = bandwidth)
              n.cause <- match(x$cause, x$states)
              plotFrame = data.frame(Pred = sort(unique(p)), 
                                     Obs = predict(pfit, newdata = data.frame(p = sort(unique(p))), 
                                                   cause = n.cause, mode = "matrix", 
                                                   times = tp, type = "cuminc"))
            }
          } else {
            nbh <- prodlim::meanNeighbors(x = p, y = jackF, 
                                          bandwidth = bw)
            plotFrame <- data.frame(Pred = nbh$uniqueX, 
                                    Obs = nbh$averageY)
          }
        }
      }
      attr(plotFrame, "bandwidth") <- bw
      plotFrame
    })
  }
  plotFrames <- lapply(modelnames, function(f) {
    getXY(f)
  })
  names(plotFrames) <- modelnames
  if (bars) {
    if ((is.logical(names[[1]]) && names[[1]] == TRUE) || 
        names[[1]] %in% c("quantiles.labels", "quantiles")) {
      qq <- attr(plotFrames[[1]], "quantiles")
      if (names[1] == "quantiles.labels") {
        pp <- seq(0, 1, 1/q)
        names <- paste0("(", sprintf("%1.0f", 
                                     100 * pp[-length(pp)]), ",", sprintf("%1.0f", 
                                                                          100 * pp[-1]), ")\n", sprintf("%1.1f", 
                                                                                                        100 * qq[-length(qq)]), " - ", sprintf("%1.1f", 
                                                                                                                                               100 * qq[-1]))
      }
      else names <- paste0(sprintf("%1.1f", 100 * 
                                     qq[-length(qq)]), " - ", sprintf("%1.1f", 
                                                                      100 * qq[-1]))
    }
  }
  out <- list(plotFrames = plotFrames, times = tp, control = control, 
              legend = legend, bars = bars, diag = diag, add = add, 
              names = names, method = method, axes = axes, percent = percent, 
              hanging = hanging, show.frequencies = show.frequencies, 
              col = col, ylim = ylim, xlim = xlim, ylab = ylab, xlab = xlab, 
              lwd = lwd, lty = lty, pch = pch, lty = lty, type = type, 
              NF = NF)
  if (method == "nne") 
    out <- c(out, list(bandwidth = sapply(plotFrames, function(x) attr(x, 
                                                                       "bandwidth"))))
  if (plot) {
    if (out$add[1] == FALSE && !out$bars[1]) {
      do.call("plot", control$plot)
    }
    if (out$diag[1] && !out$bars[1]) {
      segments(x0 = 0, y0 = 0, x1 = 1, y1 = 1, col = "gray77", 
               lwd = 2, xpd = FALSE)
    }
    if (out$diag[1] && !out$bars[1]) {
      segments(x0 = 0, y0 = 0, x1 = 1, y1 = 1, col = "gray77", 
               lwd = 2, xpd = FALSE)
    }
    if (out$bars[1]) {
      stopifnot(NF == 1)
      pf <- na.omit(plotFrames[[1]])
      Pred <- pf$Pred
      Obs <- pf$Obs
      if (is.character(legend[[1]]) || legend[[1]] == FALSE) {
        control$barplot$legend.text <- NULL
      }
      else {
        if (is.null(control$barplot$legend.text)) {
          control$barplot$legend.text <- control$legend$legend
        }
        control$barplot$args.legend <- control$legend
      }
      if (is.null(control$barplot$space)) 
        control$barplot$space <- rep(c(1, 0), length(Pred))
      PredObs <- c(rbind(Pred, Obs))
      control$barplot$height <- PredObs
      if (out$hanging) {
        control$barplot$offset <- c(rbind(0, Pred - Obs))
        minval <- min(Pred - Obs)
        if (minval < 0) 
          negY.offset <- 0.05 + seq(0, 1, 0.05)[prodlim::sindex(jump.times = seq(0, 
                                                                                 1, 0.05), eval.times = abs(minval))]
        else negY.offset <- 0
        control$barplot$ylim[1] <- min(control$barplot$ylim[1], 
                                       -negY.offset)
        control$names$y <- control$names$y - negY.offset
      }
      coord <- do.call("barplot", control$barplot)
      if (length(out$names) > 0 && (out$names[[1]] != FALSE) && 
          is.character(out$names)) {
        if (out$names[[1]] != FALSE && length(out$names) == 
            (length(coord)/2)) {
          mids <- rowMeans(matrix(coord, ncol = 2, byrow = TRUE))
          text(x = mids, y = control$names$y, out$names, 
               xpd = NA, cex = control$names$cex)
        }
      }
      if (out$hanging) {
        do.call("abline", control$abline)
      }
      if (out$show.frequencies) {
        if (out$hanging) {
          text(x = coord, cex = control$frequencies$cex, 
               pos = 3, y = (as.vector(rbind(Pred, Pred)) + 
                               rep(control$frequencies$offset, times = length(as.vector(coord))/2)), 
               paste(round(100 * c(rbind(Pred, Obs)), 0), 
                     ifelse(control$frequencies$percent, "%", 
                            ""), sep = ""), xpd = NA)
        }
        else {
          text(coord, pos = 3, c(rbind(Pred, Obs)) + 
                 control$frequencies$offset, cex = control$frequencies$cex, 
               paste(round(100 * c(rbind(Pred, Obs)), 0), 
                     ifelse(control$frequencies$percent, "%", 
                            ""), sep = ""), xpd = NA)
        }
      }
      coords <- list(xcoord = coord[, 1], ycoord = PredObs, 
                     offset = control$barplot$offset)
      out <- c(out, coords)
    }
    else {
      nix <- lapply(1:NF, function(f) {
        pf <- out$plotFrames[[f]]
        pf <- na.omit(pf)
        if (showRug) {
          do.call(graphics::rug, c(list(x = pf$Pred, 
                                        col = control$rug$col[f], lwd = control$rug$lwd[f] * 
                                          0.5)))
        }
        if (showPseudo) {
          if (!is.null(control$pseudo$density) & control$pseudo$density > 
              0) {
            control$pseudo$col <- prodlim::dimColor(control$pseudo$col[f], 
                                                    control$pseudo$density)
          }
          if ((gotcha <- match("density", names(control$pseudo), 
                               nomatch = 0)) > 0) {
            control$pseudo <- control$pseudo[-gotcha]
          }
          do.call(points, c(list(x = pframe[model == 
                                              modelnames[f]][, risk], y = pframe[model == 
                                                                                   modelnames[f]][[Rvar]]), control$pseudo))
        }
        lineF <- lapply(control$lines, function(x) x[[min(f, 
                                                          length(x))]])
        lineF$x <- pf$Pred
        lineF$y <- pf$Obs
        do.call("lines", lineF)
      })
      if (is.character(out$legend[[1]]) || out$legend[[1]] == 
          TRUE) {
        legend.coords <- do.call("legend", control$legend)
        if (!is.null(addtable2plot.DefaultArgs)) {
          if (is.null(control$addtable2plot[["x"]])) 
            control$addtable2plot[["x"]] <- legend.coords$rect$left + 
              legend.coords$rect$w
          if (is.null(control$addtable2plot[["y"]])) 
            control$addtable2plot[["y"]] <- legend.coords$rect$top - 
              legend.coords$rect$h
          do.call(plotrix::addtable2plot, control$addtable2plot)
        }
      }
    }
    if (out$axes) {
      if (out$percent) {
        if (is.null(control$axis1$labels)) {
          control$axis1$labels <- paste(100 * control$axis1$at, 
                                        "%")
        }
        if (is.null(control$axis2$labels)) {
          control$axis2$labels <- paste(100 * control$axis2$at, 
                                        "%")
        }
      }
      if (!out$bars) 
        do.call("axis", control$axis1)
      do.call("axis", control$axis2)
    }
  }
  class(out) <- "calibrationPlot"
  invisible(out)
}