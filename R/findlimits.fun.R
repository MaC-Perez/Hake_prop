"findlimits.fun" <-
  function (data, dist = 10, plot = NULL, what.out = "lines", correction = NULL) 
{
  if (is.na(match(what.out, c("lines", "window")))) 
    stop("what has to be either lines or window")
  if (!(is.null(plot) | sum(match(c("contour", "limits"), plot), 
                            na.rm = T) > 0)) 
    stop("plot has to be either NULL, contour or limits")
  what.plot <- ifelse(is.null(plot), 1, ifelse(plot == "contour", 
                                               2, 3))
  if (is.null(data$x)) {
    warning(paste("File", data, "does not contain x, applying deg2km"))
    data <- deg2km(data)
  }
  data <- data[!is.na(data$x), ]
  data <- data[!is.na(data$y), ]
  buffer.add <- c(-2 * dist, 2 * dist)
  range.x <- range(data$x) + buffer.add
  range.y <- range(data$y) + buffer.add
  data.ppp <- ppp(data$x, data$y, xrange = range.x, yrange = range.y)
  data.dist <- distmap(data.ppp)
  data.srf <- list(x = data.dist$xcol, y = data.dist$yrow, 
                   z = t(data.dist$v))
  data.lim <- contourLines(data.srf, nlevels = 1, levels = dist)
  if (length(data.lim) != 1) {
    warning("The limits do not form a single line, check that holes are well defined")
    outlist <- vector("list", length = length(data.lim))
    for (i in c(1:length(data.lim))) {
      d1 <- list(x = data.lim[[i]]$x[c(length(data.lim[[i]]$x):2)], 
                 y = data.lim[[i]]$y[c(length(data.lim[[i]]$x):2)])
      outlist[[i]] <- d1
    }
  }
  else outlist <- list(x = data.lim[[1]]$x[c(length(data.lim[[1]]$x):2)], 
                       y = data.lim[[1]]$y[c(length(data.lim[[1]]$x):2)])
  if (what.plot == 2) 
    contour(data.srf)
  if (what.out == "window") {
    outside.owin <- owin(poly = outlist)
    if (!is.null(correction)) 
      outside.owin <- intersect.owin(outside.owin, correction)
    if (what.plot == 2) {
      plot(outside.owin)
      points(data.ppp)
    }
    outside.owin
  }
  else {
    if (what.plot == 3) {
      plot(seq(range.x[1], range.x[2], length = 10), seq(range.y[1], 
                                         range.y[2], length = 10), type = "n", axes = F, 
           xlab = "Long", ylab = "Lat")
      points(data$x, data$y, pch = 19, cex = 0.7)
      box()
      if (is.na(match("x", names(outlist)))) {
        for (i in c(1:length(outlist))) {
          lines(outlist[[i]]$x, outlist[[i]]$y, col = i)
          text(outlist[[i]]$x[1], outlist[[i]]$y[1], labels=i, col=i, cex=1.2)
        }
      }
      else lines(outlist$x, outlist$y, col = 2)
    }
    invisible(outlist)
  }
}
