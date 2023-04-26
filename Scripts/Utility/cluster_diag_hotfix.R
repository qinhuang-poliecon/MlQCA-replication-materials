###------------------------------------fix 1 ------------------------------------##
### this script slightly modifies 2 internal functions in the `setMethods` package
### essentially removed "nan" results that prob come from NA or uninformative variable selections
fix_clust_diag <- function (x, y, unit, cluster, necessity = FALSE, wicons = FALSE) 
{
  X <- xtabs(x ~ unit + cluster)
  Y <- xtabs(y ~ unit + cluster)
  p <- 2
  con.ragin <- function(X, Y) {
    sum(apply(cbind(X, Y), 1, min))/sum(X)
  }
  con.pool <- function(x, y) {
    z <- cbind(as.vector(x), as.vector(y))
    sum(apply(z, 1, min))/sum(x)
  }
  con.betw <- function(X, Y) {
    unlist(lapply(1:ncol(X), function(j) con.ragin(X[, j], 
                                                   Y[, j])))
  }
  con.with <- function(X, Y) {
    unlist(lapply(1:nrow(X), function(i) con.ragin(X[i, 
                                                     ], Y[i, ])))
  }
  cvr.ragin <- function(X, Y) {
    sum(apply(cbind(X, Y), 1, min))/sum(Y)
  }
  cvr.pool <- function(x, y) {
    z <- cbind(as.vector(x), as.vector(y))
    sum(apply(z, 1, min))/sum(y)
  }
  cvr.betw <- function(X, Y) {
    J <- ncol(X)
    unlist(lapply(1:J, function(j) cvr.ragin(X[, j], Y[, 
                                                       j])))
  }
  cvr.with <- function(X, Y) {
    N <- nrow(X)
    unlist(lapply(1:N, function(i) cvr.ragin(X[i, ], Y[i, 
                                                       ])))
  }
  if (!necessity) {
    dBP <- function(X, Y) {
      J <- ncol(X)
      bc <- con.betw(X, Y)
      bc <- bc[!is.nan(bc)]
      sqrt(sum(((bc/sum(bc)) - (1/J))^p))
    }
    dWP <- function(X, Y) {
      N <- nrow(X)
      wc <- con.with(X, Y)
      wc <- wc[!is.nan(wc)]
      sqrt(sum(((wc/sum(wc)) - (1/N))^p))
    }
  }
  else {
    dBP <- function(X, Y) {
      J <- ncol(X)
      bc <- cvr.betw(X, Y)
      sqrt(sum(((bc/sum(bc)) - (1/J))^p))
    }
    dWP <- function(X, Y) {
      N <- nrow(X)
      wc <- cvr.with(X, Y)
      sqrt(sum(((wc/sum(wc)) - (1/N))^p))
    }
  }
  unit <- as.character(unit)
  cluster <- as.character(cluster)
  CNRC <- data.frame(table(cluster))
  cnrc <- paste(as.character(CNRC[, 1]), " (", as.character(CNRC[, 
                                                                 2]), ") ", sep = "")
  CNRU <- data.frame(table(unit))
  cnru <- paste(as.character(CNRU[, 1]), " (", as.character(CNRU[, 
                                                                 2]), ") ", sep = "")
  if (!necessity) {
    r1 <- con.pool(x, y)
    r2 <- con.betw(X, Y)
    r3 <- dBP(X, Y)
    r4 <- con.with(X, Y)
    r5 <- dWP(X, Y)
    r6 <- list(pooled = cvr.pool(x, y), between = cvr.betw(X, 
                                                           Y), within = cvr.with(X, Y))
    r7 <- cnrc
    r8 <- cnru
  }
  else {
    r1 <- cvr.pool(x, y)
    r2 <- cvr.betw(X, Y)
    r3 <- dBP(X, Y)
    r4 <- cvr.with(X, Y)
    r5 <- dWP(X, Y)
    r6 <- list(pooled = con.pool(x, y), between = con.betw(X, 
                                                           Y), within = con.with(X, Y))
    r7 <- cnrc
    r8 <- cnru
  }
  r <- list(POCOS = r1, BECOS = r2, dBP = r3, WICONS = r4, 
            dWP = r5, Coverages = r6, wiconsprint = wicons, cluster_ids = r7, 
            unit_ids = r8)
  class(r) <- "clusterdiagnostics"
  return(r)
}

## fixInNamespace("cluster.diagnostics", "SetMethods")
old_clust_diag <- get("cluster.diagnostics", envir = asNamespace("SetMethods"))
environment(fix_clust_diag) <- environment(old_clust_diag)
attributes(fix_clust_diag) <- attributes(old_clust_diag)  # don't know if this is really needed
assignInNamespace("cluster.diagnostics", fix_clust_diag, ns="SetMethods")

###------------------------------------fix 2 ------------------------------------##
fix_pimplot <- function (data = NULL, results, outcome, incl.tt = NULL, ttrows = c(), 
          necessity = FALSE, sol = 1, all_labels = FALSE, markers = TRUE, 
          labcol = "black", jitter = FALSE, font = "sans", fontface = "italic", 
          fontsize = 3, crisp = FALSE, consH = FALSE, ...) 
{
  dots <- list(...)
  if (length(dots) != 0) {
    if ("neg.out" %in% names(dots)) {
      print("Argument neg.out is deprecated. The negated outcome is identified automatically from the minimize solution.")
    }
    if ("use.tilde" %in% names(dots)) {
      print("Argument use.tilde is deprecated. The usage of the tilde is identified automatically from the minimize solution.")
    }
  }
  if (length(grep("~", outcome)) > 0) {
    outcome <- outcome[grep("~", outcome)]
    outcome <- gsub("\\~", "", outcome)
    outcome <- unlist(outcome)
  }
  outcome <- toupper(outcome)
  if (!necessity) {
    data <- results$tt$initial.data
    if (is.null(incl.tt)) {
      if (length(ttrows) > 0) {
        oldtt <- results$tt$tt
        newtt <- oldtt[ttrows, ]
        P <- as.data.frame(results$tt$minmat)
        P <- P[colnames(P) %in% rownames(newtt)]
        if (results$options$neg.out | length(grep("~", 
                                                  results$call$outcome)) > 0) {
          neg.out = TRUE
          P$out <- results$tt$recoded.data[, outcome]
        }
        else {
          neg.out = FALSE
          P$out <- results$tt$recoded.data[, outcome]
        }
        n_c <- ncol(P) - 1
        par(ask = F)
        aux.plot <- function(i) {
          if (all_labels) {
            fil <- rownames(P)
          }
          else {
            fil <- rownames(P)
            fil[with(P, !(P[i] > 0.5))] <- ""
          }
          if (!neg.out) {
            xx <- xy.plot(P[, i, drop = FALSE], "out", data = P, 
                    xlab = paste("Row ", colnames(P)[i]), 
                    ylab = outcome, main = "Sufficiency Plot", 
                    labcol = labcol, jitter = jitter, consH = consH, 
                    font = font, fontface = fontface, fontsize = fontsize, 
                    labs = fil, crisp = crisp, shape = if (markers == 
                                                           FALSE) {
                      19
                    }
                    else {
                      ifelse((P$out < 0.5 & P[, i, drop = FALSE] > 
                                0.5), 9, 19)
                    })
          }
          else {
            xx <- xy.plot(P[, i, drop = FALSE], "out", data = P, 
                    xlab = paste("Row ", colnames(P)[i]), 
                    ylab = paste("~", outcome), main = "Sufficiency Plot", 
                    labcol = labcol, jitter = jitter, consH = consH, 
                    font = font, fontface = fontface, fontsize = fontsize, 
                    labs = fil, crisp = crisp, shape = if (markers == 
                                                           FALSE) {
                      19
                    }
                    else {
                      ifelse((P$out < 0.5 & P[, i, drop = FALSE] > 
                                0.5), 9, 19)
                    })
          }
          return(xx)
        }
        plotList <- list()
        for (i in 1:n_c) {
          plotList[[i]] <- aux.plot(i)
        }
      }else {
        P <- pimdata(results = results, outcome = outcome, 
                     sol = sol)
        n_c <- ncol(P) - 1
        par(ask = F)
        if (results$options$neg.out | length(grep("~", 
                                                  results$call$outcome)) > 0) {
          neg.out = TRUE
        }
        else {
          neg.out = FALSE
        }
        aux.plot <- function(i) {
          if (all_labels) {
            fil <- rownames(P)
          }
          else {
            fil <- rownames(P)
            fil[with(P, !(P[i] > 0.5))] <- ""
            if (i == n_c) {
              fil <- rownames(P)
              fil[with(P, !(P[i] < 0.5))] <- ""
            }
          }
          if (!neg.out) {
            xx <- xy.plot(P[, i, drop = FALSE], "out", data = P, 
                    xlab = colnames(P[i]), ylab = outcome, 
                    main = "Sufficiency Plot", labcol = labcol, 
                    jitter = jitter, consH = consH, font = font, 
                    fontface = fontface, fontsize = fontsize, 
                    labs = fil, crisp = crisp, shape = if (markers == 
                                                           FALSE) {
                      19
                    }
                    else {
                      ifelse((P$out < 0.5 & P[, i, drop = FALSE] > 
                                0.5), 9, 19)
                    })
          }
          else {
            xx <- xy.plot(P[, i, drop = FALSE], "out", data = P, 
                    xlab = colnames(P)[i], ylab = paste("~", 
                                                        outcome), main = "Sufficiency Plot", 
                    labcol = labcol, jitter = jitter, consH = consH, 
                    font = font, fontface = fontface, fontsize = fontsize, 
                    labs = fil, crisp = crisp, shape = if (markers == 
                                                           FALSE) {
                      19
                    }
                    else {
                      ifelse((P$out < 0.5 & P[, i, drop = FALSE] > 
                                0.5), 9, 19)
                    })
          }
          return(xx)
        }
        plotList <- list()
        for (i in 1:n_c) {
          plotList[[i]] <- aux.plot(i)
        }
      }
    }
    else {
      oldtt <- results$tt$tt
      suppressWarnings(oldtt$incl <- as.numeric(oldtt$incl))
      if (length(incl.tt) > 1) {
        paste("You introduced more than one inclusion cut for Truth Table rows. Please introduce only one!")
      }
      else {
        newtt <- oldtt[which(oldtt$incl > incl.tt), 
        ]
        P <- as.data.frame(results$tt$minmat)
        P <- P[colnames(P) %in% rownames(newtt)]
        if (results$options$neg.out | length(grep("~", 
                                                  results$call$outcome)) > 0) {
          neg.out = TRUE
          P$out <- results$tt$recoded.data[, outcome]
        }
        else {
          neg.out = FALSE
          P$out <- results$tt$recoded.data[, outcome]
        }
        n_c <- ncol(P) - 1
        par(ask = F)
        aux.plot <- function(i) {
          if (all_labels) {
            fil <- rownames(P)
          }
          else {
            fil <- rownames(P)
            fil[with(P, !(P[i] > 0.5))] <- ""
          }
          if (!neg.out) {
            xy.plot(P[, i, drop = FALSE], "out", data = P, 
                    xlab = paste("Row ", colnames(P)[i]), 
                    ylab = outcome, main = "Sufficiency Plot", 
                    labcol = labcol, jitter = jitter, consH = consH, 
                    font = font, fontface = fontface, fontsize = fontsize, 
                    labs = fil, crisp = crisp, shape = if (markers == 
                                                           FALSE) {
                      19
                    }
                    else {
                      ifelse((P$out < 0.5 & P[, i, drop = FALSE] > 
                                0.5), 9, 19)
                    })
          }
          else {
            xy.plot(P[, i, drop = FALSE], "out", data = P, 
                    xlab = paste("Row ", colnames(P)[i]), 
                    ylab = paste("~", outcome), main = "Sufficiency Plot", 
                    labcol = labcol, jitter = jitter, consH = consH, 
                    font = font, fontface = fontface, fontsize = fontsize, 
                    labs = fil, crisp = crisp, shape = if (markers == 
                                                           FALSE) {
                      19
                    }
                    else {
                      ifelse((P$out < 0.5 & P[, i, drop = FALSE] > 
                                0.5), 9, 19)
                    })
          }
        }
        for (i in 1:n_c) {
          print(aux.plot(i))
        }
      }
    }
  }
  else {
    if (is.null(data)) 
      stop("For analyses of necessity you need to provide the name of the dataframe!")
    P <- results$coms
    if (results$options$neg.out) {
      neg.out = TRUE
      P$out <- 1 - data[, outcome]
    }
    else {
      neg.out = FALSE
      P$out <- data[, outcome]
    }
    n_c <- ncol(P) - 1
    par(ask = F)
    aux.plot <- function(i) {
      if (all_labels) {
        fil <- rownames(P)
      }
      else {
        fil <- rownames(P)
        fil[with(P, !(P[, "out"] > 0.5))] <- ""
      }
      if (!neg.out) {
        xx <- xy.plot(P[, i, drop = FALSE], "out", data = P, 
                xlab = colnames(P)[i], ylab = outcome, necessity = TRUE, 
                main = "Necessity Plot", labcol = labcol, 
                jitter = jitter, font = font, consH = consH, 
                fontface = fontface, fontsize = fontsize, 
                labs = fil, crisp = crisp, shape = if (markers == 
                                                       FALSE) {
                  19
                }
                else {
                  ifelse((P$out > 0.5 & P[, i, drop = FALSE] < 
                            0.5), 9, 19)
                })
      }
      else {
        xx <- xy.plot(P[, i, drop = FALSE], "out", data = P, 
                xlab = colnames(P)[i], ylab = paste("~", outcome), 
                necessity = TRUE, main = "Necessity Plot", 
                labcol = labcol, jitter = jitter, font = font, 
                consH = consH, fontface = fontface, fontsize = fontsize, 
                labs = fil, crisp = crisp, shape = if (markers == 
                                                       FALSE) {
                  19
                }
                else {
                  ifelse((P$out > 0.5 & P[, i, drop = FALSE] < 
                            0.5), 9, 19)
                })
      }
      return(xx)
    }
    plotList <- list()
    for (i in 1:n_c) {
      plotList[[i]] <- aux.plot(i)
    }
  }
  return(plotList)
}


old_pimplot <- get("pimplot", envir = asNamespace("SetMethods"))
environment(fix_pimplot) <- environment(old_pimplot)
attributes(fix_pimplot) <- attributes(old_pimplot)  # don't know if this is really needed
assignInNamespace("pimplot", fix_pimplot, ns="SetMethods")


###------------------------------------fix 3 ------------------------------------##
fix_QCAfit <- function (x, y, cond.lab = NULL, necessity = TRUE, neg.out = FALSE, 
          product = FALSE, sol = 1, ttrows = c(), consH = FALSE) 
{
  if (is(x, "QCA_min")) {
    if (necessity == TRUE) 
      stop("You cannot calculate parameters of fit for necessity for a qca sufficient solution")
    if (!is.character(y)) 
      stop("When using qca object, the outcome must be of type character. \n                                 Please specify the outcome using its exact name in the dataframe. ")
    if (length(ttrows) > 0) {
      results <- x
      dt <- results$tt$initial.data
      oldtt <- results$tt$tt
      newtt <- oldtt[ttrows, ]
      P <- as.data.frame(results$tt$minmat)
      P <- P[colnames(P) %in% rownames(newtt)]
      if (results$options$neg.out | length(grep("~", results$call$outcome)) > 
          0) {
        P$out <- 1 - dt[, y]
      }
      else {
        P$out <- dt[, y]
      }
      nc <- ncol(P) - 1
      a <- data.frame(matrix(ncol = 4, nrow = nc))
      row.names(a) <- colnames(P[, 1:nc])
      colnames(a) <- c("Cons.Suf", "Cov.Suf", "PRI", "Cons.Suf(H)")
      for (i in 1:nc) {
        a[i, ] <- QCAfit(P[, i], P$out, cond.lab = cond.lab, 
                         necessity = necessity, neg.out = neg.out, 
                         product = product, sol = sol, consH = TRUE)
      }
    }
    else {
      X <- pimdata(results = x, outcome = y, sol = sol)
      nc <- ncol(X) - 1
      a <- data.frame(matrix(ncol = 4, nrow = nc))
      row.names(a) <- colnames(X[, 1:nc])
      colnames(a) <- c("Cons.Suf", "Cov.Suf", "PRI", "Cons.Suf(H)")
      for (i in 1:nc) {
        a[i, ] <- QCAfit(X[, i], X$out, cond.lab = cond.lab, 
                         necessity = necessity, neg.out = neg.out, 
                         product = product, sol = sol, consH = TRUE)
      }
    }
    if (consH == FALSE) {
      a <- a[, 1:length(a) - 1]
    }
    return(a)
  }
  else {
    x <- as.matrix(x)
    if (ncol(x) > 1) {
      nx <- 1 - x
      colnames(nx) <- paste("~", colnames(nx), sep = "")
      x <- cbind(x, nx)
    }
    v <- matrix(NA, length(x[, 1]), length(x[1, ]))
    out <- matrix(NA, length(x[1, ]), 8)
    if (neg.out == TRUE) {
      y <- 1 - y
    }
    for (i in 1:length(x[1, ])) {
      v[, i] <- pmin(x[, i], y, na.rm = T)
      out[i, 1] <- sum(v[, i], na.rm = T)/sum(x[, i], na.rm = T)
      out[i, 2] <- sum(v[, i], na.rm = T)/sum(y, na.rm = T)
      out[i, 3] <- sum(v[, i], na.rm = T)/sum(y, na.rm = T)
      out[i, 4] <- sum(v[, i], na.rm = T)/sum(x[, i], na.rm = T)
      out[i, 5] <- sum(1 - x[, i], na.rm = T)/sum(1 - v[, i], na.rm = T)
      p1 <- sum(pmin(x[, i], y, na.rm = T), na.rm = T) - 
        sum(pmin(x[, i], y, 1 - y, na.rm = T), na.rm = T)
      p2 <- sum(x[, i], na.rm = T) - sum(pmin(x[, i], y, 1 - y, na.rm = T), na.rm = T)
      p3 <- sum(v[, i] + sqrt(pmax((x[, i] - y), 0, na.rm = T) * 
                                x[, i]), na.rm = T)
      out[i, 6] <- p1/p2
      out[i, 7] <- sum(v[, i], na.rm = T)/p3
      out[i, 8] <- out[i, 1] * out[i, 6]
    }
    if (product == TRUE) {
      suf <- matrix(out[, c(1:2, 6:8)], nrow = ncol(x))
      colnames(suf) <- c("Cons.Suf", "Cov.Suf", "PRI", 
                         "Cons.Suf(H)", "PRODUCT")
      suf <- format(suf, digits = 3)
      storage.mode(suf) <- "numeric"
    }
    else {
      suf <- matrix(out[, c(1:2, 6:7)], nrow = ncol(x))
      colnames(suf) <- c("Cons.Suf", "Cov.Suf", "PRI", 
                         "Cons.Suf(H)")
      suf <- format(suf, digits = 3)
      storage.mode(suf) <- "numeric"
    }
    nec <- matrix(out[, 3:5], nrow = ncol(x))
    colnames(nec) <- c("Cons.Nec", "Cov.Nec", "RoN")
    nec <- format(nec, digits = 3)
    storage.mode(nec) <- "numeric"
    if (ncol(x) > 1) {
      rownames(suf) <- colnames(x)
      rownames(nec) <- colnames(x)
    }
    else {
      rownames(suf) <- cond.lab
      rownames(nec) <- cond.lab
    }
    if (consH == FALSE) {
      suf <- suf[, 1:ncol(suf) - 1]
    }
    if (necessity == FALSE) {
      return(suf)
    }
    else {
      return(nec)
    }
  }
}

old_QCAfit <- get("QCAfit", envir = asNamespace("SetMethods"))
environment(fix_QCAfit) <- environment(old_QCAfit)
attributes(fix_QCAfit) <- attributes(old_QCAfit)  # don't know if this is really needed
assignInNamespace("QCAfit", fix_QCAfit, ns="SetMethods")
