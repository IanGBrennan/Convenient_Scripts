# following:
# http://www.fromthebottomoftheheap.net/2011/06/12/additive-modelling-and-the-hadcrut3v-global-mean-temperature-series/

## read in the empirical data next
emp.vj <- as.data.frame(read.csv("BGB.Agamidae.UPDATED.Empirical.vj.with.CIs.txt", header=T, sep="\t"))

emp.vj
# Set the ages for the windows you ran, and trim the data set
simulation.dates <- seq(0, 20, 0.1) #here I've limited it to 20 million years, because that's our focal depth
all.emp.sim <- cbind.data.frame(emp.vj$lowerCI, emp.vj$upperCI, emp.vj$empirical.ratio.vj)

# Adjust the data table, and snip it down to ~20 MYA
emp <- NULL
emp <- all.emp.sim[0:length(simulation.dates), c("emp.vj$lowerCI","emp.vj$upperCI","emp.vj$empirical.ratio.vj")]
all.emp.sim.ci <- data.frame(x=(simulation.dates), y=emp)
plot(y.emp.vj.empirical.ratio.vj ~ x, data = all.emp.sim.ci, type = "o")
plot(y.emp.vj.empirical.ratio.vj ~ x, data = all.emp.sim.ci, type = "p")

# Run a couple GAMMs to get the flavor
m1 <- gamm(y.emp.vj.empirical.ratio.vj ~ s(x, k = 10), data = all.emp.sim.ci)
m2 <- gamm(y.emp.vj.empirical.ratio.vj ~ s(x, k = 6), data = all.emp.sim.ci)
m3 <- gamm(y.emp.vj.empirical.ratio.vj ~ s(x, k = 7), data = all.emp.sim.ci)

# Have a look at the summary if it interests you:
summary(m1$gam)

# Investigate the GAMM Diagnosis
with(all.emp.sim.ci, tsDiagGamm(m1, timevar=x, observed=y.emp.vj.empirical.ratio.vj))

# Plot the initial model to have a squizz
plot(m1$gam, residuals = T, pch = 19, cex = 0.75)

# Simulate dates as time slices, then predict values along the estimated GAM line
pdat <- with(all.emp.sim.ci, data.frame(x = seq(min(x), max(x),length = 200)))
p1 <- predict(m1$gam, newdata=pdat)
lines(p1~x, data=pdat, col = "red")


########################################################
## This is a loop for exploring the affect of changing
## the GAMM parameters, as a standard, use k=10
########################################################
total.models <- NULL
for (i in 20:35) {
  model <- gamm(y.emp.vj.empirical.ratio.vj ~ s(x, k=i), data = all.emp.sim.ci)
  total.models[[i]] <- model$lme
}
i=20
model.sel(total.models[[i]], 
          total.models[[i+1]], 
          total.models[[i+2]],
          total.models[[i+3]],
          total.models[[i+4]],
          total.models[[i+5]],
          total.models[[i+6]],
          total.models[[i+7]],
          total.models[[i+8]],
          total.models[[i+9]],
          total.models[[i+10]],
          total.models[[i+11]],
          total.models[[i+12]],
          total.models[[i+13]],
          total.models[[i+14]],
          total.models[[i+15]])

colnames(total.models) <- c("k", "logLik")
total.models <- as.data.frame(total.models)
best <- subset(total.models, logLik==max(total.models[,2]))
print(best)
best.k <- best$k
best.gamm <- gamm(y.emp.vj.empirical.ratio.vj ~ s(x, k = best.k), data = all.emp.sim.ci)
plot(best.gamm$gam, residuals = T, pch = 19, cex = 0.75)
#####################################################


#############################################
## Functions for derivatives of GAM models ##
#############################################
Deriv <- function(mod, n = 200, eps = 1e-7, newdata) {
  if(isTRUE(inherits(mod, "list")))
    mod <- mod$gam
  m.terms <- attr(terms(mod), "term.labels")
  if(missing(newdata)) {
    newD <- sapply(model.frame(mod)[, m.terms, drop = FALSE],
                   function(x) seq(min(x), max(x), length = n))
    names(newD) <- m.terms
  } else {
    newD <- newdata
  }
  X0 <- predict(mod, data.frame(newD), type = "lpmatrix")
  newD <- newD + eps
  X1 <- predict(mod, data.frame(newD), type = "lpmatrix")
  Xp <- (X1 - X0) / eps
  Xp.r <- NROW(Xp)
  Xp.c <- NCOL(Xp)
  ## dims of bs
  bs.dims <- sapply(mod$smooth, "[[", "bs.dim") - 1
  # number of smooth terms
  t.labs <- attr(mod$terms, "term.labels")
  nt <- length(t.labs)
  ## list to hold the derivatives
  lD <- vector(mode = "list", length = nt)
  names(lD) <- t.labs
  for(i in seq_len(nt)) {
    Xi <- Xp * 0
    want <- grep(t.labs[i], colnames(X1))
    Xi[, want] <- Xp[, want]
    df <- Xi %*% coef(mod)
    df.sd <- rowSums(Xi %*% mod$Vp * Xi)^.5
    lD[[i]] <- list(deriv = df, se.deriv = df.sd)
    ## Xi <- Xp * 0 ##matrix(0, nrow = Xp.r, ncol = Xp.c)
    ## J <- bs.dims[i]
    ## Xi[,(i-1) * J + 1:J + 1] <- Xp[,(i-1) * J + 1:J +1]
    ## df <- Xi %*% coef(mod)
    ## df.sd <- rowSums(Xi %*% mod$Vp * Xi)^.5
    ## lD[[i]] <- list(deriv = df, se.deriv = df.sd)
  }
  class(lD) <- "Deriv"
  lD$gamModel <- mod
  lD$eps <- eps
  lD$eval <- newD - eps
  return(lD)
}

confint.Deriv <- function(object, term, alpha = 0.05, ...) {
  l <- length(object) - 3
  term.labs <- names(object[seq_len(l)])
  if(missing(term))
    term <- term.labs
  Term <- match(term, term.labs)
  ##term <- term[match(term, term.labs)]
  if(any(miss <- is.na(Term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  ## if(is.na(term))
  ##     stop("'term' not a valid model term.")
  res <- vector(mode = "list", length = length(term))
  names(res) <- term
  residual.df <- length(object$gamModel$y) - sum(object$gamModel$edf)
  tVal <- qt(1 - (alpha/2), residual.df)
  ## tVal <- qt(1 - (alpha/2), object$gamModel$df.residual)
  for(i in seq_along(term)) {
    upr <- object[[term[i]]]$deriv + tVal * object[[term[i]]]$se.deriv
    lwr <- object[[term[i]]]$deriv - tVal * object[[term[i]]]$se.deriv
    res[[term[i]]] <- list(upper = drop(upr), lower = drop(lwr))
  }
  res$alpha = alpha
  res
}

signifD <- function(x, d, upper, lower, eval = 0) {
  miss <- upper > eval & lower < eval
  incr <- decr <- x
  want <- d > eval
  incr[!want | miss] <- NA
  want <- d < eval
  decr[!want | miss] <- NA
  list(incr = incr, decr = decr)
}

plot.Deriv <- function(x, alpha = 0.05, polygon = TRUE,
                       sizer = FALSE, term, eval = 0, lwd = 3,
                       col = "lightgrey", border = col,
                       ylab, xlab, ...) {
  l <- length(x) - 3
  ## get terms and check specified (if any) are in model
  term.labs <- names(x[seq_len(l)])
  if(missing(term))
    term <- term.labs
  Term <- match(term, term.labs)
  if(any(miss <- is.na(Term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  if(all(is.na(Term)))
    stop("All terms in 'term' not found in model.")
  l <- sum(!miss)
  nplt <- n2mfrow(l)
  ## tVal <- qt(1 - (alpha/2), x$gamModel$df.residual)
  residual.df <- length(x$gamModel$y) - sum(x$gamModel$edf)
  tVal <- qt(1 - (alpha/2), residual.df)
  if(missing(ylab))
    ylab <- expression(italic(hat(f)*"'"*(x)))
  if(missing(xlab)) {
    xlab <- attr(terms(x$gamModel), "term.labels")[Term]
    names(xlab) <- xlab
  }
  layout(matrix(seq_len(l), nrow = nplt[1], ncol = nplt[2]))
  CI <- confint(x, term = term, alpha = alpha)
  for(i in seq_along(term)) {
    ## for(i in seq_len(l)) {
    upr <- CI[[term[i]]]$upper
    lwr <- CI[[term[i]]]$lower
    ylim <- range(upr, lwr)
    plot(x$eval[,term[i]], x[[term[i]]]$deriv, type = "n",
         ylim = ylim, ylab = ylab, xlab = xlab[term[i]], ...)
    if(isTRUE(polygon)) {
      polygon(c(x$eval[,term[i]], rev(x$eval[,term[i]])),
              c(upr, rev(lwr)), col = col, border = border)
    } else {
      lines(x$eval[,term[i]], upr, lty = "dashed")
      lines(x$eval[,term[i]], lwr, lty = "dashed")
    }
    abline(h = 0, ...)
    if(isTRUE(sizer)) {
      lines(x$eval[,term[i]], x[[term[i]]]$deriv, lwd = 1)
      S <- signifD(x[[term[i]]]$deriv, x[[term[i]]]$deriv, upr, lwr,
                   eval = eval)
      lines(x$eval[,term[i]], S$incr, lwd = lwd, col = "blue")
      lines(x$eval[,term[i]], S$decr, lwd = lwd, col = "red")
    } else {
      lines(x$eval[,term[i]], x[[term[i]]]$deriv, lwd = 2)
    }
  }
  layout(1)
  invisible(x)
}
#############################################

# Now let's visualize the trends in the data
m1.d <- Deriv(m1, n=200)
plot(m1.d, sizer=T, alpha=0.05)

plot(y.emp.vj.empirical.ratio.vj ~ x, data = all.emp.sim.ci, type = "p")
lines(p1 ~ x, data = pdat)
CI <- confint(m1.d, alpha = 0.05)
S <- signifD(p1, m1.d$x$deriv, CI$x$upper, CI$x$lower,
             eval = 0)
lines(S$incr ~ x, data = pdat, lwd = 3, col = "blue")
lines(S$decr ~ x, data = pdat, lwd = 3, col = "red")
#










## look at autocorrelation in residuals:
acf(resid(m1$lme, type = "normalized"))
## ...wait... look at the plot, only then do...
pacf(resid(m1$lme, type = "normalized"))
## seems like like sharp cut-off in ACF and PACF - AR terms probably best

m2 <- gamm(Annual ~ s(Year, k = 30), data = gtemp,
           correlation = corARMA(form = ~ Year, p = 1))







## load custom functions
tmp <- tempfile()
download.file("https://github.com/gavinsimpson/random_code/raw/master/derivFun.R", tmp)
source(tmp)

tmp <- tempfile()
download.file("https://github.com/gavinsimpson/random_code/raw/master/tsDiagGamm.R", tmp)
source(tmp)



## Global temperatures
URL <- url("http://www.cru.uea.ac.uk/cru/data/temperature/HadCRUT3v-gl.dat")
tmp <- tempfile()
clim <- download.file("http://www.cru.uea.ac.uk/cru/data/temperature/HadCRUT3v-gl.dat", tmp)
load(clim)

gtemp <- read.table(clim, fill = TRUE)
## Don't need the even rows
gtemp <- gtemp[-seq(2, nrow(gtemp), by = 2), ]
## set the Year as rownames
rownames(gtemp) <- gtemp[,1]
## Add colnames
colnames(gtemp) <- c("Year", month.abb, "Annual")
## Data for 2011 incomplete so work only with 1850-2010 data series
gtemp <- gtemp[-nrow(gtemp), ]
## Plot the data
ylab <- expression(Temperature~Anomaly~(1961-1990)~degree*C)
plot(Annual ~ Year, data = gtemp, type = "o", ylab = ylab)

plot()
