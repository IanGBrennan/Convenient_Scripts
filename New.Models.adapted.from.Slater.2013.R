require(mvtnorm)
#############################################################################################################################

timeshiftTree <- function(phy, breakPoint, endRate) {
		if(is.ultrametric(phy)) {	
    	times <- branching.times(phy);
    } 
    if(!is.ultrametric(phy)) {
    	times <- BranchingTimesFossil(phy)[1:phy$Nnode];
    }
        for (i in 1:length(phy$edge.length)) {
        bl <- phy$edge.length[i]
        age = times[which(names(times) == phy$edge[i, 1])]
        if ((age - bl) < breakPoint) 
            phy$edge.length[i] <- (age - min(age, breakPoint)) * 
                1 + (min(age, breakPoint) - (age - bl)) * endRate
    }
    phy
}

#############################################################################################################################

releaseTree <- function (phy, alpha, breakPoint, endRate = 1) {
	if(is.ultrametric(phy)) {	
    	times <- branching.times(phy);
    } 
        if(!is.ultrametric(phy)) {
    	times <- BranchingTimesFossil(phy)[1: phy$Nnode];
    }
        names(times) <- (as.numeric(names(times)))
    Tmax <- max(times)
    phy2 <- phy
    for (i in 1:length(phy$edge.length)) {
        bl <- phy$edge.length[i]
        age = times[which(names(times) == phy$edge[i, 1])]
        if (age > breakPoint) {
			t1 = max(times) - age
			t2 = t1 + min(bl, (age-breakPoint))
		        	ou.length <- (1/(2 * alpha)) * exp(-2 * alpha * 
            (Tmax - t2)) * (1 - exp(-2 * alpha * t2)) - (1/(2 * 
            alpha)) * exp(-2 * alpha * (Tmax - t1)) * (1 - exp(-2 * 
            alpha * t1)) 
                          bm.length <-  (bl - min(bl, (age-breakPoint))) * endRate;
			phy2$edge.length[i] = ou.length + bm.length;
		} else {
			phy2$edge.length[i] <- phy$edge.length[i] * endRate;
		}	
	}
	return(phy2);
}	


#############################################################################################################################

fitContinuous_paleo <- function (phy, data, data.names = NULL, model = c("BM", "OU", 
    "lambda", "kappa", "delta", "EB", "white", "trend","timeshift", "release", 
    "releaseradiate", "radiate.constrain", "SRC", "TRC", "altTRC", 
    "BM1.OU.BM2", "BM1.OU.BM1", "TRC.time", "BM1.OU.BM1.time", "Burst.Bind",
    "Burst.Bind.Time", "Burst.Burst", "SRC_ME", "TRC_ME", "BMOU", "BMOUi", "BMOU_ME", "BMOUi_ME"), shift.time = 30, shift.time2 = 30, bounds = NULL, 
    meserr = NULL) 
{
    model <- match.arg(model)
    td <- treedata(phy, data, data.names, sort = T)
    ntax = length(td$phy$tip.label)
    if (is.null(meserr)) {
        me = td$data
        me[] = 0
        meserr = me
    }
    else if (length(meserr) == 1) {
        me = td$data
        me[] = meserr
        meserr = me
    }
    else if (is.vector(meserr)) {
        if (!is.null(names(meserr))) {
            o <- match(rownames(td$data), names(meserr))
            if (length(o) != ntax) 
                stop("meserr is missing some taxa from the tree")
            meserr <- as.matrix(meserr[o, ])
        }
        else {
            if (length(meserr) != ntax) 
                stop("No taxon names in meserr, and the number of taxa does not match the tree")
            me <- td$data
            me[] = meserr
            meserr = me
        }
    }
    else {
        if (!is.null(rownames(meserr))) {
            o <- match(rownames(td$data), rownames(meserr))
            meserr = meserr[o, ]
        }
        else {
            if (sum(dim(meserr) != dim(td$data)) != 0) 
                stop("No taxon names in meserr, and the number of taxa does not match the tree")
            print("No names in meserr; assuming that taxa are in the same order as tree")
        }
    }
    ds <- list()
    ds$tree <- td$phy
    cat("Fitting ", model, "model:\n")
    bounds.default <- matrix(c(1e-10, 1000, # beta
                               1e-07, 1,   # lambda
                               1e-06, 1,   # kappa
                               1e-05, 2.999999, # delta
                               1e-10, 100, # alpha
                               -3, -0.001, # a
                               1e-10, 100, # nv
                               -100, 100, # mu
                               1e-10, 100,
                               #0.1, max(nodeHeights(ds$tree))
                               4, 15, # scalar
                               0.01, 0.1), # measurement error (MErr)
                             nrow = 11, ncol = 2, byrow = TRUE)
    rownames(bounds.default) <- c("beta", "lambda", "kappa", 
        "delta", "alpha", "a", "nv", "mu", "scalar", "shift.timing", "MErr")
    colnames(bounds.default) <- c("min", "max")
    if (is.null(bounds)) {
        bounds <- bounds.default
    }
    else {
        if (class(bounds) != "list") {
            stop("Please specify user defined parameter bounds as a list()")
        }
        else {
            specified <- !c(is.null(bounds$beta), is.null(bounds$lambda), 
                is.null(bounds$kappa), is.null(bounds$delta), 
                is.null(bounds$alpha), is.null(bounds$a), is.null(bounds$nv), 
                is.null(bounds$mu))
            bounds.user <- matrix(c(bounds$beta, bounds$lambda, 
                bounds$kappa, bounds$delta, bounds$alpha, bounds$a, 
                bounds$nv, bounds$mu, bounds$shift.timing), nrow = sum(specified), 
                ncol = 2, byrow = TRUE)
            rownames(bounds.user) <- c("beta", "lambda", "kappa", 
                "delta", "alpha", "a", "nv", "mu", "shift.timing", "MErr")[specified]
            colnames(bounds.user) <- c("min", "max")
            bounds <- bounds.default
            bounds[specified, ] <- bounds.user
        }
    }
    ds$bounds <- data.frame(t(bounds))
    ds$model <- model
    ds$shift.time <- shift.time
    ds$shift.time2<- shift.time2;
    result <- list()
    for (i in 1:ncol(td$data)) {
        ds$data = td$data[, i]
        ds$meserr = meserr[, i]
        result[[i]] <- fitContinuousModel_paleo(ds, print = print)
        if (!is.null(colnames(td$data))) 
            names(result)[i] <- colnames(td$data)[i]
        else names(result)[i] <- paste("Trait", i, sep = "")
    }
    result
}


#############################################################################################################################

fitContinuousModel_paleo <- function (ds, print = TRUE) 
{
    bounds <- ds$bounds
    model <- ds$model
    n <- length(ds$data)
    beta.start <- var(ds$data)/max(branching.times(ds$tree))
    out <- NULL
    y <- ds$data
    tree <- ds$tree
    meserr <- ds$meserr
    shift.time <- ds$shift.time
    shift.time2<- ds$shift.time2;
    n <- length(y)
    if (model == "BM") {
        k <- 2
        vcv <- vcv.phylo(tree)
        start = log(beta.start)
        lower = log(bounds[1, "beta"])
        upper = log(bounds[2, "beta"])
        foo <- function(x) {
            vv <- exp(x) * vcv
            diag(vv) <- diag(vv) + meserr^2
            mu <- phylogMean(vv, y)
            mu <- rep(mu, n)
            -dmvnorm(y, mu, vv, log = T)
        }
        o <- optim(foo, p = start, lower = lower, upper = upper, 
            method = "L")
   		root.state <- as.numeric(ml.root(tree = tree, model = model , y = y, meserr, params = o, shift.time));
        results <- list(lnl = -o$value, root.state = root.state, beta = exp(o$par))
    }
    else if (model == "lambda") {
        k <- 3
        start = log(c(beta.start, 0.5))
        lower = log(bounds[1, c("beta", "lambda")])
        upper = log(bounds[2, c("beta", "lambda")])
        foo <- function(x) {
            vcv <- vcv.phylo(tree)
            index <- matrix(TRUE, n, n)
            diag(index) <- FALSE
            vcv[index] <- vcv[index] * exp(x[2])
            vv <- exp(x[1]) * vcv
            diag(vv) <- diag(vv) + meserr^2
            mu <- phylogMean(vv, y)
            mu <- rep(mu, n)
            -dmvnorm(y, mu, vv, log = T)
        }
        o <- optim(foo, p = start, lower = lower, upper = upper, 
            method = "L")
        results <- list(lnl = -o$value, beta = exp(o$par[1]), 
            lambda = exp(o$par[2]))
    }
    else if (model == "kappa") {
        k <- 3
        start = log(c(beta.start, 0.5))
        lower = log(bounds[1, c("beta", "kappa")])
        upper = log(bounds[2, c("beta", "kappa")])
        foo <- function(x) {
            t <- transform(tree, "kappa", exp(x[2]))
            vcv <- vcv.phylo(t)
            vv <- exp(x[1]) * vcv
            diag(vv) <- diag(vv) + meserr^2
            mu <- phylogMean(vv, y)
            mu <- rep(mu, n)
            -dmvnorm(y, mu, vv, log = T)
        }
        o <- optim(foo, p = start, lower = lower, upper = upper, 
            method = "L")
        results <- list(lnl = -o$value, beta = exp(o$par[1]), 
            lambda = exp(o$par[2]))
    }
    else if (model == "delta") {
        k <- 3
        start = log(c(beta.start, 0.5))
        lower = log(bounds[1, c("beta", "delta")])
        upper = log(bounds[2, c("beta", "delta")])
        foo <- function(x) {
            t <- transform(tree, "delta", exp(x[2]))
            vcv <- vcv.phylo(t)
            vv <- exp(x[1]) * vcv
            diag(vv) <- diag(vv) + meserr^2
            mu <- phylogMean(vv, y)
            mu <- rep(mu, n)
            -dmvnorm(y, mu, vv, log = T)
        }
        o <- optim(foo, p = start, lower = lower, upper = upper, 
            method = "L")
        results <- list(lnl = -o$value, beta = exp(o$par[1]), 
            delta = exp(o$par[2]))
    }
    else if (model == "white") {
        k <- 2
        start = c(mean(y), log(var(y)))
        lower = c(-Inf, log(bounds[1, "nv"]))
        upper = c(Inf, log(bounds[2, "nv"]))
        lnl.noise <- function(p, x, se) {
            root <- p[1]
            vs <- exp(p[2])
            n <- length(x)
            VV <- diag(vs, nrow = n)
            diag(VV) <- diag(VV) + se^2
            M <- rep(root, n)
            -dmvnorm(x, M, VV, log = TRUE)
        }
        o <- optim(start, fn = lnl.noise, x = y, se = meserr, 
            lower = lower, upper = upper, method = "L")
        results <- list(lnl = -o$value, mean = o$par[1], nv = exp(o$par[2]))
    }
    else if (model == "trend") {
        k <- 3
        vcv <- vcv.phylo(tree)
        ww <- lm(y ~ diag(vcv))
        p0 <- c(phylogMean(vcv, y), var(y)/max(branching.times(tree)), 
            coef(ww)[2])
        if (is.na(p0[3])) {
            p0[3] <- 0
            if (is.ultrametric(tree)) 
                cat("WARNING: Cannot estimate a trend with an ultrametric tree; lnl will be the same as the BM model")
        }
        lower = c(-Inf, log(bounds[1, "beta"]), bounds[1, "mu"])
        upper = c(Inf, log(bounds[2, "beta"]), bounds[2, "mu"])
        lnl.BMtrend <- function(p, vcv, x, se) {
            root <- p[1]
            vs <- exp(p[2])
            ms <- p[3]
            VV <- vs * vcv
            diag(VV) <- diag(VV) + se^2
            n <- length(x)
            M <- root + ms * diag(vcv)
            -dmvnorm(x, M, VV, log = TRUE)
        }
        o <- optim(p0, fn = lnl.BMtrend, vcv = vcv, x = y, se = meserr, 
            lower = lower, upper = upper, method = "L")
        names(o$par) <- NULL
        results <- list(lnl = -o$value, mean = o$par[1], beta = exp(o$par[2]), 
            mu = o$par[3])
    }
    else if (model == "OU") {
        k <- 3
        start = log(c(beta.start, 0.5))
        lower = log(bounds[1, c("beta", "alpha")])
        upper = log(bounds[2, c("beta", "alpha")])
		vmat <- vcv(tree)
        foo <- function(x) {
            vv <- ouMatrix(vmat, exp(x[2])) * exp(x[1])
            diag(vv) <- diag(vv) + meserr^2
            mu <- phylogMean(vv, y)
            mu <- rep(mu, n)
            -dmvnorm(y, mu, vv, log = T)
        }
        outTries <- list()
        start = c(log(beta.start), -50)
        outTries[[1]] <- optim(foo, p = start, lower = lower, 
            upper = upper, method = "L")
        tv <- var(y)
        start = log(c(tv * 2000, 1000))
        outTries[[2]] <- optim(foo, p = start, lower = lower, 
            upper = upper, method = "L")
        for (i in 1:10) {
            while (1) {
                lower = c(runif(2, min = -20, max = -1))
                upper = lower + runif(2, min = 0, max = 10)
                start = c(runif(1, min = lower[1], max = upper[1]), 
                  runif(1, min = lower[2], max = upper[2]))
                te <- try(outTries[[i + 2]] <- optim(foo, p = start, 
                  lower = lower, upper = upper, method = "L"), 
                  silent = T)
                if (class(te) != "try-error") 
                  break
            }
        }
        atry <- -5:4
        stry <- log(tv * 2 * exp(atry))
        for (i in 1:10) {
            while (1) {
                lower = c(-20, -20)
                upper = c(10, 10)
                start = c(stry[i], atry[i])
                te <- try(outTries[[i + 12]] <- optim(foo, p = start, 
                  lower = lower, upper = upper, method = "L"), 
                  silent = T)
                if (class(te) != "try-error") 
                  break
            }
        }
        ntries <- 22
        ltry <- numeric(ntries)
        lsol <- matrix(nrow = ntries, ncol = 2)
        for (j in 1:ntries) {
            ltry[j] <- outTries[[j]]$value
            lsol[j, ] <- exp(outTries[[j]]$par)
        }
        ltd <- ltry - min(ltry)
        b <- min(which(ltry == min(ltry)))
        gc <- which(ltd < 0.01)
        us <- lsol[gc, 1]
        usc <- sum((us - min(us)) > 0.01)
        out <- outTries[[b[1]]]
        if (usc > 1) {
            out$message = "Warning: likelihood surface is flat."
        }
        if (out$convergence != 0) {
            out$message = "Warning: may not have converged to a proper solution."
        }
         root.state <- as.numeric(ml.root(tree = tree, model = model , y = y, meserr, params = out, shift.time));

        results <- list(lnl = -out$value, root.state = root.state, beta = exp(out$par[1]), 
            alpha = exp(out$par[2]), convergence = out$convergence, 
            message = out$message, k = k)
    }
    else if (model == "EB") {
        k <- 3
        start = c(log(beta.start), 0.01)
        lower = c(log(bounds[1, "beta"]), bounds[1, "a"])
        upper = c(log(bounds[2, "beta"]), bounds[2, "a"])
        foo <- function(x) {
            t <- rescale(tree, "EB", a = x[2])
            vcv <- vcv.phylo(t)
            vv <- exp(x[1]) * vcv
            diag(vv) <- diag(vv) + meserr^2
            mu <- phylogMean(vv, y)
            mu <- rep(mu, n)
            -dmvnorm(y, mu, vv, log = T)
        }
        o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
        root.state <- as.numeric(ml.root(tree = tree, model = model , y = y, meserr, params = o, shift.time));
        results <- list(lnl = -o$value,root.state = root.state, beta = exp(o$par[1]), 
            a = o$par[2])
    }
    else if (model == "Burst.Burst") {
      k <- 4
      start = c(log(c(beta.start, 0.1)), 1, 0.01)
      lower = c(log(bounds[1, c("beta", "beta")]), bounds[1, c("a", "a")])
      upper = c(log(bounds[2, c("beta", "beta")]), bounds[2, c("a", "a")])
      r.m <- split.vcv(tree, shift.time);
      foo <- function(x) {
        t.init1 <- rescale(tree, "EB", a = x[3])
        t1 <- rescale(t.init1, "depth", depth=max(nodeHeights(tree)))
        t.init2 <- rescale(tree, "EB", a = x[4])
        t2 <- rescale(t.init2, "depth", depth=max(nodeHeights(tree)))
        release.mat1 <- split.vcv(t1, shift.time)
        release.mat2 <- split.vcv(t2, shift.time)
        eb1.mat <- release.mat1[[1]] * exp(x[1])
        eb2.mat <- release.mat2[[2]] * exp(x[2])
        vv <- eb1.mat + eb2.mat
        mu <- phylogMean(vv, y)
        mu <- rep(mu, n)
        -dmvnorm(y, mu, vv, log = T)
      }
      o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
      #ml.ouMatrix <- r.m[[2]] * exp(o$par[1])
      ml.vcv <- ((split.vcv(rescale(rescale(tree, "EB", a=o$par[3]), "depth", depth=max(nodeHeights(tree))))[[1]]) * o$par[1]) 
      + ((split.vcv(rescale(rescale(tree, "EB", a=o$par[4]), "depth", depth=max(nodeHeights(tree))))[[2]]) * o$par[2])
      #ml.vcv <- exp(o$par[1]) * ((ouMatrix(ml.ouMatrix, exp(o$par[2]))) + (r.m[[1]]))
      root.state <- phylogMean(ml.vcv, y)
      #root.state <- as.numeric(ml.root(tree = tree, model = model , y = y, meserr, params = o, shift.time));
      results <- list(lnl = -o$value, root.state = root.state, beta1 = exp(o$par[1]),  
                      beta2 = exp(o$par[2]), a1 = o$par[3], a2 = o$par[4]);
    }
    else if (model == "Burst.Bind") {
        k <- 4
        start = c(log(c(beta.start, 0.1)), 1, 0.01)
        lower = c(log(bounds[1, c("beta", "alpha")]), bounds[1, "a"])
        upper = c(log(bounds[2, c("beta", "alpha")]), bounds[2, "a"])
        r.m <- split.vcv(tree, shift.time);
        foo <- function(x) {
          t.init <- rescale(tree, "EB", a = x[3])
          t <- rescale(t.init, "depth", depth=max(nodeHeights(tree)))
          release.mat1 <- split.vcv(t.init, shift.time)
          #release.mat1 <- split.vcv(t.init, (shift.time/max(nodeHeights(tree)))*max(nodeHeights(t.init)))
          release.mat2 <- split.vcv(tree, shift.time)
          eb.mat <- release.mat1[[1]] * exp(x[1])
          bm2ou <- release.mat2[[2]] * exp(x[1])
          ou.mat <- ouMatrix(bm2ou, exp(x[2]))
          vv <- eb.mat + ou.mat
          mu <- phylogMean(vv, y)
          mu <- rep(mu, n)
          -dmvnorm(y, mu, vv, log = T)
        }
        o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
        ml.ouMatrix <- r.m[[2]] * exp(o$par[1])
        ml.vcv <- ((exp(o$par[1])) * (split.vcv((rescale(tree, "EB", a=o$par[3])), shift.time)[[1]])) + ((ouMatrix(ml.ouMatrix, exp(o$par[2]))))
        #ml.vcv <- exp(o$par[1]) * ((ouMatrix(ml.ouMatrix, exp(o$par[2]))) + (r.m[[1]]))
        root.state <- phylogMean(ml.vcv, y)
        #root.state <- as.numeric(ml.root(tree = tree, model = model , y = y, meserr, params = o, shift.time));
        results <- list(lnl = -o$value, root.state = root.state, beta = exp(o$par[1]),  
                      alpha = exp(o$par[2]), a = o$par[3]);
    }
    
    else if (model == "Burst.Bind.Time") {
        k <- 5
        start = c(log(c(beta.start, 0.1)), 1, 0.01)
        lower = c(log(bounds[1, c("beta", "alpha")]), bounds[1, c("a", "shift.timing")])
        upper = c(log(bounds[2, c("beta", "alpha")]), bounds[2, c("a", "shift.timing")])
        r.m <- split.vcv(tree, shift.time);
        foo <- function(x) {
          t.init <- rescale(tree, "EB", a = x[3])
          t <- rescale(t.init, "depth", depth=max(nodeHeights(tree)))
          release.mat1 <- split.vcv(t, x[4])
          release.mat2 <- split.vcv(tree, x[4])
          eb.mat <- release.mat1[[1]] * exp(x[1])
          bm2ou <- release.mat2[[2]] * exp(x[1])
          ou.mat <- ouMatrix(bm2ou, exp(x[2]))
          vv <- eb.mat + ou.mat
          mu <- phylogMean(vv, y)
          mu <- rep(mu, n)
          -dmvnorm(y, mu, vv, log = T)
        }
        o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
        ml.ouMatrix <- r.m[[2]] * exp(o$par[1])
        ml.vcv <- ((exp(o$par[1])) * (split.vcv((rescale(tree, "EB", a=o$par[3])), o$par[4])[[1]])) + ((ouMatrix(ml.ouMatrix, exp(o$par[2]))))
        #ml.vcv <- exp(o$par[1]) * ((ouMatrix(ml.ouMatrix, exp(o$par[2]))) + (r.m[[1]]))
        root.state <- phylogMean(ml.vcv, y)
        #root.state <- as.numeric(ml.root(tree = tree, model = model , y = y, meserr, params = o, shift.time));
        results <- list(lnl = -o$value, root.state = root.state, beta = exp(o$par[1]),  
                        alpha = exp(o$par[2]), a = o$par[3], shift.time = o$par[4]);
    }
    else if (model == "timeshift") {
    	k <- 3;
    	cat("shift point is at ", shift.time, "million years before present", "\n")
    	start = c(log(0.1), 1)
        lower = c(log(bounds[1, "beta"]), bounds[1, "nv"])
        upper = c(log(bounds[2, "beta"]), bounds[2, "nv"])
        foo <- function(x) {
            t <- timeshiftTree(phy = tree, breakPoint = shift.time, endRate = x[2]);
            vcv <- vcv(t)
            vv <- exp(x[1]) * vcv
            diag(vv) <- diag(vv) + meserr^2
            mu <- phylogMean(vv, y)
            mu <- rep(mu, n)
            -dmvnorm(y, mu, vv, log = T)
        }
   
        o <- optim(foo, p = start, lower = lower, upper = upper, method = "L")
           		root.state <- as.numeric(ml.root(tree = tree, model = model , y = y, meserr, params = o, shift.time));
		results <- list(lnl = -o$value, root.state = root.state, beta = exp(o$par[1]), 
            rate2 = exp(o$par[1]) * o$par[2])
     }
     
    else if (model == "release") {
    	k <- 3;
    	cat("shift point is at ", shift.time, "million years before present", "\n")
      	start = log(c(beta.start, 0.01))
        lower = log(bounds[1, c("beta", "alpha")])
        upper = log(bounds[2, c("beta", "alpha")])
        
        release.mat <- split.vcv(tree, shift.time);

        # foo <- function(x) {
            # t <- releaseTree(phy = tree, alpha = exp(x[2]), breakPoint = shift.time);
            # vcv <- vcv.phylo(t)
            # vv <- exp(x[1]) * vcv
            # diag(vv) <- diag(vv) + meserr^2
            # mu <- phylogMean(vv, y)
            # mu <- rep(mu, n)
            # -dmvnorm(y, mu, vv, log = T)
        # }
        foo <- function(x) {
          	ou.mat <- ouMatrix(release.mat[[1]], exp(x[2]))
          	bm.mat <- release.mat[[2]]
          
          	vv <- exp(x[1]) * (ou.mat + bm.mat)
          	diag(vv) <- diag(vv) + meserr^2
          	mu <- phylogMean(vv, y)
          	mu <- rep(mu, n)
          	-dmvnorm(y, mu, vv, log = T)
          }
        o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
        ml.vcv <- exp(o$par[1]) * ((ouMatrix(release.mat[[1]], exp(o$par[2]))) + release.mat[[2]])
        root.state <- phylogMean(ml.vcv, y)
        results <- list(lnl = -o$value, root.state = root.state, beta = exp(o$par[1]), 
            alpha = exp(o$par[2]));
     }
     
    else if (model == "releaseradiate") {
    	k <- 4;
    	cat("release point is at ", shift.time, "million years before present", "\n")
      	start = c(log(c(beta.start, 0.05)), 1)
        lower = c(log(bounds[1, c("beta", "alpha", "scalar")]))
        upper = c(log(bounds[2, c("beta", "alpha","scalar")]))
        release.mat <- split.vcv(tree, shift.time);
      foo <- function(x) {
			    ou.mat <- ouMatrix(release.mat[[1]], exp(x[2]))
			    bm.mat <- release.mat[[2]] * exp(x[3])
			    vv <- exp(x[1]) * (ou.mat + bm.mat) #same set-up as the 'o' object
			    diag(vv) <- diag(vv) + meserr^2
			    mu <- phylogMean(vv, y)
			    mu <- rep(mu, n)
			    -dmvnorm(y, mu, vv, log = T)
        }
        o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
        ml.vcv <- exp(o$par[1]) * ((ouMatrix(release.mat[[1]], exp(o$par[2]))) + (exp(o$par[3]) * release.mat[[2]]))
        root.state <- phylogMean(ml.vcv, y)
        results <- list(lnl = -o$value, root.state = root.state, beta = exp(o$par[1]),alpha = exp(o$par[2]), post.shift.scalar = exp(o$par[3]));
    } 
    
    else if (model == "radiate.constrain") {
      k <- 4;
      cat("release point is at ", shift.time, "million years before present", "\n")
      start = c(log(c(beta.start, 0.05)), 1)
      lower = c(log(bounds[1, c("beta", "alpha", "scalar")]))
      upper = c(log(bounds[2, c("beta", "alpha","scalar")]))
      release.mat <- split.vcv(tree, shift.time);
      foo <- function(x) {
        ou.mat <- ouMatrix(release.mat[[2]], exp(x[2]))
        bm.mat <- release.mat[[1]] * exp(x[3])
        vv <- exp(x[1]) * (ou.mat + bm.mat)
        diag(vv) <- diag(vv) + meserr^2
        mu <- phylogMean(vv, y)
        mu <- rep(mu, n)
        -dmvnorm(y, mu, vv, log = T)
      }
      o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
      ml.vcv <- exp(o$par[1]) * ((ouMatrix(release.mat[[2]], exp(o$par[2]))) + (exp(o$par[3]) * release.mat[[1]]))
      root.state <- phylogMean(ml.vcv, y)
      results <- list(lnl = -o$value, root.state = root.state, beta = exp(o$par[1]),alpha = exp(o$par[2]), post.shift.scalar = exp(o$par[3]));
    }
    
    else if (model == "SRC" || model == "BMOU") {
      k <- 3;
      cat("applying alpha parameter of OU process starting at", shift.time, "million years before present", "\n")
      start = log(c(beta.start, 0.01))
      lower = log(bounds[1, c("beta", "alpha")])
      upper = log(bounds[2, c("beta", "alpha")])
      release.mat <- split.vcv(tree, shift.time);
      
      # foo <- function(x) {
      # t <- releaseTree(phy = tree, alpha = exp(x[2]), breakPoint = shift.time);
      # vcv <- vcv.phylo(t)
      # vv <- exp(x[1]) * vcv
      # diag(vv) <- diag(vv) + meserr^2
      # mu <- phylogMean(vv, y)
      # mu <- rep(mu, n)
      # -dmvnorm(y, mu, vv, log = T)
      # }
      foo <- function(x) {
        ou.mat <- ouMatrix(release.mat[[2]], exp(x[2]))
        bm.mat <- release.mat[[1]]
        
        vv <- exp(x[1]) * (ou.mat + bm.mat)
        diag(vv) <- diag(vv) + meserr^2
        mu <- phylogMean(vv, y)
        mu <- rep(mu, n)
        -dmvnorm(y, mu, vv, log = T)
      }
      o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
      ml.vcv <- exp(o$par[1]) * ((ouMatrix(release.mat[[2]], exp(o$par[2]))) + release.mat[[1]])
      root.state <- phylogMean(ml.vcv, y)
      results <- list(lnl = -o$value, root.state = root.state, beta = exp(o$par[1]), 
                      alpha = exp(o$par[2]));
    }
    # SingleRateConstraint (BMOU) model that estimates the measurement error
    else if (model == "SRC_ME" || model == "BMOU_ME") {
      k <- 4;
      cat("applying alpha parameter of OU process starting at", shift.time, "million years before present", "\n")
      start = c(log(c(beta.start, 0.01)), 0.1)
      lower = c(log(bounds[1, c("beta", "alpha")]),bounds[1,"MErr"])
      upper = c(log(bounds[2, c("beta", "alpha")]),bounds[2,"MErr"])
      release.mat <- split.vcv(tree, shift.time);
      
      # foo <- function(x) {
      # t <- releaseTree(phy = tree, alpha = exp(x[2]), breakPoint = shift.time);
      # vcv <- vcv.phylo(t)
      # vv <- exp(x[1]) * vcv
      # diag(vv) <- diag(vv) + meserr^2
      # mu <- phylogMean(vv, y)
      # mu <- rep(mu, n)
      # -dmvnorm(y, mu, vv, log = T)
      # }
      foo <- function(x) {
        ou.mat <- ouMatrix(release.mat[[2]], exp(x[2]))
        bm.mat <- release.mat[[1]]
        
        vv <- exp(x[1]) * (ou.mat + bm.mat)
        diag(vv) <- diag(vv) + x[3]^2
        mu <- phylogMean(vv, y)
        #mu <- mean(y)
        mu <- rep(mu, n)
        -dmvnorm(y, mu, vv, log = T)
      }
      o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
      ml.vcv <- exp(o$par[1]) * ((ouMatrix(release.mat[[2]], exp(o$par[2]))) + release.mat[[1]])
      diag(ml.vcv) <- diag(ml.vcv) + o$par[3]^2
      root.state <- phylogMean(ml.vcv, y)
      #root.state <- mean(y)
      results <- list(lnl = -o$value, root.state = root.state, beta = exp(o$par[1]), 
                      alpha = exp(o$par[2]), MErr = o$par[3]);
    }
    
    # TwoRateConstraint (BMOUi) model that estimates the measurement error
    else if (model == "TRC_ME" || model == "BMOUi_ME") {
      k <- 5; #parameters are logLik, beta, alpha, post-shift-scalar
      cat("shift to OU process (sig.sq + alpha) begins at", shift.time, "million years before present", "\n")
      start = c(log(c(beta.start, 0.05)), 1, 0.1)
      lower = c(log(bounds[1, c("beta", "alpha", "scalar")]), bounds[1,"MErr"])
      upper = c(log(bounds[2, c("beta", "alpha", "scalar")]), bounds[2,"MErr"])
      #the above (upper/lower) are the values which will get optimized
      #here, x[1]=beta, x[2]=alpha, x[3]=scalar
      release.mat <- split.vcv(tree, shift.time);
      foo <- function(x) {
        bm2ou <- release.mat[[2]] * exp(x[3]) #apply the post-shift-scalar to the rate of the OU part of the tree/matrix
        ou.mat <- ouMatrix(bm2ou, exp(x[2])) #transform the second half of the tree by an OU process with alpha constraint
        bm.mat <- release.mat[[1]] #allow trait evolution to follow BM for the first part of the tree/matrix
        vv <- exp(x[1]) * (ou.mat + bm.mat) #combine the two together 
        diag(vv) <- diag(vv) + x[4]^2
        mu <- phylogMean(vv, y)
        mu <- rep(mu, n)
        -dmvnorm(y, mu, vv, log = T)
      }
      o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
      ml.ouMatrix <- release.mat[[2]] * exp(o$par[3])
      ml.vcv <- exp(o$par[1]) * ((ouMatrix(ml.ouMatrix, exp(o$par[2]))) + (release.mat[[1]]))
      diag(ml.vcv) <- diag(ml.vcv) + o$par[4]^2
      root.state <- phylogMean(ml.vcv, y)
      results <- list(lnl = -o$value, root.state = root.state, beta = exp(o$par[1]), alpha = exp(o$par[2]), 
                      post.shift.scalar = exp(o$par[3]), post.shift.sigma = exp(o$par[1])/exp(o$par[3]), ME = o$par[4]);
    }
    
    # correcting the "radiate.constrain" model (this applies the p-s-scalar, then transforms by alpha)
    else if (model == "TRC.time") {
      k <- 5; #parameters are logLik, beta, alpha, post-shift-scalar, shift.time
      cat("shift to OU process (sig.sq + alpha) begins at", shift.time, "million years before present", "\n")
      start = c(log(c(beta.start, 0.05)), 1, 7)
      lower = c(log(bounds[1, c("beta", "alpha", "scalar")]), bounds[1,"shift.timing"])
      upper = c(log(bounds[2, c("beta", "alpha", "scalar")]), bounds[2,"shift.timing"])
      #the above (upper/lower) are the values which will get optimized
      #here, x[1]=beta, x[2]=alpha, x[3]=scalar
      foo <- function(x) {
        release.mat <- split.vcv(tree, x[4]);
        bm2ou <- release.mat[[2]] * exp(x[3]) #apply the post-shift-scalar to the rate of the OU part of the tree/matrix
        ou.mat <- ouMatrix(bm2ou, exp(x[2])) #transform the second half of the tree by an OU process with alpha constraint
        bm.mat <- release.mat[[1]] #allow trait evolution to follow BM for the first part of the tree/matrix
        vv <- exp(x[1]) * (ou.mat + bm.mat) #combine the two together 
        diag(vv) <- diag(vv) + meserr^2
        mu <- phylogMean(vv, y)
        mu <- rep(mu, n)
        -dmvnorm(y, mu, vv, log = T)
      }
      o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
      release.matrix <- split.vcv(tree, o$par[4])
      ml.ouMatrix <- release.matrix[[2]] * exp(o$par[3])
      ml.vcv <- exp(o$par[1]) * ((ouMatrix(ml.ouMatrix, exp(o$par[2]))) + (release.matrix[[1]]))
      root.state <- phylogMean(ml.vcv, y)
      results <- list(lnl = -o$value, root.state = root.state, beta = exp(o$par[1]), alpha = exp(o$par[2]), post.shift.scalar = exp(o$par[3]), best.shift.time=o$par[4]);
    }
    
    # correcting the "radiate.constrain" model (this applies the p-s-scalar, then transforms by alpha)
    else if (model == "TRC" || model == "BMOUi") {
      k <- 4; #parameters are logLik, beta, alpha, post-shift-scalar
      cat("shift to OU process (sig.sq + alpha) begins at", shift.time, "million years before present", "\n")
      start = c(log(c(beta.start, 0.05)), 1)
      lower = c(log(bounds[1, c("beta", "alpha", "scalar")]))
      upper = c(log(bounds[2, c("beta", "alpha", "scalar")]))
      #the above (upper/lower) are the values which will get optimized
      #here, x[1]=beta, x[2]=alpha, x[3]=scalar
      release.mat <- split.vcv(tree, shift.time);
      foo <- function(x) {
        bm2ou <- release.mat[[2]] * exp(x[3]) #apply the post-shift-scalar to the rate of the OU part of the tree/matrix
        ou.mat <- ouMatrix(bm2ou, exp(x[2])) #transform the second half of the tree by an OU process with alpha constraint
        bm.mat <- release.mat[[1]] #allow trait evolution to follow BM for the first part of the tree/matrix
        vv <- exp(x[1]) * (ou.mat + bm.mat) #combine the two together 
        diag(vv) <- diag(vv) + meserr^2
        mu <- phylogMean(vv, y)
        mu <- rep(mu, n)
        -dmvnorm(y, mu, vv, log = T)
      }
      o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
      ml.ouMatrix <- release.mat[[2]] * exp(o$par[3])
      ml.vcv <- exp(o$par[1]) * ((ouMatrix(ml.ouMatrix, exp(o$par[2]))) + (release.mat[[1]]))
      root.state <- phylogMean(ml.vcv, y)
      results <- list(lnl = -o$value, root.state = root.state, beta = exp(o$par[1]), alpha = exp(o$par[2]), post.shift.scalar = exp(o$par[3]), post.shift.sigma = exp(o$par[1])/exp(o$par[3]));
    }
    
    # Alternative TRC, applies alpha, then scalar
    else if (model == "altTRC") {
      k <- 4; #parameters are logLik, beta, alpha, post-shift-scalar
      cat("shift to OU process (sig.sq + alpha) begins at", shift.time, "million years before present", "\n")
      start = c(log(c(beta.start, 0.05)), 1)
      lower = c(log(bounds[1, c("beta", "alpha", "scalar")]))
      upper = c(log(bounds[2, c("beta", "alpha", "scalar")]))
      #the above (upper/lower) are the values which will get optimized
      #here, x[1]=beta, x[2]=alpha, x[3]=scalar
      release.mat <- split.vcv(tree, shift.time);
      foo <- function(x) {
        bm2ou <- ouMatrix(release.mat[[2]], exp(x[2])) #transform the second half of the tree by an OU process with alpha constraint
        ou.mat <- bm2ou * exp(x[3]) #apply the post-shift-scalar to the rate of the OU part of the tree/matrix
        bm.mat <- release.mat[[1]] #allow trait evolution to follow BM for the first part of the tree/matrix
        vv <- exp(x[1]) * (ou.mat + bm.mat) #combine the two together 
        diag(vv) <- diag(vv) + meserr^2
        mu <- phylogMean(vv, y)
        mu <- rep(mu, n)
        -dmvnorm(y, mu, vv, log = T)
      }
      o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
      ml.ouMatrix <- release.mat[[2]] * exp(o$par[3])
      ml.vcv <- exp(o$par[1]) * ((ouMatrix(ml.ouMatrix, exp(o$par[2]))) + (release.mat[[1]]))
      root.state <- phylogMean(ml.vcv, y)
      results <- list(lnl = -o$value, root.state = root.state, beta = exp(o$par[1]), alpha = exp(o$par[2]), post.shift.scalar = exp(o$par[3]));
    }
    
    # trying a three-rate model (BM1-OU-BM2)
    else if (model == "BM1.OU.BM2") {
      k <- 5; #parameters are logLik, beta (sig.sq BM1), alpha, post-shift-scalar (OU), post-shift-scalar (BM2)
      cat("evolving under a BM process until shift to OU process (sig.sq + alpha) begins at", shift.time, "million years before present",
          "it then reverts to a BM process at", shift.time2, "million years before present", "\n")
      start = c(log(c(beta.start, 0.05)), 1,1)
      lower = c(log(bounds[1, c("beta", "alpha", "scalar", "scalar")]))
      upper = c(log(bounds[2, c("beta", "alpha", "scalar", "scalar")]))
      #the above (upper/lower) are the values which will get optimized
      #here, x[1]=beta, x[2]=alpha, x[3]=scalar
      release.mat <- split.3.vcv(tree, shift.time, shift.time2);
      foo <- function(x) {
        bm2ou <- release.mat[[2]] * exp(x[3]) #apply the post-shift-scalar to the rate of the OU part of the tree/matrix
        ou.mat <- ouMatrix(bm2ou, exp(x[2])) #transform the second half of the tree by an OU process with alpha constraint
        bm.mat1 <- release.mat[[1]] #allow trait evolution to follow BM for the first part of the tree/matrix
        bm.mat2 <- release.mat[[3]] * exp(x[4]) #allow trait evolution to follow BM for the third part of the tree/matrix
        vv <- exp(x[1]) * (ou.mat + bm.mat1 + bm.mat2) #combine the three matrices together, multiply by the rate
        diag(vv) <- diag(vv) + meserr^2
        mu <- phylogMean(vv, y)
        mu <- rep(mu, n)
        -dmvnorm(y, mu, vv, log = T)
      }
      o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
      ml.ouMatrix <- release.mat[[2]] * exp(o$par[3])
      ml.vcv <- exp(o$par[1]) * ((ouMatrix(ml.ouMatrix, exp(o$par[2]))) + (release.mat[[1]]) + (release.mat[[3]]))
      root.state <- phylogMean(ml.vcv, y)
      results <- list(lnl = -o$value, root.state = root.state, beta = exp(o$par[1]), alpha = exp(o$par[2]), post.shift.scalar = exp(o$par[3]), scalar2 = exp(o$par[4]));
    }
    
    # trying a three-rate model (BM1-OU-BM1)
    else if (model == "BM1.OU.BM1") {
      k <- 4; #parameters are logLik, beta (sig.sq BM1), alpha, post-shift-scalar (OU)
      cat("evolving under a BM process until shift to OU process (sig.sq + alpha) begins at", shift.time, "million years before present",
          "it then reverts to a BM process at", shift.time2, "million years before present", "\n")
      start = c(log(c(beta.start, 0.05)), 1,1)
      lower = c(log(bounds[1, c("beta", "alpha", "scalar", "scalar")]))
      upper = c(log(bounds[2, c("beta", "alpha", "scalar", "scalar")]))
      #the above (upper/lower) are the values which will get optimized
      #here, x[1]=beta, x[2]=alpha, x[3]=scalar
      release.mat <- split.3.vcv(tree, shift.time, shift.time2);
      foo <- function(x) {
        bm2ou <- release.mat[[2]] * exp(x[3]) #apply the post-shift-scalar to the rate of the OU part of the tree/matrix
        ou.mat <- ouMatrix(bm2ou, exp(x[2])) #transform the second half of the tree by an OU process with alpha constraint
        bm.mat1 <- release.mat[[1]] #allow trait evolution to follow BM for the first part of the tree/matrix
        bm.mat2 <- release.mat[[3]] #allow trait evolution to follow BM for the third part of the tree/matrix
        vv <- exp(x[1]) * (ou.mat + bm.mat1 + bm.mat2) #combine the three matrices together, multiply by the rate
        diag(vv) <- diag(vv) + meserr^2
        mu <- phylogMean(vv, y)
        mu <- rep(mu, n)
        -dmvnorm(y, mu, vv, log = T)
      }
      o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
      ml.ouMatrix <- release.mat[[2]] * exp(o$par[3])
      ml.vcv <- exp(o$par[1]) * ((ouMatrix(ml.ouMatrix, exp(o$par[2]))) + (release.mat[[1]]) + (release.mat[[3]]))
      root.state <- phylogMean(ml.vcv, y)
      results <- list(lnl = -o$value, root.state = root.state, beta = exp(o$par[1]), alpha = exp(o$par[2]), post.shift.scalar = exp(o$par[3]));
    }
    
    # trying a three-rate model (BM1-OU-BM1)
    else if (model == "BM1.OU.BM1.time") {
      k <- 6; #parameters are logLik, beta (sig.sq BM1), alpha, post-shift-scalar (OU)
      cat("evolving under a BM process until shift to OU process (sig.sq + alpha) begins at", shift.time, "million years before present",
          "it then reverts to a BM process at", shift.time2, "million years before present", "\n")
      start = c(log(c(beta.start, 0.05)), 1, 1, 10, 10)
      lower = c(log(bounds[1, c("beta", "alpha", "scalar", "scalar")]), bounds[1,c("shift.timing", "shift.timing")])
      upper = c(log(bounds[2, c("beta", "alpha", "scalar", "scalar")]), bounds[2,c("shift.timing", "shift.timing")])
      #the above (upper/lower) are the values which will get optimized
      #here, x[1]=beta, x[2]=alpha, x[3]=scalar
      foo <- function(x) {
        release.mat <- split.3.vcv(tree, x[5], x[6]);
        bm2ou <- release.mat[[2]] * exp(x[3]) #apply the post-shift-scalar to the rate of the OU part of the tree/matrix
        ou.mat <- ouMatrix(bm2ou, exp(x[2])) #transform the second half of the tree by an OU process with alpha constraint
        bm.mat1 <- release.mat[[1]] #allow trait evolution to follow BM for the first part of the tree/matrix
        bm.mat2 <- release.mat[[3]] #allow trait evolution to follow BM for the third part of the tree/matrix
        vv <- exp(x[1]) * (ou.mat + bm.mat1 + bm.mat2) #combine the three matrices together, multiply by the rate
        diag(vv) <- diag(vv) + meserr^2
        mu <- phylogMean(vv, y)
        mu <- rep(mu, n)
        -dmvnorm(y, mu, vv, log = T)
      }
      o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
      release.matrix <- split.3.vcv(tree, o$par[5], o$par[6]);
      ml.ouMatrix <- release.matrix[[2]] * exp(o$par[3])
      ml.vcv <- exp(o$par[1]) * ((ouMatrix(ml.ouMatrix, exp(o$par[2]))) + (release.matrix[[1]]) + (release.matrix[[3]]))
      root.state <- phylogMean(ml.vcv, y)
      results <- list(lnl = -o$value, root.state = root.state, beta = exp(o$par[1]), alpha = exp(o$par[2]), 
                      post.shift.scalar = exp(o$par[3]), shift.timing1 = o$par[5], shift.timing2 = o$par[6]);
    }
    
    
    results$aic <- 2 * k - 2 * results$lnl
    results$aicc <- 2 * k * (n - 1)/(n - k - 2) - 2 * results$lnl
    results$k <- k
    return(results)
}

#############################################################################################################################


phylogMean <- function (phyvcv, data) 
{
    o <- rep(1, length(data))
    ci <- solve(phyvcv)
    m1 <- solve(t(o) %*% ci %*% o)
    m2 <- t(o) %*% ci %*% data
    return(m1 %*% m2)
}

#############################################################################################################################

BranchingTimesFossil <- function (phy, tol = .Machine$double.eps^0.5) 
{
    if (class(phy) != "phylo") 
        stop("object \"phy\" is not of class \"phylo\"")
    phy2 <- phy
    phy <- new2old.phylo(phy)
    tmp <- as.numeric(phy$edge)
    nb.tip <- max(tmp)
    nb.node <- -min(tmp)
    xx <- as.numeric(rep(NA, nb.tip + nb.node))
    names(xx) <- as.character(c(-(1:nb.node), 1:nb.tip))
    xx["-1"] <- 0
    for (i in 2:length(xx)) {
        nod <- names(xx[i])
        ind <- which(as.numeric(phy$edge[, 2]) == nod)
        base <- phy$edge[ind, 1]
        xx[i] <- xx[base] + phy$edge.length[ind]
    }
    bt <- abs(xx - max(xx));
	
	for(i in 1:length(bt)) {
		
		if(bt[i]<.Machine$double.eps^0.5) bt[i] <- 0; 	}
	
	names(bt) <- c(seq(nb.tip+1, nb.tip+nb.node), phy$tip.label)
	
	
	return(bt);
}

#############################################################################################################################
ml.root <- function(tree, model, y, meserr, params, shift.time, shift.time2) {
	
	if(model == "BM") {
		vcv <- vcv.phylo(tree)
		vv <- exp(params$par) * vcv;
		diag(vv) <- diag(vv) + meserr^2;
		return(phylogMean(vv, y));
	} 
	if(model == "OU") {
		vcvOrig <- vcv.phylo(tree);
		vcv <- ouMatrix(vcvOrig, exp(params$par[2]));
		vv <- exp(params$par[1]) * vcv;
		diag(vv) <- diag(vv) + meserr^2;
		return(phylogMean(vv, y));
	}
	#if(model == "EB") {
	#	t <- transform(tree, "EB", a = params$par[2]);
  #       vcv <- vcv.phylo(t);
  #         vv <- exp(params$par[1]) * vcv;
  #         diag(vv) <- diag(vv) + meserr^2;
  #          return(phylogMean(vv, y))
	#}
  #if(model == "Burst.Bind") {
  #  t <- transform(tree, "EB", a = params$par[4]);
  #  vcv <- vcv.phylo(t);
  #  vv <- exp(params$par[1]) * vcv;
  #  diag(vv) <- diag(vv) + meserr^2;
  #  return(phylogMean(vv, y))
  #}
	if(model == "timeshift") {
            t <- timeshiftTree(phy = tree, breakPoint = shift.time, endRate = params$par[2]);
            vcv <- vcv.phylo(t)
            vv <- exp(params$par[1]) * vcv
            diag(vv) <- diag(vv) + meserr^2
            return(phylogMean(vv, y));	
	} else if(model == "release") {
            t <- releaseTree(phy = tree, alpha = exp(params$par[2]), breakPoint = shift.time);
            vcv <- vcv.phylo(t)
            vv <- exp(params$par[1]) * vcv
            diag(vv) <- diag(vv) + meserr^2
           	return(phylogMean(vv, y));
	} else if(model == "releaseradiate") {
            t <- releaseTree(phy = tree, alpha = exp(params$par[2]), breakPoint = shift.time, endRate = params$par[3]);
            vcv <- vcv.phylo(t)
            vv <- exp(params$par[1]) * vcv
            diag(vv) <- diag(vv) + meserr^2
            mu <- phylogMean(vv, y);
	} else if(model == "SRC" || model == "SRC_ME" || model == "BMOU" || model = "BMOU_ME") {
	          t <- releaseTree(phy = tree, alpha = exp(params$par[2]), breakPoint = shift.time);
	          vcv <- vcv.phylo(t)
	          vv <- exp(params$par[1]) * vcv
	          diag(vv) <- diag(vv) + meserr^2
	          return(phylogMean(vv, y));
	} else if(model == "Burst.Bind") {
	          t <- releaseTree(phy = tree, alpha = exp(params$par[2]), breakPoint = shift.time);
	          vcv <- vcv.phylo(t)
	          vv <- exp(params$par[1]) * vcv
	          diag(vv) <- diag(vv) + meserr^2
	          return(phylogMean(vv, y));
	} else if(model == "Burst.Bind.Time") {
	          t <- releaseTree(phy = tree, alpha = exp(params$par[2]), breakPoint = params$par[4]);
	          vcv <- vcv.phylo(t)
	          vv <- exp(params$par[1]) * vcv
	          diag(vv) <- diag(vv) + meserr^2
	          return(phylogMean(vv, y));
	} else if(model == "radiate.constrain") {
	          t <- releaseTree(phy = tree, alpha = exp(params$par[2]), breakPoint = shift.time, endRate = params$par[3]);
	          vcv <- vcv.phylo(t)
	          vv <- exp(params$par[1]) * vcv
	          diag(vv) <- diag(vv) + meserr^2
	          mu <- phylogMean(vv, y)
	} else if(model == "TRC" || model == "TRC_ME" || model == "BMOUi" || model == "BMOUi_ME") {
	          t <- releaseTree(phy = tree, alpha = exp(params$par[2]), breakPoint = shift.time, endRate = params$par[3]);
	          vcv <- vcv.phylo(t)
	          vv <- exp(params$par[1]) * vcv
	          diag(vv) <- diag(vv) + meserr^2
	          mu <- phylogMean(vv, y)
	} else if(model == "TRC.time") {
	          t <- releaseTree(phy = tree, alpha = exp(params$par[2]), breakPoint = shift.time, endRate = params$par[3]);
	          vcv <- vcv.phylo(t)
	          vv <- exp(params$par[1]) * vcv
	          diag(vv) <- diag(vv) + meserr^2
	          mu <- phylogMean(vv, y)
	} else if(model == "BM1.OU.BM2") {
	          t <- releaseTree(phy = tree, alpha = exp(params$par[2]), breakPoint = shift.time, endRate = params$par[3]);
	          vcv <- vcv.phylo(t)
	          vv <- exp(params$par[1]) * vcv
	          diag(vv) <- diag(vv) + meserr^2
	          mu <- phylogMean(vv, y)
	} else if(model == "BM1.OU.BM1") {
	  t <- releaseTree(phy = tree, alpha = exp(params$par[2]), breakPoint = shift.time, endRate = params$par[3]);
	  vcv <- vcv.phylo(t)
	  vv <- exp(params$par[1]) * vcv
	  diag(vv) <- diag(vv) + meserr^2
	  mu <- phylogMean(vv, y)
	} else if(model == "BM1.OU.BM1.time") {
	  t <- releaseTree(phy = tree, alpha = exp(params$par[2]), breakPoint = shift.time, endRate = params$par[3]);
	  vcv <- vcv.phylo(t)
	  vv <- exp(params$par[1]) * vcv
	  diag(vv) <- diag(vv) + meserr^2
	  mu <- phylogMean(vv, y)
	}
}## end function

#########################################################################################
ouMatrix <- function (vcvMatrix, alpha) 
{
    vcvDiag <- diag(vcvMatrix)
    diagi <- matrix(vcvDiag, nrow = length(vcvDiag), ncol = length(vcvDiag))
    diagj <- matrix(vcvDiag, nrow = length(vcvDiag), ncol = length(vcvDiag), 
        byrow = T)
    Tij = diagi + diagj - (2 * vcvMatrix)
    vcvRescaled = (1/(2 * alpha)) * exp(-alpha * Tij) * (1 - 
        exp(-2 * alpha * vcvMatrix))
    return(vcvRescaled)
}

#########################################################################################

ACDC.prior <- function (phy, decrease.max = 1e-05, increase.max = 1e+05) 
{
    if (is.ultrametric(phy)) {
        max.bt <- max(branching.times(phy))
    }
    else {
        max.bt <- max(BranchingTimesFossil(phy))
    }
    prior.min <- log(decrease.max)/max.bt
    prior.max <- log(increase.max)/max.bt
    return(list(min = prior.min, max = prior.max))
}
#########################################################################################
mrca.of.pair <- function (phy, tip1, tip2) 
{
    if (is.na(match(tip1, phy$tip.label))) 
        stop(paste(tip1, " is not a valid taxon in this tree"))
    if (is.null(tip2)) {
        bb <- which(phy$tip.label == tip1)
        mrca <- parent.node(phy, bb)
    }
    else {
        if (is.na(match(tip2, phy$tip.label))) 
            stop(paste(tip2, " is not a valid taxon in this tree"))
        nn <- phy$Nnode
        nt <- length(phy$tip.label)
        minLeaves <- nt
        mrca <- NULL
        for (i in 1:nn) {
            leaves <- node.leaves(phy, i + nt)
            if (tip1 %in% leaves & tip2 %in% leaves) {
                ll <- length(leaves)
                if (ll < minLeaves) {
                  mrca <- i + nt
                  minLeaves <- ll
                }
            }
        }
        if (is.null(mrca)) 
            mrca <- length(phy$tip.label) + 1
    }
    return(mrca)
}


#### write a function to obtain the variance covariance matrix for an OU process with non-ultrametric trees. Follows Eq. 2 and Eq. 6 in Hansen 1997

ou.fossil.vcv <- function(distance, vcv.mat, alpha, sigma_sq) {

	## following Hansen 1997
	## first need the phylogenetic distance matrix between tips in units of time
	
	N <- nrow(distance)
	## and then the covariance (root to tip distance) in units of time
	taxa <- rownames(vcv.mat)
	times <- diag(vcv.mat) ## gets the root to tip distances for later
	
	# following Hansen, we need the exponent of -alpha* the time separating pairs of species -- the exponential decay with separation time (Hansen 1997 p1344)
	 
	dist.conv <- exp(-alpha * distance)

	# we also need the variance of the common ancestor
	
	cov.conv <- 1-exp(-2*alpha*vcv.mat)
	
	## now take the product of these two and the equilibrium variance (sigmasq / 2* alpha) 
	vcv.ou <- (sigma_sq/(2*alpha)) *(cov.conv* as.vector(dist.conv))
	
	## finally we need to replace the diagonals 
	diag(vcv.ou) <- (sigma_sq / (2*alpha)) * (1 - exp(-2*alpha * times ))
	rownames(vcv.ou) <- colnames(vcv.ou)  <-taxa
	## and return..
	return(vcv.ou)
	}
	

####
split.vcv <- function(phy, time) {
	
	mat <- vcv(phy)
	mat1 <- mat
	mat2 <- mat
	n <- nrow(mat)
	root <- max(mat)
	shift.from.root <- root - time

make.mat.1 <- function(x, shift.from.root) {
		if(x == 0) {
			return(0)
		}
		if(x<shift.from.root) {
			return(x)
		} 
		if(x > shift.from.root){
			return(shift.from.root)
		}
}

make.mat.2 <- function(x, time, shift.from.root) {
		if(x == 0) {
			return(0)
		}
		
		if(x < shift.from.root) {
				return(0)
				} else{
					return(x-shift.from.root)
					}
}


mat1 <- matrix(sapply(mat1, make.mat.1, shift.from.root =  shift.from.root), nrow = n, ncol = n, byrow=T)

diag1 <- diag(mat);
diag.foo <- function(x) {
		if(x<time) {
			return(x)
		} 
		if(x > time){
		
		} 
}
mat2 <- matrix(sapply(mat2, make.mat.2, time =  time, shift.from.root = shift.from.root), nrow = n, ncol = n, byrow=T)

rownames(mat1) <- rownames(mat2) <- rownames(mat)
colnames(mat1) <- colnames(mat2) <- colnames(mat)

return(list(mat1 = mat1, mat2 = mat2))

}



####
split.3.vcv <- function(phy, time1, time2) {
  
  mat <- vcv(phy)
  mat1 <- mat
  mat2 <- mat
  mat3 <- mat
  n <- nrow(mat)
  root <- max(mat)
  shift1.from.root <- root - time1
  shift2.from.root <- root - time2
  
make.mat.1 <- function(x, shift1.from.root) {
    if(x == 0) { #this makes all 0 = 0
      return(0)
    }
    if(x<shift1.from.root) { #this keeps all values that are below s1.f.r (older than the shift1)
      return(x)
    } 
    if(x > shift1.from.root){ #this makes all values above s1.f.r = s1.f.r (equal to the shift1 age)
      return(shift1.from.root)
    }
}
  
make.mat.2 <- function(x, time1, time2, shift1.from.root, shift2.from.root) {
    if(x == 0) { #this keeps all 0 = 0
      return(0)
    }
    if(x > shift2.from.root) {
      return(shift2.from.root-shift1.from.root)
    }
    if((shift2.from.root > x) & (x > shift1.from.root)) { #this makes all values that are above s2.f.r = 0 (makes younger than shift2 = 0)
    return(x-shift1.from.root)
    }else{
      return(0)
    }
}
  
make.mat.3 <- function(x, time2, shift2.from.root) {
    if(x == 0) {
      return(0)
    }
    
    if(x < shift2.from.root) {
      return(0)
    } 
    else{
      return(x-shift2.from.root)
    }
}
  
mat1 <- matrix(sapply(mat1, make.mat.1, shift1.from.root =  shift1.from.root), nrow = n, ncol = n, byrow=T)
  
diag1 <- diag(mat);
diag.foo <- function(x) {
    if(x<time) {
      return(x)
    } 
    if(x > time){
      
    } 
}
mat2 <- matrix(sapply(mat2, make.mat.2, time1 =  time1, time2 = time2, shift1.from.root = shift1.from.root, shift2.from.root = shift2.from.root), nrow = n, ncol = n, byrow=T)
mat3 <- matrix(sapply(mat3, make.mat.3, time2 =  time2, shift2.from.root = shift2.from.root), nrow = n, ncol = n, byrow=T)
  
rownames(mat1) <- rownames(mat2) <- rownames(mat3) <- rownames(mat)
colnames(mat1) <- colnames(mat2) <- rownames(mat3) <- colnames(mat)
  
return(list(mat1 = mat1, mat2 = mat2, mat3 = mat3))
  
}
####################################################################
chop.tree <- function (tr, timepoint = 10, return_pieces = TRUE) 
{
  tr_table = prt(tr, printflag = FALSE, get_tipnames = FALSE)
  tr_table
  TF_exists_more_recently_than_10mya = tr_table$time_bp < timepoint
  labels_for_tips_existing_more_recently_than_10mya = tr_table$label[TF_exists_more_recently_than_10mya == 
                                                                       TRUE]
  edge_times_bp = get_edge_times_before_present(tr)
  edges_start_earlier_than_10mya = edge_times_bp[, 1] > timepoint
  edges_end_later_than_10mya = edge_times_bp[, 2] <= timepoint
  edges_to_chainsaw = edges_start_earlier_than_10mya + edges_end_later_than_10mya == 
    2
  nodes_to_chainsaw = tr$edge[, 2][edges_to_chainsaw]
  numtips = length(tr$tip.label)
  tree_to_chainsaw = tr
  if (return_pieces == TRUE) {
    return_pieces_list = as.list(rep(NA, length(nodes_to_chainsaw)))
    return_pieces_basenames = as.list(rep(NA, length(nodes_to_chainsaw)))
    chopTable = NULL
  }
  chainsaw_table = NULL
  for (i in 1:length(nodes_to_chainsaw)) {
    if (nodes_to_chainsaw[i] <= numtips) {
      if (return_pieces == TRUE) {
        return_pieces_list[[i]] = timepoint
        tmp_tipname = tr$tip.label[nodes_to_chainsaw[i]]
        return_pieces_basenames[[i]] = tmp_tipname
      }
    }
    #else{
    #  nodes_to_chainsaw[i] <- paste("Node", nodes_to_chainsaw[i], sep=".")
    #}
    
    else {
      tmp_subtree = extract.clade(tr, nodes_to_chainsaw[i])
      branchlength_below_subtree_LCA_node = timepoint - 
        get_max_height_tree(tmp_subtree)
      tmp_subtree$root.edge = branchlength_below_subtree_LCA_node
      if (return_pieces == TRUE) {
        return_pieces_list[[i]] = tmp_subtree
        tmp_labels_merge = paste(tmp_subtree$tip.label, 
                                 collapse = ",", sep = "")
        tmp_labels_split = strsplit(tmp_labels_merge, 
                                    split = ",")[[1]]
        new_labels = sort(tmp_labels_split)
        basename_after_cutting = paste(new_labels, collapse = ",", 
                                       sep = "")
        return_pieces_basenames[[i]] = basename_after_cutting
      }
      tmp_number_of_tips = length(tmp_subtree$tip.label)
      numtips_to_drop = tmp_number_of_tips - 1
      tmp_labels = tmp_subtree$tip.label
      labels_to_drop = tmp_labels[1:numtips_to_drop]
      ordered_labels_to_make_into_new_name = sort(tmp_labels)
      #name_new_tip = paste(ordered_labels_to_make_into_new_name, 
      #                     collapse = ",", sep = "")
      name_new_tip = paste("Node", nodes_to_chainsaw[i], sep=".")
      label_kept_num = length(tmp_labels)
      label_kept = tmp_labels[label_kept_num]
      new_label = name_new_tip
      tree_to_chainsaw$tip.label[tree_to_chainsaw$tip.label == 
                                   label_kept] = new_label
      tree_to_chainsaw = drop.tip(tree_to_chainsaw, labels_to_drop)
    }
  }
  tree_to_chainsaw_table = prt(tree_to_chainsaw, printflag = FALSE)
  tree_to_chainsaw_table_tips_TF_time_bp_LT_10my = tree_to_chainsaw_table$time_bp < 
    timepoint
  tmp_edge_lengths = tree_to_chainsaw_table$edge.length[tree_to_chainsaw_table_tips_TF_time_bp_LT_10my]
  times_bp_for_edges_to_chainsaw = tree_to_chainsaw_table$time_bp[tree_to_chainsaw_table_tips_TF_time_bp_LT_10my]
  adjustment = times_bp_for_edges_to_chainsaw - timepoint
  revised_tmp_edge_lengths = tmp_edge_lengths + adjustment
  tree_to_chainsaw_table$edge.length[tree_to_chainsaw_table_tips_TF_time_bp_LT_10my] = revised_tmp_edge_lengths
  ordered_nodenames = get_nodenums(tree_to_chainsaw)
  parent_branches = get_indices_where_list1_occurs_in_list2(ordered_nodenames, 
                                                            tree_to_chainsaw$edge[, 2])
  NA_false = is.not.na(tree_to_chainsaw_table$edge.length)
  tree_to_chainsaw$edge.length[parent_branches[NA_false]] = tree_to_chainsaw_table$edge.length[NA_false]
  if (return_pieces == TRUE) {
    chainsaw_result = NULL
    chainsaw_result$tree_to_chainsaw = tree_to_chainsaw
    chainsaw_result$return_pieces_list = return_pieces_list
    chainsaw_result$return_pieces_basenames = return_pieces_basenames
    class(chainsaw_result) = "chainsaw_result"
    return(chainsaw_result)
  }
  else {
    return(tree_to_chainsaw)
  }
}
