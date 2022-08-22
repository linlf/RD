## Data should be in a long format with columns for study identifier,
##   treatment identifier, events and total arm sample size

## Loading the "lme4" package for implementing (generalized) linear mixed-effects models
library("lme4")

## Function for avoiding bootstrap resampled meta-analyses with warnings and errors
##   when implementing the bivariate random-effects model for deriving the CI
tryCatch.W.E <- function(expr){
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
  warning = w.handler), warning = W)
}

## The logit function
logit <- function(x){
  if(x < 0.001) x <- 0.001
  if(x > 0.999) x <- 0.999
  out <- log(x/(1 - x))
  return(out)
}
logit <- Vectorize(logit)

## The expit function (the inverse function of the logit)
expit <- function(x) 1/(1 + exp(-x))
expit <- Vectorize(expit)


## The function for estimating the RD at different baseline risks.
## It was used to generate Figure 4 in the main content.
## The input arguments for this function are as follows:
##   sid and tid: the study and treatment IDs, respectively;
##   e and n: the event count and sample size, respectively;
##   data: the data object;
##   link: the link function for the bivariate random-effects model,
##     with the logit link as the default;
##   alpha: the significance level (the default is 0.05), leading to (1 - alpha)x100% CI;
##   b.iter: the number of iterations for the bootstrap resampling to derive the CI;
##   from and to: the range of baseline risks in the generated plot;
##   ylim: (optional) the range of the y-axis (risk difference) in the generated plot;
##   p0: (optional) the baseline risks specified by users at which
##     the numeric RD estimates with CIs will be recorded;
##   seed: the random number generator for reproducing the bootstrap resampling.
## The output from this function includes:
##   a plot showing the risk difference against the baseline risk;
##   parameter estimates (parameter.estimate) consisting of mu0, mu1, sd0, sd1, and rho,
##     representing the overall means of logit risks, between-study standard deviations
##     in the two groups, and the correlation coefficient between logit risk
##     in group 1 and logit risk in group 2;
##   numbers of warnings, errors, and not a number (NaN) during
##     the bootstrap resampling (bootstrap.message);
##   estimates of risk difference (RD.est) at the values assigned to p0 (if given).
RD.cond <- function(sid, tid, e, n, data, link = "logit",
  alpha = 0.05, b.iter = 1000, from = 0, to = 1, ylim, p0, seed = 1234){
  ## assessing the validity of the input objects
  if(missing(sid)) stop("need to specify study IDs.")
  if(missing(tid)) stop("need to specify treatment IDs (0/1).")
  if(missing(e)) stop("need to specify event counts.")
  if(missing(n)) stop("need to specify sample sizes.")
  if(!missing(data)){
    sid <- eval(substitute(sid), data, parent.frame())
    tid <- eval(substitute(tid), data, parent.frame())
    e <- eval(substitute(e), data, parent.frame())
    n <- eval(substitute(n), data, parent.frame())
  }
  if((length(sid) != length(n)) | (length(tid) != length(n)) | (length(e) != length(n))){
    stop("the input data have different dimensions.")
  }
  if(alpha < 0 | alpha > 1) stop("alpha should be between 0 and 1.")
  if(!missing(p0)){
    if(any(p0 < 0) | any(p0 > 1)) stop("p0 must be between 0 and 1.")
  }

  ## implementing the bivariate random-effects model
  ord <- order(sid)
  sid <- sid[ord]
  tid <- tid[ord]
  e <- e[ord]
  n <- n[ord]
  data <- data.frame(sid = sid, tid = tid, cid = 1 - tid, e = e, n = n)
  rslt <- glmer(cbind(e, n - e) ~ factor(tid) + (cid + tid - 1 | sid),
    data = data, family = binomial(link = link))
  smry <- summary(rslt)
  mu0 <- smry$coefficients[1,1]
  mu1 <- smry$coefficients[2,1] + mu0
  sd0 <- attr(smry$varcor[[1]], "stddev")[1]
  sd1 <- attr(smry$varcor[[1]], "stddev")[2]
  rho <- attr(smry$varcor[[1]], "correlation")[1,2]
  names(sd0) <- names(sd1) <- NULL
  param.est <- c(mu0 = mu0, mu1 = mu1, sd0 = sd0, sd1 = sd1, rho = rho)

  RD <- function(mu0, mu1, sd0, sd1, rho, p0){
    out <- expit(mu1 + rho*sd1/sd0*(logit(p0) - mu0)) - p0
    return(out)
  }
  RD <- Vectorize(RD)
  RD.est <- function(p0){RD(mu0, mu1, sd0, sd1, rho, p0)}

  ## performing the bootstrap resampling for deriving the CI
  set.seed(seed)
  b.idx <- 0
  b.w <- b.e <- b.nan <- 0
  mu0.b <- mu1.b <- sd0.b <- sd1.b <- rho.b <- rep(NA, b.iter)
  while(b.idx < b.iter){
    idx.temp <- sample(length(unique(sid)), replace = TRUE)
    idx.temp <- c(rbind(idx.temp*2 - 1, idx.temp*2))
    data.temp <- data[idx.temp,]
    data.temp$sid <- rep(1:length(unique(sid)), each = 2)
    suppressMessages(out.W.E <- tryCatch.W.E(
      rslt.temp <- glmer(cbind(e, n - e) ~ factor(tid) + (cid + tid - 1 | sid),
        data = data.temp, family = binomial(link = "logit"))))
    suppressMessages(out.W.E2 <- tryCatch.W.E(smry.temp <- summary(rslt.temp)))
    if(is.null(out.W.E$warning) & !inherits(out.W.E$value, "error") &
      is.null(out.W.E2$warning) & !inherits(out.W.E2$value, "error")){
      mu0.b.temp <- smry.temp$coefficients[1,1]
      mu1.b.temp <- smry.temp$coefficients[2,1] + mu0.b.temp
      sd0.b.temp <- attr(smry.temp$varcor[[1]], "stddev")[1]
      sd1.b.temp <- attr(smry.temp$varcor[[1]], "stddev")[2]
      rho.b.temp <- attr(smry.temp$varcor[[1]], "correlation")[1,2]
      if(!is.null(out.W.E$warning) | !is.null(out.W.E2$warning)){
        b.w <- b.w + 1
      }
      if(inherits(out.W.E$value, "error") | inherits(out.W.E2$value, "error")){
        b.e <- b.e + 1
      }
      if(!is.nan(mu0.b.temp) & !is.nan(mu1.b.temp) &
        !is.nan(sd0.b.temp) & !is.nan(sd1.b.temp) & !is.nan(rho.b.temp)){
        b.idx <- b.idx + 1
        mu0.b[b.idx] <- mu0.b.temp
        mu1.b[b.idx] <- mu1.b.temp
        sd0.b[b.idx] <- sd0.b.temp
        sd1.b[b.idx] <- sd1.b.temp
        if(rho.b.temp < -0.999) rho.b.temp <- -0.999
        if(rho.b.temp > 0.999) rho.b.temp <- 0.999
        rho.b[b.idx] <- rho.b.temp
      }else{
        b.nan <- b.nan + 1
      }
    }
  }
  b.msg <- c("Warning" = b.w, "Error" = b.e, "NaN" = b.nan)

  ## deriving the CI based on the bootstrap resampled meta-analyses
  RD.ci <- function(p0){
    RD.b <- RD(mu0.b, mu1.b, sd0.b, sd1.b, rho.b, p0)
    lb <- quantile(RD.b, probs = 0.025, na.rm = TRUE)
    ub <- quantile(RD.b, probs = 0.975, na.rm = TRUE)
    names(lb) <- names(ub) <- NULL
    out <- list(lb = lb, ub = ub)
    return(out)
  }

  RD.ci.lb <- function(p0){RD.ci(p0)$lb}
  RD.ci.ub <- function(p0){RD.ci(p0)$ub}
  RD.ci.lb <- Vectorize(RD.ci.lb)
  RD.ci.ub <- Vectorize(RD.ci.ub)

  ## generating the plot
  xx <- seq(max(c(from - (to - from)/10, 0)), min(c(to + (to - from)/10, 1)), length.out = 1000)
  yy1 <- RD.ci.lb(xx)
  yy2 <- RD.ci.ub(xx)
  if(missing(ylim)) ylim <- c(min(yy1), max(yy2))
  plot(1, type = "n", xlab = "Baseline risk", ylab = "Risk difference",
    xlim = c(from, to), ylim = ylim)
  polygon(c(xx, rev(xx)), c(yy1, rev(yy2)), border = NA,
    col = adjustcolor("blue", alpha.f = 0.4))
  curve(RD.est, from = max(c(from - (to - from)/10, 0)), to = min(c(to + (to - from)/10, 1)),
    lwd = 2, add = TRUE)
  risk <- e/n
  risk.b <- risk[tid == 0]
  risk.t <- risk[tid == 1]
  RD.studies <- risk.t - risk.b
  points(risk.b, RD.studies, pch = 16, col = adjustcolor("black", alpha.f = 0.4))

  ## organizing the output
  out <- list(parameter.estimate = param.est, bootstrap.message = b.msg)
  if(!missing(p0)){
    out$RD.est <- cbind(p0 = p0, RD.est = RD.est(p0),
      RD.CI.lb = RD.ci.lb(p0), RD.CI.ub = RD.ci.ub(p0))
  }
  return(out)
}

## Estimating risk difference and CI for specific baseline risks (p0) and visualization
RD_CI <- RD.cond(sid, tid, e, n, data = data,  b.iter = 1000, from = 0.1, to = 0.45,
  ylim = c(0, 0.6), p0 = 0.3)
RD_CI