#' Makes predictions from an eadrm object
#'
#' Similar to other predict methods, this function predicts fitted
#' values (and optionally confidence intervals) for a eadrm object.
#'
#' @param object Fitted model of class eadrm.
#' @param newx Vector of new concentration values at which predictions
#' are to be made. Defaults to the concentration values that were used
#' to fit the model.
#' @param ci.obj Output from eadrm.ci that is used to compute
#' confidence intervals. Defaults to NULL, in which case no confidence
#' intervals are computed.
#' @param ... Additional arguments passed to or from other methods.
#' Currently ignored.
#' @return If no confidence intervals are requested, a vector
#' of predicted responses for each concentration in newx is returned.
#' Otherwise returns a list of three vectors yhat.med, yhat.l95, and
#' yhat.u95, which correspond to the median predicted response and
#' lower/upper 95\% confidence bounds for each concentration in newx.
#' @seealso \code{\link{eadrm}}, \code{\link{eadrm.ci}}
#' @importFrom stats median quantile
#' @export
#' @examples
#' ea.fit <- eadrm(CarboA$y, CarboA$x)
#' predict(ea.fit)
predict.eadrm <- function(object, newx=object$xvals, ci.obj=NULL, ...) {
    if (is.null(ci.obj)) {
        params <- as.numeric(object$params)
        if (object$Model=="h3") {
            out <- (params[1])/(1+exp(params[3]*(log(newx)-log(params[2]))))
        } else if (object$Model=="h4") {
            out <- params[1] - (params[1]-params[2])/
                (1+(newx/params[3])^(params[4]))
        } else if (object$Model=="h5") {
            out <- params[2]+(params[1]-params[2])/
                ((1+exp(params[4]*(log(newx)-log(params[3]))))^params[5])
        } else if (object$Model=="e") {
            out <- params[1]*exp(params[2]*newx)
        }
        return(out)
    }
    else {
        boot.out <- matrix(nrow=length(newx), ncol=ncol(ci.obj$replicate.mat))
        for (i in 1:ncol(ci.obj$replicate.mat)) {
            object.new <- object
            object.new$params <- ci.obj$replicate.mat[,i]
            boot.out[,i] <- predict.eadrm(object.new, newx)
        }
        yhat.med <- apply(boot.out, 1, median)
        yhat.l95 <- apply(boot.out, 1, quantile, probs=0.025, type=8)
        yhat.u95 <- apply(boot.out, 1, quantile, probs=0.975, type=8)
        return(list(yhat.med=yhat.med, yhat.l95=yhat.l95, yhat.u95=yhat.u95))
    }
}

#' Calculates EC50 (or some other specified EC value)
#'
#' Calculates the concentration that induces a response corresponding
#' to a specific proportion between the baseline and maximum. It is most
#' commonly used to compute EC50.
#'
#' @param eadrm.obj Fitted eadrm model object.
#' @param ec.in A value between 0 and 1 corresponding to the desired
#' proportion. Defaults to 0.5, in which case EC50 is computed.
#' @param ci.obj Output from eadrm.ci that is used to compute
#' confidence intervals. Defaults to NULL, in which case no confidence
#' intervals are computed.
#' @return If no confidence intervals are requested, it returns the
#' concentration corresponding to the requested proportion.
#' Otherwise returns a list of three values ec.med, ec.l95, and
#' ec.u95, which correspond to the median concentration and the
#' corresponding lower/upper 95\% confidence bounds.
#' @seealso \code{\link{eadrm}}, \code{\link{eadrm.ci}}
#' @importFrom stats uniroot
#' @export
#' @examples
#' ea.fit <- eadrm(CarboA$y, CarboA$x)
#' calc.ec(ea.fit)
calc.ec <- function(eadrm.obj, ec.in=0.5, ci.obj=NULL) {
    if (is.null(ci.obj)) {
        rootf <- function(x, eadrm.in, resp.in) {
            return(predict.eadrm(eadrm.in, newx=x)-resp.in)
        }
        if (eadrm.obj$Model=="h3") {
            min.y <- 0
            max.y <- as.numeric(eadrm.obj$params)[1]
        } else if (eadrm.obj$Model=="h4" | eadrm.obj$Model=="h5") {
            min.y <- as.numeric(eadrm.obj$params)[2]
            max.y <- as.numeric(eadrm.obj$params)[1]
            if (max.y<min.y) {
                junk <- max.y
                max.y <- min.y
                min.y <- junk
            }
        } else if (eadrm.obj$Model=="e") {
            if (as.numeric(eadrm.obj$params)[2]<0) {
                min.y <- 0
                max.y <- as.numeric(eadrm.obj$params)[1]
            } else {
                min.y <- as.numeric(eadrm.obj$params)[1]
                max.y <- predict.eadrm(eadrm.obj, newx=max(eadrm.obj$xvals))
            }
        }
        y.target <- min.y + ec.in*(max.y - min.y)
        return(uniroot(rootf, c(0,1), extendInt="yes", eadrm.in=eadrm.obj,
                       resp.in=y.target)$root)
    }
    else {
        boot.out <- rep(NA, ncol=ncol(ci.obj$replicate.mat))
        for (i in 1:ncol(ci.obj$replicate.mat)) {
            eadrm.obj.new <- eadrm.obj
            eadrm.obj.new$params <- ci.obj$replicate.mat[,i]
            boot.out[i] <- calc.ec(eadrm.obj.new, ec.in)
        }
        ec.med <- median(boot.out)
        ec.l95 <- quantile(boot.out, probs=0.025, type=8)
        ec.u95 <- quantile(boot.out, probs=0.975, type=8)
        return(list(ec.med=ec.med, ec.l95=ec.l95, ec.u95=ec.u95))
    }
}

#' Finds the dose that corresponds to a particular level of the response
#'
#' Calculates the concentration that induces a particular level of the
#' response.
#'
#' @param eadrm.obj Fitted eadrm model object.
#' @param response The desired response level.
#' @param ci.obj Output from eadrm.ci that is used to compute
#' confidence intervals. Defaults to NULL, in which case no confidence
#' intervals are computed.
#' @return If no confidence intervals are requested, it returns the
#' concentration corresponding to the specified response.
#' Otherwise returns a list of three values ed.med, ed.l95, and
#' ed.u95, which correspond to the median concentration and the
#' corresponding lower/upper 95\% confidence bounds.
#' @seealso \code{\link{eadrm}}, \code{\link{eadrm.ci}}
#' @importFrom stats uniroot
#' @export
#' @examples
#' ea.fit <- eadrm(CarboA$y, CarboA$x)
#' calc.ec(ea.fit)
calc.ed <- function(eadrm.obj, response, ci.obj=NULL) {
    if (is.null(ci.obj)) {
        rootf <- function(x, eadrm.in, resp.in) {
            return(predict.eadrm(eadrm.in, newx=x)-resp.in)
        }
        return(uniroot(rootf, c(0,1), extendInt="yes", eadrm.in=eadrm.obj,
                       resp.in=response)$root)
    }
    else {
        boot.out <- rep(NA, ncol=ncol(ci.obj$replicate.mat))
        for (i in 1:ncol(ci.obj$replicate.mat)) {
            eadrm.obj.new <- eadrm.obj
            eadrm.obj.new$params <- ci.obj$replicate.mat[,i]
            boot.out[i] <- calc.ed(eadrm.obj.new, response)
        }
        ed.med <- median(boot.out)
        ed.l95 <- quantile(boot.out, probs=0.025, type=8)
        ed.u95 <- quantile(boot.out, probs=0.975, type=8)
        return(list(ed.med=ed.med, ed.l95=ed.l95, ed.u95=ed.u95))
    }
}

#' Plot an eadrm object
#'
#' Plots the data used to fit an eadrm object as well as the fitted
#' dose-response curve
#'
#' @param x Object of class eadrm to plot.
#' @param ... Additional arguments to plot. Currently ignored.
#' @seealso \code{\link{eadrm}}
#' @importFrom graphics lines plot points
#' @export
#' @examples
#' ea.fit <- eadrm(CarboA$y, CarboA$x)
#' plot(ea.fit)
plot.eadrm <- function(x, ...) {
    x.min <- max(c(0, min(x$xvals)-0.1))
    x.max <- max(x$xvals)+0.1
    x.minr <- max(c(0, min(x$xvals)))
    x.maxr <- max(x$xvals)
    x1 <- seq(x.min, x.max, by=0.00001)
    x2 <- predict.eadrm(x, newx=x1)
    plot(c(x$xvals, x1[x1>=x.minr&x1<=x.maxr]),
         c(x$yvals, x2[x1>=x.minr&x1<=x.maxr]), xlab="Dose",
         ylab="Response", type="n")
    points(x$xvals, x$yvals)
    lines(x1, x2)
}

#' Computes confidence intervals for an eadrm model fit
#'
#' Calculates confidence intervals for an eadrm model fit by
#' repeatedly fitting the model to the same data set and examining
#' the distribution of the coefficients.
#'
#' @param obs A vector of response values (y-values).
#' @param xvals A vector of doses (x-values).
#' @param model Type of dose-response model to fit. Possible values
#' include "h3", "h4", and "h5" (corresponding to 3-parameter,
#' 4-parameter, and 5-parameter log-logistic models, respectively) and
#' "e" (corresponding to an exponential model). Defaults to "h4".
#' @param B Number of replicate models to fit. Defaults to 1000.
#' @param ... Additional parameters for the eadrm function.
#' @section Details:
#' This function calls the \code{\link{eadrm}} function B times with
#' the same parameters and records the model coefficients for each
#' iteration of the model. Confidence intervals for the coefficients
#' are calculated by examining the quantiles of the distribution
#' of the coefficients over the B iterations. A matrix of the
#' coefficients for each iteration is also calculated. This matrix
#' can be used to compute confidence intervals for predicted values
#' and estimates of EC50.
#' @return A list containing the following elements:
#' \describe{
#' \item{med.est:}{A vector of the median values of the coefficients
#' across the B iterations}
#' \item{l95.est,u95est:}{Vectors of the lower/upper 95\% confidence
#' bounds for the coefficients across the B iterations}
#' \item{replicate.mat:}{A p x B matrix, where p is the number of
#' coefficients in the model. Each column of B corresponds to the
#' coefficients for one fitted model.}
#' }
#' @seealso \code{\link{eadrm}}, \code{\link{predict.eadrm}},
#' \code{\link{calc.ec}}, \code{\link{calc.ed}}
#' @export
#' @examples
#' \donttest{ea.ci <- eadrm.ci(CarboA$y, CarboA$x)}
eadrm.ci <- function(obs, xvals, model='h4', ..., B=1000) {
    first.boot <- as.numeric(eadrm(obs, xvals, model, ...)$params)
    cur.boot <- matrix(nrow=length(first.boot), ncol=B)
    cur.boot[,1] <- first.boot
    for (i in 2:B) {
        cur.boot[,i] <- as.numeric(eadrm(obs, xvals, model, ...)$params)
    }
    if (model=="h3") {
        cur.names <- c("EMAX","EC50","W")
    } else if (model=="h4") {
        cur.names <- c("EMAX","EMIN","EC50","W")
    } else if (model=="h5") {
        cur.names <- c("EMAX","EMIN","EC50","W","f")
    } else if (model=="e") {
        cur.names <- c("B","K")
    }
    med.est <- apply(cur.boot, 1, median)
    l95.est <- apply(cur.boot, 1, quantile, probs=0.025, type=8)
    u95.est <- apply(cur.boot, 1, quantile, probs=0.975, type=8)
    rownames(cur.boot) <- cur.names
    names(med.est) <- cur.names
    names(l95.est) <- cur.names
    names(u95.est) <- cur.names
    return(list(median.est=med.est, l95.ci=l95.est, u95.ci=u95.est,
                replicate.mat=cur.boot))
}
