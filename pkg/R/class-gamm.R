# gamm/gamm4 support

`gamm` <-
function(formula, random = NULL, ..., lme4 = inherits(random, "formula")) {
	pkg <- if(lme4) "gamm4" else "mgcv"
    if (!require(pkg, character.only = TRUE)) stop("'gamm' requires package '", pkg, "' to be installed")
	FUN <- if(lme4) gamm4::gamm4 else mgcv::gamm
	cl <- match.call(FUN)
	cl$lme4 <- lme4
	#cl[[1]] <- call("::", as.name("MuMIn"), as.name("gamm"))
	structure(c(FUN(formula, random, ...), list(call=cl)), class=c(if(lme4) "gamm4", "gamm", "list"))
}

`update.gamm` <-
function(object, ...) {
# or, if call is as attribute: object$call <- attr(object, "call")
	update.default(object, ...)
}

`print.gamm` <-
function(x, ...) {
	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
	cat("--- \n")
	print(x[[if(inherits(x, "gamm4")) "mer" else "lme"]])
	cat("--- \n")
	print(x$gam)
	invisible(x)
}


`logLik.gamm` <-
function (object, ...)
	logLik(object[[if(is.null(object$lme)) "mer" else "lme"]], ...)

`formula.gamm` <-
function (x, ...) formula(x$gam, ...)

# XXX: compatibility with R < 2.13.0
if (exists("nobs", mode = "function", where = "package:stats", inherits = FALSE)) {
	`nobs.gamm` <-
	function (object, ...)  stats:::nobs.glm(object$gam, ...)
} else {
	`nobs.gamm` <-
	function (object, ...) nobs.glm(object$gam, ...)
}

`coeffs.gamm` <-
function (model) coef(model$gam)

`getAllTerms.gamm` <-
function (x, ...) getAllTerms(x$gam, ...)

`tTable.gamm` <-
function (model, ...) tTable(model$gam, ...)