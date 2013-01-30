`DebugPrint` <- function(x) { cat(deparse(substitute(x)), "= \n") ; print(x) }
`srcc` <- function() {
	ret <- eval(expression(source("clipboard", local = TRUE)), .GlobalEnv)
	return(if(ret$visible) ret$value else invisible(ret$value))
}

#if (!exists("getElement", mode = "function", where = "package:base", inherits = FALSE)) {
`getElement` <- function (object, name) {
    if (isS4(object))
		if (.hasSlot(object, name)) slot(object, name) else NULL
    else object[[name, exact = TRUE]]
}
#}

# cbind list of data.frames omitting duplicated column (names)
`cbindDataFrameList` <-
function(x) {
	dfnames <- unlist(lapply(x, colnames), use.names = FALSE)
	uq <- !duplicated(dfnames)
	res <- do.call("cbind", x)[,uq]
	colnames(res) <- dfnames[uq]
	return(res)
}

# same for rbind, check colnames and add NA's when any are missing
`rbindDataFrameList` <-
function(x) {
	all.colnames <- unique(unlist(lapply(x, colnames), use.names = FALSE))
	x <- lapply(x, function(y) {
		y[all.colnames[!(all.colnames %in% colnames(y))]] <- NA
		return(y[all.colnames])
	})
	return(do.call("rbind", x))
}

`videntical` <-
function(x) all(vapply(x[-1L], identical, logical(1L), x[[1L]]))

# Check class for each object in a list
`linherits` <- function(x, whats) {
	as.logical(vapply(x, inherits, integer(length(whats)), names(whats),
		which=TRUE)) == whats
}

# substitute has(a, !b, ...) for !is.na(a) & is.na(b) ..., in expression
`.substHas` <- function(e) {
	if(is.expression(e)) e <- e[[1L]]
	n <- length(e)
	if(n == 1L) return(e)
	if(e[[1L]] != "has") {
		for(i in 1L:n) e[[i]] <- .substHas(e[[i]])
		return(e)
	}
	res <- NULL
	for(i in seq.int(2L, n)) {
		ex <- if(length(e[[i]]) == 2L && e[[i]][[1L]] == "!")
			call("is.na", e[[i]][[2L]]) else
			call("!", call("is.na", e[[i]]))
		res <- if(i == 2L) ex else call("&", res, ex)
	}
	res <- call("(", res)
	return(res)
}

# substitute function calls in 'e'. 'name' is replaced by 'fun.to'.
`.substFun` <- function(e, name, fun.to, ignore.I = TRUE) {
	if(is.expression(e)) e <- e[[1L]]
	n <- length(e)
	if(n == 1L && !is.call(e)) return(e)
	if(ignore.I && e[[1L]] == "I") return(e)
	if(n != 1L) for(i in 2L:n) e[[i]] <- .substFun(e[[i]], name, fun.to, ignore.I = ignore.I)
	if(e[[1L]] == name) e[[1L]] <- as.name(fun.to)
	return(e)
}

# substitute function calls in 'e'. 'func' must take care of the substitution job.
`.substFun4Fun` <- function(e, name, func = identity, ...) {
	if(is.expression(e)) e <- e[[1L]]
	n <- length(e)
	if(n == 0L) return(e) else if (n == 1L) {
		if (!is.call(e)) return(e)
	} else for(i in 2L:n) e[i] <- list(.substFun4Fun(e[[i]], name, func, ...))
	if(e[[1L]] == name) e <- func(e, ...)
	return(e)
}

# evaluate 'expr' in 'env' after adding variables passed as ...
.evalExprIn <- function(expr, env, enclos, ...) {
	list2env(list(...), env)
	eval(expr, envir = env, enclos = enclos)
}

# substitute names for varName[1], varName[2], ... in expression
`.subst4Vec` <- function(expr, names, varName, n = length(names), fun = "[") {
	eval(call("substitute", expr,
		env = structure(lapply(seq_len(n), function(i) call(fun, varName, i)), names = names)),
		envir = NULL)
}

# tries to make a list of element names
`.makeListNames` <- function(x) {
	nm <- names(x)
	lapply(seq_along(x), function(i) {
		if(is.null(nm) || nm[i] == "") {
			switch(mode(x[[i]]),
				call = deparse(x[[i]], control = NULL),
				symbol =, name = as.character(x[[i]]),
				NULL =, logical =, numeric =, complex =, character = x[[i]], i
				)
		} else nm[i]
	})
}


# test if dependency chain is satisfied: x[n] can be TRUE only if x[1:n] are also TRUE
`.subset_dc` <- function(...) {
	n <- length(x <- c(...))
	if(n > 1L) all(x[-n] >= x[-1L]) else TRUE
}

# vectorized version of .subset_do (used within subset.model.selection)
`.subset_vdc` <- function(...) apply(cbind(..., deparse.level = 0L), 1L, .subset_dc)


`prettyEnumStr` <- function(x, sep = ", ", sep.last = gettext(" and "), quote = TRUE) {
	n <- length(x)
	if(is.function(quote))
		x <- quote(x) else {
			if(identical(quote, TRUE)) quote <- '"'
			if(is.character(quote)) x <- paste(quote, x, quote, sep = "")
		}
	paste(x, if(n > 1L) c(rep(sep, n - 2L), sep.last, "") else NULL,
		collapse = "", sep = "")
}

# `splitList` <- function (x, k) {
    # n <- length(x)
    # ret <- unname(split.default(x, findInterval(seq_len(n), seq(0L, n +
        # 1L, length = k + 1L))))
	# if(k > n) ret <- c(ret, vector(k - n, mode = "list"))
	# ret
# }


`.parallelPkgCheck` <- function(quiet = FALSE) {
	# all this is to trick the R-check
	if(!("snow" %in% loadedNamespaces())) {
		if(getRversion() < "2.14.0") {
			if(length(.find.package("snow", quiet = TRUE)))
				do.call("require", list("snow"))
		} else if(length(.find.package("parallel", quiet = TRUE)))
			do.call("require", list("parallel", quiet = TRUE))
	}
	if(!exists("clusterCall", mode = "function")) {
		if(quiet) return(FALSE) else
			stop("cannot find function 'clusterCall'")
	} else return(TRUE)
}

`clusterVExport` <- local({
   `getv` <- function(obj)
		for (i in names(obj)) assign(i, obj[[i]], envir = as.environment(1L))
	function(cluster, ...) {
		Call <- match.call()
		Call$cluster <- NULL
		Call <- Call[-1L]
		vars <- list(...)
		vnames <- names(vars)
		#if(!all(sapply(Call, is.name))) warning("at least some elements do not have syntactic name")
		if(is.null(vnames)) {
			names(vars) <- vapply(Call, deparse, character(1L), control = NULL,
				nlines = 1L)
		} else if (any(vnames == "")) {
			names(vars) <- ifelse(vnames == "", vapply(Call, deparse,
				character(1L), control = NULL, nlines = 1L), vnames)
		}
		get("clusterCall")(cluster, getv, vars)
		# clusterCall(cluster, getv, vars)
	}
})

# test if 'x' can be updated (in current environment or on a cluster)
# level is 0/FALSE - no checking, 1 - check if variables and functions exist,
# >1 - reevaluate x and compare with original 
`testUpdatedObj` <- function(cluster = NA, x, call = .getCall(x),
	level = 1L, exclude = "subset") {
	
	if(isTRUE(level)) level <- 2L

	if (level > 0L) {
		xname <- deparse(substitute(x))
		doParallel <- inherits(cluster, "cluster")
		if(doParallel) {
			clusterCall <- get("clusterCall")
			whereStr <- gettext(" in the cluster nodes' environment")
			csapply <- function(...) clusterCall(cluster, "sapply", ...)
		} else {
			whereStr <- ""
			csapply <- function(...) sapply(...)
		}
		if(is.null(call)) stop(gettextf("'%s' has no call component", xname))
		call.orig <- call
		if(!is.null(call$data)) {
			# get rid of formulas, as they are evaluated within 'data'
			call <- call[!sapply(call, function(x) "~" %in% all.names(x))]
			call[exclude] <- NULL
		}	
		
		v <- all.vars(call, functions = FALSE)
		if(!all(z <- unlist(csapply(v, "exists", where = 1L)))) {
			z <- unique(names(z[!z]))
			stop(sprintf(ngettext(length(z), "variable %s not found%s",
				"variables %s not found%s"), prettyEnumStr(z, quote = "'"), whereStr))
			}
		vfun <- all.vars(call, functions = TRUE)
		if(!all(z <- unlist(csapply(vfun[!(vfun %in% v)], "exists",
			mode = "function", where = 1L)))) {
			zz <- unique(names(z[!z]))
			stop(sprintf(ngettext(length(zz), "function %s not found%s",
				"functions %s not found%s"), prettyEnumStr(zz, quote = "'"), whereStr))
			}
		if(level > 1L && !missing(x)) {
			if(doParallel) {
				# XXX: Import: clusterCall
				if(!all(vapply(lapply(clusterCall(cluster, eval, call.orig), all.equal, x), isTRUE, TRUE)))
					stop(gettextf("'%s' evaluated on the cluster nodes differs from the original one",
				xname))
			} else if (!isTRUE(all.equal(x, update(x))))
				stop(gettextf("updated '%s' differ(s) from the original one", xname))
		}
	}
}

`tryCatchWE` <- function (expr) {
	Warnings <- NULL
	list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
		warning = function(w) {
			Warnings <<- c(Warnings, list(w))
			invokeRestart("muffleWarning")
		}), warnings = Warnings)
}

# like apply(, 2) but returns a list (does not do any checking)
`applyrns` <- function (X, FUN, ...) {
	n <- nrow(X)
	ret <- vector(n, mode = "list")
	for(i in seq_len(n)) if(!is.null(z <- FUN(X[i, ], ...))) ret[[i]] <- z
	ret
}