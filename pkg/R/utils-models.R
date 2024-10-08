fixLogLik <-
function(ll, object) {
	if(is.null(attr(ll, "nall")) && is.null(attr(ll, "nobs")))
		attr(ll, "nobs") <- nobs(object)
	ll
}

`.getLik` <-
function(x) {
    if(isGEE(x)) {
		list(logLik = quasiLik, name = "qLik")
	} else {
		list(logLik = logLik, name = "logLik")
    }
}


`.getRank` <-
function(rank = NULL, rank.args = NULL, object = NULL, ...) {
	rank.args <- c(rank.args, list(...))

	if(is.null(rank)) {
		x <- NULL # just not to annoy R check
		IC <- as.function(c(alist(x =, do.call("AICc", list(x)))))
		attr(IC, "call") <- call("AICc", as.name("x"))
		class(IC) <- c("function", "rankFunction")
		return(IC)
	} else if(inherits(rank, "rankFunction") && length(rank.args) == 0L) {
		return(rank)
	} else if(is.list(rank) && length(rank) == 1L && is.function(rank[[1L]])) {
        srank <- names(rank)[1L]
        rank <- rank[[1L]]
    } else {
        srank <- substitute(rank, parent.frame())
        if(srank == "rank") srank <- substitute(rank)
    }

	rank <- match.fun(rank)
	ICName <- switch(mode(srank), call = as.name("IC"), character = as.name(srank), name=, srank)
	ICarg <- c(list(as.name("x")), rank.args)
	ICCall <- as.call(c(ICName, ICarg))
	IC <- as.function(c(alist(x =), list(substitute(do.call("rank", ICarg),
		list(ICarg = ICarg)))))

	if(!is.null(object)) {
		test <- IC(object)
		if (!is.numeric(test) || length(test) != 1L)
			stop("'rank' should return numeric vector of length 1")
	}

	attr(IC, "call") <- ICCall
	class(IC) <- c("function", "rankFunction")
	IC
}


# Like `regmatches`, but for captured groups, and simplistic.
# Only does results of `regexpr`.
.matches <-
function (x, m) { 
	cs <- attr(m, "capture.start")
	cl <- attr(m, "capture.length")
	rval <- array(NA_character_, dim = dim(cs))
	for(i in seq_along(x))
		rval[i, ] <- substring(x[i], cs[i, ], cs[i, ] + cl[i, ] - 1L)
	rval
}

# sorts alphabetically interaction components in model term names
# if 'peel', tries to remove coefficients wrapped into function-like syntax
# (this is meant mainly for 'unmarkedFit' models with names such as "psi(a:b:c)")
# This unwrapping is done only if ALL names are suitable for it.
# FIXME: this function would also "fix" strings like "log(b:a)" (but not 
# "log(`b:a`)") - but this is unlikely to be a valid model term.
`fixCoefNames` <-
function(x, peel = TRUE) {
	if(!length(x)) return(x)
	ia <- grep(":", x, fixed = TRUE)
	if(!length(ia)) return(structure(x, order = rep.int(1L, length(x))))
	
	ixi <- x[ia]
	if(peel) {
		# peel only when ALL items are prefixed/wrapped
		# for pscl::hurdle. Cf are prefixed with count_|zero_
		if(peel <- all(startsWith(ixi, c("count_", "zero_")))) {
			pos <- regexpr("_", ixi, fixed = TRUE)
			peelpfx <- substring(ixi, 1L, pos)
			peelsfx <- ""
			ixi <- substring(ixi, pos + 1L)
		} else {
			# unmarkedFit with its phi(...), lambda(...) etc...
			if(peel <- all(grepl("^[a-zA-Z]{2,5}\\(.+\\)$", x, perl = TRUE))) {
				# only if 'XXX(...)', i.e. exclude 'XXX():YYY()' or such
                # assumes coefficient types are ascii letters only 
                # and of length >= 2
				m <- regexpr("^(([a-zA-Z]{2,5})\\(((?:[^()]*|(?1))*)\\))$", ixi, perl = TRUE, useBytes = TRUE)
                cptgrps <- .matches(ixi, m)
				if(peel <- all(nzchar(cptgrps[, 2L]))) {
					peelpfx <- paste0(cptgrps[, 2L], "(")
					peelsfx <- ")"
					ixi <- cptgrps[, 3L]
				}
			}
		}
	}
	# replace {...}, [...], (...), ::, and ::: with placeholders
	m <- gregexpr("(?:\\{(?:[^\\{\\}]*|(?0))*\\}|\\[(?:[^\\[\\]]*|(?0))*\\]|\\((?:[^()]*|(?0))*\\)|(?>:::?))", ixi, perl = TRUE)
	xtpl <- ixi
	regmatches(xtpl, m) <- lapply(m, function(x) {
		if((ml <- attr(x, "match.length"))[1L] == -1L) return(character(0L))
		sapply(ml, function(n) paste0(rep("_", n), collapse = ""))
	})
	
	# split by ':' and sort
	splits <- gregexpr(":", xtpl, fixed = TRUE)
	ixi <- mapply(function(x, p) {
		if(p[1L] == -1) return(x)
		paste0(base::sort(substring(x, c(1L, p + 1L), c(p - 1L, nchar(x)))), collapse = ":")
	}, ixi, splits, USE.NAMES = FALSE, SIMPLIFY = TRUE)
	
	if(peel) ixi <- paste0(peelpfx, ixi, peelsfx)

	x[ia] <- ixi
	ord <- rep.int(1L, length(x))
	ord[ia] <- vapply(splits, length, 0L) + 1L
	attr(x, "order") <- ord
	x
}


## like 'strsplit', but ignores split characters within quotes and matched
## parentheses
expr.split <-
function(x, split = ":",
	paren.open = c("(", "[", "{"), paren.close = c(")", "]", "}"),
	quotes = c("\"", "'", "`"), esc = "\\",
	prepare = NULL) {

	## error checking:
	#if(length(paren.open) != length(paren.close))
	#	stop("'paren.open' and 'paren.close' are not of the same length")
	#if(any(test <- vapply(c('paren.open', 'paren.close', 'quotes', 'esc', 'split'), function(x, frame) {
	#	any(nchar(get(x, frame, inherits = FALSE)) != 1L)
	#}, FALSE, frame = sys.frame()))) {
	#	stop(sprintf(ngettext(sum(test), "argument %s is not single character",
	#		"arguments %s are not single character") , prettyEnumStr(names(test)[test])),
	#		 domain = "R-MuMIn")
	#}

	x0 <- x
	if(is.function(prepare)) x <- prepare(x)
	m <- length(x)
	n <- nchar(x)
	res <- vector("list", m)
	for(k in 1L:m) {
		pos <- integer(0L)
		inquote <- ch <- ""
		inparen <- integer(3L)
		for(i in seq.int(n[k])) {
			chprv <- ch
			ch <- substr(x[k], i, i)
			if(nzchar(inquote)) { # in quotes
				if(chprv == esc && ch == esc) ch <- " " else
					if(chprv != esc && ch == inquote)	inquote <- ""
			} else {
				inparen[j] <- inparen[j <- (inparen != 0L) & (ch == paren.close)] - 1L
				if(ch %in% quotes)
					inquote <- ch else if (any(j <- (ch == paren.open)))
					inparen[j] <- inparen[j] + 1L else if (all(inparen == 0L) && ch == split)
					pos <- c(pos, i)
			}
		}
		res[[k]] <- substring(x0[k], c(1L, pos + 1L), c(pos - 1L, n[k]))
	}
	res
}

getResponseFormula <-
function(f) {
	f <- if(!is.null(tf <- attr(f, "terms"))) {
		formula(tf)
	} else formula(f)
	if((length(f) == 2L) || (is.call(f[[2L]]) && f[[2L]][[1L]] == "~"))
		0 else f[[2L]]
}


#Tries to find out whether the models are fitted to the same data
checkIsModelDataIdentical <-
function(models, error = TRUE) {

	cl <- sys.call(sys.parent())
	err <-  if (error) 	function(x) stop(simpleError(x, cl))
		else function(x) warning(simpleWarning(x, cl))
	res <- TRUE

	responses <- lapply(models, function(x) getResponseFormula(formula(x)))

 	if(!all(vapply(responses[-1L], "==", FALSE, responses[[1L]]))) {
		err("response differs between models")
		res <- FALSE
	}

	# XXX: need to compare deparse'd 'datas' due to ..1 bug(?) in which dotted
	#  arguments (..1 etc) passed by lapply are not "identical"
	datas <- vapply(lapply(models, function(x) get_call(x)$data), asChar, "")

	# XXX: when using only 'nobs' - seems to be evaluated first outside of MuMIn
	# namespace which e.g. gives an error in glmmML - the glmmML::nobs method
	# is faulty.
	nresid <- vapply(models, function(x) nobs(x), 1) # , nall=TRUE

	if(!all(sapply(datas[-1L], identical, datas[[1L]])) ||
		!all(nresid == nresid[[1L]])) { # better than 'nresid[-1L] == nresid[[1L]]'
		# XXX: na.action checking here
		err("models are not all fitted to the same data")
		res <- FALSE
	}
	invisible(res)
}

.checkNaAction <-
function(x, cl = get_call(x),
		 naomi = c("na.omit", "na.exclude"), what = "model",
		 envir = parent.frame()) {
	naact <- NA_character_
	msg <- NA_character_

	# handles strings, symbols and calls (let's naively assume no one tries to pass
	# anything else here)
	.getNAActionString <- function(x) {
		if(is.symbol(x)) {
			x <- as.character(x)
		} else if(is.call(x)) {
			x <- eval(x, envir)
			if(is.symbol(x)) x <- as.character(x)
		}
		return(x)
	}
	# TEST:
	#.checkNaAction(list(call = as.call(alist(fun, na.action = getOption("na.action", default = na.fail)))))
	#.checkNaAction(list(call = as.call(alist(fun, na.action = na.fail))))
	#.checkNaAction(list(call = as.call(alist(fun, na.action = na.omit))))
	
	if (!is.null(cl$na.action)) {
		naact <- .getNAActionString(cl$na.action)
		if (naact %in% naomi)
			msg <- sprintf("%s uses 'na.action' = \"%s\"", what, naact)
	} else {
		naact <- formals(eval(cl[[1L]], envir))$na.action
		if (missing(naact)) {
			naact <- getOption("na.action")
			if(is.function(naact)) {
				statsNs <- getNamespace("stats")
				for(i in naomi) if(identical(get(i, envir = statsNs, inherits = FALSE), naact,
					ignore.environment = TRUE)) {
					naact <- i
					break
					}
			}

			naact <- .getNAActionString(naact)
			if (is.character(naact) && (naact %in% naomi))
				msg <- sprintf("%s's 'na.action' argument is not set and options('na.action') is \"%s\"",
					what, naact)
		} else if (!is.null(naact)) {
			naact <- .getNAActionString(naact)
			if (naact %in% naomi)
				msg <- sprintf("%s uses the default 'na.action' = \"%s\"", what, naact)
		}
	}
	res <- is.na(msg)
	attr(res, "na.action") <- naact
	attr(res, "message") <- msg
	res
}

`abbreviateTerms` <-
function(x, minlength = 4, minwordlen = 1,
	capwords = FALSE, deflate = FALSE) {
	if(!length(x)) return(x)
	if(deflate) dx <-
		#gsub("([\\(,]) *\\w+ *= *(~ *(1 *[\\+\\|]?)?)? *", "\\1", x, perl = TRUE)
		gsub("([,\\(\\[]|^)( *~ *)(1 *([\\|\\+] *)?)?", "\\1",
			gsub("([\\(,]) *\\w+ *= *", "\\1", x, perl = TRUE), perl = TRUE)
		else dx <- x

	#.DebugPrint(x)
	s <- strsplit(dx, "(?=[\\W_])", perl = TRUE)
	# remove I(...):
	s <- lapply(s, function(z) {
		z <- if((n <- length(z)) > 3L && all(z[c(1L, 2L, n)] == c("I", "(", ")")))
			z[3L:(n - 1L)] else z
		z[z != " "]
	})
	v <- unique(unlist(s, use.names = FALSE))
	i <- grep("[[:alpha:]]", v, perl = FALSE)
	av <- v
	if(length(i)) {
		tb <- rbindDataFrameList(lapply(s, function(x)
			as.data.frame(rbind(c(table(x))), stringsAsFactors = TRUE)))
		tb[is.na(tb)] <- 0L

		if(length(v) > length(i)) minlength <-
			minlength - max(c(0L, apply(tb[, v[-i], drop = FALSE], 1L,
			"*", nchar(colnames(tb[, v[-i], drop = FALSE])))))
		n <- min(minlength / rowSums(tb[, v[i], drop = FALSE]))
		if(deflate) {
			repl1 <- c("TRUE" = "T", "FALSE" = "F", "NULL" = "")
			for(j in seq_along(repl1)) av[av == names(repl1)[j]] <- repl1[j]
		}
		av[i] <- abbreviate(av[i], max(n, minwordlen))
		if(capwords) av[i] <- paste0(toupper(substring(av[i], 1L, 1L)),
				tolower(substring(av[i], 2L)))
	}
	for(j in seq_along(s)) s[[j]] <- paste(av[match(s[[j]], v)], collapse = "")
	names(av) <- v
	structure(unlist(s), names = x, variables = av[i])
}

`modelDescr` <-
function(models, withModel = FALSE, withFamily = TRUE,
	withArguments = TRUE, remove.cols = c("formula", "random", "fixed", "model",
	"data", "family", "cluster", "model.parameters"),
	remove.pattern = "formula$", ...) {

	if(withModel) {
		allTermsList <- lapply(models, function(x) {
			tt <- getAllTerms(x)
			rtt <- attr(tt, "random.terms")
			if(!is.null(rtt))
				rtt[i] <- paste0("(", rtt[i <- !grepl("^[a-z]*\\(.+\\)$", rtt)], ")")
			c(tt, rtt)
		})
		allTerms <- unique(unlist(allTermsList))
		abvtt <- abbreviateTerms(allTerms, ...)
		variables <- attr(abvtt, "variables")
		abvtt <- gsub("\\(1 \\| (\\S+)(?: %in%.*)?\\)", "(\\1)", abvtt, perl = TRUE)
		abvtt <- sapply(allTermsList, function(x) paste(abvtt[match(x, allTerms)],
			collapse = "+"))
	} else abvtt <- variables <- NULL

	if(withFamily) {
		fam <- sapply(models, function(x) tryCatch(unlist(family(x)[c("family",
			"link")]), error = function(e) character(2L)))

		f <- fam[1L, ]
		f[is.na(f)] <- ""
		#f <- vapply(strsplit(f, "(", fixed = TRUE), "[", "", 1L)
		#f[f == "Negative Binomial"] <- "negative.binomial"
		#fam <- cbind(fam, unlist(MASS::negative.binomial(1.345)[c("family", "link")]))
	    f <- sub("(?:\\((.*)\\))?$", "(\\1", f)
		f <- paste0(f, ifelse(substring(f, nchar(f)) == "(", "", ","), fam[2, ], ")")
		fam <- f		
		
		#fam[2L, fam[2L, ] ==
			#vapply(unique(f),
				#function(x) {
					#rval <- if(is.na(x)) NA_character_ else formals(get(x))$link[1L]
					#if(!is.character(rval)) NA_character_ else rval
					#}, FUN.VALUE = "")[f]] <- NA_character_
		#j <- !is.na(fam[2L,])
		#fnm <- fam[1L, j]
		#fnm <- ifelse(substring(fnm, nchar(fnm)) != ")",
		#paste0(fnm, "("), paste0(substring(fnm, 1, nchar(fnm) - 1),
				#", "))
		#fam[1L, j] <- paste0(fnm, fam[2L, j], ")")
	}

	if(withArguments) {
		
		cl <- lapply(models, get_call)
		haveNoCall <-  vapply(cl, is.null, FALSE)
		cl[haveNoCall] <- lapply(cl[haveNoCall], function(x) call("none", formula = NA))
 		arg <- lapply(cl, function(x) sapply(x[-1L], function(argval)
			switch(mode(argval), character = , logical = argval,
			numeric = signif(argval, 3L), asChar(argval))))
		arg <- rbindDataFrameList(lapply(lapply(arg, t), as.data.frame))
		
		i <- !(colnames(arg) %in% remove.cols)
		i[i][grep(remove.pattern[1L], colnames(arg)[i])] <- FALSE

		arg <- cbind(class = as.factor(sapply(lapply(models, class), "[", 1L)),
			arg[, i, drop = FALSE])
		reml <-	rep(NA, length(models))
		if(!is.null(arg$method)) {
			reml <- ((arg$class == "lme" &
				is.na(arg$method)) | arg$method == "REML")
			arg$method  <- NULL
		}
		if(!is.null(arg$REML)) reml <- ifelse(is.na(arg$REML), reml, arg$REML == "TRUE")
		arg$REML <- as.factor(reml)

		arg <- as.matrix(arg)
		arg[is.na(arg) | arg == "NULL"] <- ""
		arg <- arg[, apply(arg, 2L, function(x) length(unique(x))) != 1L, drop = FALSE]
		if(ncol(arg)) arg <- gsub("([\"'\\s]+|\\w+ *=)","", arg, perl = TRUE)
	}
	ret <- as.data.frame(cbind(model = abvtt, family = if(withFamily) fam else NULL,
		arg, deparse.level = 0L), stringsAsFactors = TRUE)
	attr(ret, "variables") <- variables
	ret
}


# TODO: add theta
# TODO: glmmTMB family
# FIXME: family.betareg is not binomial
#negbin(theta = .1, link = "log")$getTheta()
#ocat(theta=1,link="identity",R=NULL)$getTheta(1)
#ziP(theta = 1, link = "identity",b=0)

family2char <-
function(x, fam = x$family, link = x$link) {
	if(startsWith(fam, "Negative Binomial")) {
		theta <- as.numeric(strsplit(fam, "[\\(\\)]")[[1L]][2L])
		paste0("negative.binomial", "(", theta, ",", link, ")")
	} else if (startsWith(fam, "Tweedie(")) {
		paste0(substr(fam, 1L, nchar(fam) - 1L), ",link=", link, ")")
	} else {
		switch(fam, "scaled t" = "scat",
			   "Beta regression" = if("putTheta" %in% names(x))
					"betar" else "beta.regression",
			   "zero inflated Poisson" = "ziP",
			   "Cox PH" = "cox.ph",
			   "Ordered Categorical" = "ocat",
			   "Multivariate normal" = "mvn",
			   beta = "beta_family",
			   )
		paste0(fam, "(", link, ")")
	}
}


`commonCallStr` <-
function(models, calls = lapply(models, get_call)) {
	
	x <- lapply(calls, as.list)
	alln <- unique(unlist(lapply(x, names)))
	uniq <- vector("list", length(alln))
	names(uniq) <- alln
	uniq[[1L]] <- lapply(x, "[[", 1L)
	for(i in alln[-1]) uniq[[i]] <- lapply(x, "[[", i)
	uniq <- rapply(uniq, classes = "formula", function(x) {
		environment(x)  <- .GlobalEnv
		x
	}, how = "replace")
	uniq <- lapply(uniq, unique)
	nu <- sapply(uniq, length)
	strvarious <- "<*>"
	
	rval <- lapply(uniq, '[[', 1L)
	j <- sapply(rval, inherits, "formula") & nu > 1L
	for(i in which(j)) {
		response <- getResponseFormula(rval[[i]])
		rval[[i]] <- if(identical(response, 0))
			call("~", as.name(sprintf("__%d-rhsform__", nu[i]))) else
			call("~", getResponseFormula(rval[[i]]), as.name(sprintf("__%d-rhsform__", nu[i])))
	}
	
	notj <- !j & nu > 1
	rval[notj] <- paste0("<", nu[notj], " unique values>")
	if(nu[1L] > 1) rval[[1L]] <- paste(sapply(uniq[[1L]], asChar), collapse = "|")
		
	rval <- paste(deparse(rval[[1L]], control = NULL),
		"(", paste(names(rval[-1L]), "=", rval[-1L], collapse = ", "), ")", sep = "")
	
	rval <- gsub("`__(\\d+)-rhsform__`", "<\\1 unique rhs>", rval, perl = TRUE)
	rval
}

updateDeps <-
function(expr, deps) {
	ret <- list()
	env <- sys.frame(sys.nframe())
	expr <- exprapply0(expr, "dc", function(z) {
		v <- vapply(as.list(z[-1L]), asChar, "")
		n <- length(v)
		k <- match(v, colnames(deps))
		for(i in 2L:n) deps[k[1L:(i - 1L)], k[i]] <- TRUE
		assign("deps", deps, envir = env, inherits = FALSE)
		TRUE
	})
	list(deps = deps, expr = expr)
}

# tests if smooth terms for variables in gam/gamm models have the same 'k'
testSmoothKConsistency <-
function(models) {

	# method 2: guess from coefficient names:
	if(inherits(models, "model.selection")) {

		x <- lapply(attr(models, "coefTables"), rownames)
		# XXX: add label 'te(x1,x2)'

		res <- unlist(unname(lapply(x, function(x) {
			s <- grep("^(s|t[ei2])\\(.+\\)\\.\\d+$", x, perl = TRUE)
			if(length(s) != 0L) {
				m <- regexpr("^(?:s|t[ei2])\\((.+)\\)\\.\\d+$",  x[s], perl = TRUE)
				cst <- attr(m, "capture.start")[, 1L]
				y <- substring(x[s], cst, cst + attr(m, "capture.length")[, 1L] - 1L)
				tapply(y, y, length)
			} else NULL
		})), recursive = FALSE)
		names(res) <- sapply(lapply(expr.split(names(res), ","), sort),
			paste0, collapse = ",")
	} else {
		# use information stored in gam objects:
		.getSmoothK <-
		function(x) {
			if(inherits(x, "gamm") || (is.list(x) &&
				(length(x) >= 2L) && identical(names(x)[2L], "gam"))) {
				x <- x[[2L]]
			} else if(!inherits(x, "gam")) return(NULL)
			n <- length(x$smooth)
			rval <- vector("list", n)
			nmv <- character(n)
			for(i in seq_len(n)) {
				y <- x$smooth[[i]]
				if(is.null(y$margin)) {
					nmv[i] <- y$term
					rval[[i]] <- y$bs.dim
				} else {
					nm1 <- vapply(y$margin, `[[`, "", "term")
					o <- order(nm1)
					nmv[i] <- paste0(nm1[o], collapse = ",")
					rval[[i]] <- sapply(y$margin, `[[`, "bs.dim")[o]
				}
				nmv[i] <- paste0(sub("\\(.*", "", y$label), "(", nmv[i], ")")
			}
			names(rval) <- nmv
			rval
		}
		res <- unlist(unname(lapply(models, .getSmoothK)), recursive = FALSE)
	}

	if(!is.null(res)) {
		res <- vapply(split(res, names(res)), function(x) {
			k1 <- x[[1L]]
			for(i in 1L:length(x)) if(!identical(x[[i]], k1)) return(TRUE)
			return(FALSE)
		}, FALSE)


		if(any(res))
			warning("smooth term dimensions differ between models for variables ",
				prettyEnumStr(names(res)[res], quote = "'"),
				". Related coefficients are incomparable."
			)
	}
	invisible()
}
