## TODO: chunk size for evaluate = FALSE

`pdredge` <-
function(global.model, cluster = NULL,
	beta = c("none", "sd", "partial.sd"),
	evaluate = TRUE,
	rank = "AICc", fixed = NULL, m.lim = NULL, m.min, m.max, subset,
	trace = FALSE, varying, extra, ct.args = NULL, 
    deps = attr(allTerms0, "deps"),
    check = FALSE, ...) {

    .Deprecated("dredge")
	
    cl <- match.call()
    cl[[1L]] <- as.symbol("dredge")
    return(eval.parent(cl))
}

expr <- body(.dredge.parallel)
vdef <- list()
MuMIn:::exprApply(expr, "<-", function(e) {
	if(is.symbol(e[[2]])) 
		vdef[[length(vdef) + 1L]] <<- e[[2]]
	e
}) -> .
vdef <- sapply(vdef, deparse)
vall <- all.vars(expr, functions = FALSE, unique = TRUE)
vundef <- setdiff(vall, vdef)



fdef <- list()
MuMIn::exprApply(expr, "function", function(expr1) {
	#print(names(expr1[[2]]))
	vdef1 <- list()
	MuMIn::exprApply(expr1[[3]], "<-", function(e) {
		print(e)
		if(is.symbol(e[[2]])) vdef1[[length(vdef1) + 1L]] <<- e[[2]]
	})
	vdef1 <- sapply(vdef1, deparse)
	formals <- names(expr1[[2]])
	vall <- all.vars(expr1, functions = FALSE, unique = TRUE)
	vundef <- setdiff(vall, c(vdef1, formals))
	fdef[[length(fdef) + 1L]]  <<- list(formals = formals,
		vdef = vdef1, vall = vall, vundef = vundef) 
	expr1
}) -> .


lapply(fdef, "[[", "vundef")
	
	
}) -> .

MuMIn:::.exprapply

`.dredge.parallel` <-
function() {
	
  #cat(deparse(paste(  scan("clipboard", ""))))
	  #, collapse = " <- "))
	  

#    # Define
#	IC <- ICName <- allTerms <- allTerms0 <- applyExtras <- argsOptions <- betaMode <-
#	check <- cluster <- comb.seq <- comb.sfx <- ct.args <- deps <- doParallel <-
#	evaluate <- extraNames <- extraResult <- fvarying <- global.model <-
#	gmCall <- gmCoefNames <- gmEnv <- gmNobs <- hasSubset <- interceptLabel <-
#	lik <- m.max <- m.min <- nExtra <- nIntercepts <- nVariants <- nVars <-
#	nVarying <- ncomb <- resultChunkSize <- rvNcol <- ssEnv <- strbeta <-
#	subsetExpr <- termsOrder <- variants <- variantsFlat <- varyingNames <-
#	vlen <- NULL
	  
	#  for(a in c("IC", "ICName", "allTerms", "allTerms0", "applyExtras", "argsOptions", "betaMode",
	#  "check", "cluster", "comb.seq", "comb.sfx", "ct.args", "deps",
	#  "doParallel", "evaluate", "extraNames", "extraResult", "fvarying",
	#  "global.model", "gmCall", "gmCoefNames", "gmEnv", "gmNobs", "hasSubset",
	#  "interceptLabel", "lik", "m.max", "m.min", "nExtra", "nIntercepts",
	#  "nVariants", "nVars", "nVarying", "ncomb", "resultChunkSize", "rvNcol",
	#  "ssEnv", "strbeta", "subsetExpr", "termsOrder", "variants",
	#  "variantsFlat", "varyingNames", "vlen" ))
	#	assign(a, get(a, parent.frame()))
	# 

###PAR
	qlen <- 25L
	# Imports: clusterCall, clusterApply
	#doParallel <- isTRUE(evaluate) && inherits(cluster, "cluster")
	if(doParallel) {
		.parallelPkgCheck() # XXX: workaround to avoid importing from 'parallel'
		clusterCall <- get("clusterCall")
		clusterApply <- get("clusterApply")
		clusterCall(cluster, "require", .packageName, character.only = TRUE)
		.getRow <- function(X) clusterApply(cluster, X, fun = ".pdredge_process_model")
	} else {
		.getRow <- function(X) lapply(X, pdredge_process_model, envir = props)
		clusterCall <- function(...) NULL
		message("Not using cluster.")
	}
###PAR

###PAR
	# parallel: check whether the models would be identical:
	if(doParallel && check) testUpdatedObj(cluster, global.model, gmCall, level = check)
###PAR


	# BEGIN parallel
	qi <- 0L
	queued <- vector(qlen, mode = "list")
	props <- list(
				gmEnv = gmEnv,
				IC = IC,
				# beta = beta,
				# allTerms = allTerms,
				nExtra = nExtra,
				matchCoefCall = as.call(c(list(
					as.name("matchCoef"), as.name("fit1"),
					all.terms = allTerms, beta = betaMode,
					allCoef = TRUE), ct.args))
				# matchCoefCall = as.call(c(alist(matchCoef, fit1, all.terms = Z$allTerms,
				#   beta = Z$beta, allCoef = TRUE), ct.args))
		)
	if(nExtra) {
		props$applyExtras <- applyExtras
		props$extraResultNames <- names(extraResult)
	}
	props <- as.environment(props)

	if(doParallel) {
		message("DEBUG")
		clusterVExport(cluster,   pdredge_props = props,
								  .pdredge_process_model = pdredge_process_model
								  )
		clusterCall(cluster, eval, call("options", options("na.action")), env = 0L)
	}
	# END parallel

	retColIdx <- if(nVarying) -nVars - seq_len(nVarying) else TRUE

	dotrace <- if(trace == 1L) {
		dotrace <- function()  {
			cat(iComb, ": "); print(clVariant)
			utils::flush.console()
		}
	} else if(trace > 1L) {
		progressBar <- .progbar(max = ncomb, title = "\"dredge\" working...")
		on.exit(.closeprogbar(progressBar))
		function() progressBar(value = iComb,
			title = sprintf("dredge: %d of ca. %.0f subsets", k, (k / iComb) * ncomb))
	} else function() {}


	warningList <- list()

	iComb <- -1L
	while((iComb <- iComb + 1L) < ncomb) {
		varComb <- iComb %% nVariants
		jComb <- (iComb - varComb) / nVariants
		#print(c(iComb, jComb, ncomb, varComb + 1L))
		if(varComb == 0L) {
			isok <- TRUE

			comb <- c(as.logical(intToBits(jComb)[comb.seq]), comb.sfx)
			nvar <- sum(comb) - nIntercepts

			# !!! POSITIVE condition for 'pdredge', NEGATIVE for 'dredge':
			if((nvar >= m.min && nvar <= m.max) &&
				formula_margin_check(comb, deps) &&
				switch(hasSubset,
					# 1 - no subset, 2 - matrix, 3 - expression
					TRUE,                                    # 1
					all(subset[comb, comb], na.rm = TRUE),   # 2
					evalExprInEnv(subsetExpr, env = ssEnv, enclos = parent.frame(),
						comb = comb, `*nvar*` = nvar),		 # 3
					TRUE
					)
				) {

				newArgs <- makeArgs(global.model, allTerms[comb], argsOptions) #comb
				formulaList <- if(is.null(attr(newArgs, "formulaList"))) newArgs else
					attr(newArgs, "formulaList")

				if(!is.null(attr(newArgs, "problems"))) {
					print.warnings(structure(vector(mode = "list",
						length = length(attr(newArgs, "problems"))),
							names = attr(newArgs, "problems")))
				} # end if <problems>
				cl <- gmCall
				cl[names(newArgs)] <- newArgs
			} else isok <- FALSE # end if <subset, m.max >= nvar >= m.min>
		} #  end if(jComb != prevJComb)

		if(isok) {
			## --- Variants ---------------------------
			clVariant <- cl
			isok2 <- TRUE
			if(nVarying) {
				cvi <- variants[varComb + 1L, ]
				isok2 <- (hasSubset != 4L) || evalExprInEnv(subsetExpr, env = ssEnv,
					enclos = parent.frame(), comb = comb, `*nvar*` = nvar,
					cVar = variantsFlat[cvi])
				clVariant[varyingNames] <- fvarying[cvi]
			}

			if(isok2) {
				if(evaluate) {
					dotrace()

					qi <- qi + 1L
					queued[[(qi)]] <- list(call = clVariant, id = iComb)
				} else { # if !evaluate
					k <- k + 1L # all OK, add model to table
					rvlen <- length(ord)
					if(k > rvlen) {
						nadd <- min(resultChunkSize, ncomb - rvlen)
						#message(sprintf("extending result from %d to %d", rvlen, rvlen + nadd))
						addi <- seq.int(rvlen + 1L, length.out = nadd)
						calls[addi] <- vector("list", nadd)
						ord[addi] <- integer(nadd)
					}
					calls[[k]] <- clVariant
					ord[k] <- iComb
				}
			}
		} # if isok

		#if(evaluate && qi && (qi + nvariants > qlen || iComb == ncomb)) {
		if(evaluate && qi && (qi > qlen || (iComb + 1) == ncomb)) {
			qseq <- seq_len(qi)
			qresult <- .getRow(queued[qseq])
			utils::flush.console()

			if(!all(vapply(qresult, function(x) is.list(x) && "value" %in% names(x), FALSE)))
				stop("some results returned from cluster node(s) are malformed or NULL. \n",
					"This should not happen and indicates problems with ",
					"the cluster node", domain = "R-MuMIn")
			haveProblems <- logical(qi)

			nadd <- sum(sapply(qresult, function(x) inherits(x$value, "condition")
				+ length(x$warnings)))
			wi <- length(warningList)
			if(nadd) warningList <- c(warningList, vector(nadd, mode = "list"))

			# DEBUG: print(sprintf("Added %d warnings, now is %d", nadd, length(warningList)))

			for (i in qseq)
				for(cond in c(qresult[[i]]$warnings,
					if(inherits(qresult[[i]]$value, "condition"))
						list(qresult[[i]]$value))) {
						wi <- wi + 1L
						warningList[[wi]] <- if(is.null(conditionCall(cond)))
							queued[[i]]$call else conditionCall(cond)
						if(inherits(cond, "error")) {
							haveProblems[i] <- TRUE
							msgsfx <- "(model %d skipped)"
						} else
							msgsfx <- "(in model %d)"
						names(warningList)[wi] <- paste(conditionMessage(cond),
							 gettextf(msgsfx, queued[[i]]$id))
						attr(warningList[[wi]], "id") <- queued[[i]]$id
				}

			withoutProblems <- which(!haveProblems)
			qrows <- lapply(qresult[withoutProblems], "[[", "value")
			qresultLen <- length(qrows)
			rvlen <- nrow(rval)

			if(retNeedsExtending <- k + qresultLen > rvlen) {
				nadd <- min(max(resultChunkSize, qresultLen), ncomb - rvlen)
				rval <- rbind(rval, matrix(NA_real_, ncol = rvNcol, nrow = nadd),
					deparse.level = 0L)
				addi <- seq.int(rvlen + 1L, length.out = nadd)
				coefTables[addi] <- vector("list", nadd)
				calls[addi] <- vector("list", nadd)
				ord[addi] <- integer(nadd)
			}
			qseqOK <- seq_len(qresultLen)
			for(m in qseqOK) rval[k + m, retColIdx] <- qrows[[m]]
			ord[k + qseqOK] <- vapply(queued[withoutProblems], "[[", 1L, "id")
			calls[k + qseqOK] <- lapply(queued[withoutProblems], "[[", "call")
			coefTables[k + qseqOK] <- lapply(qresult[withoutProblems], "[[", "coefTable")
			k <- k + qresultLen
			qi <- 0L
		}
	} ### for (iComb ...)

	if(k == 0L) {
		if(length(warningList)) print.warnings(warningList)
		stop("the result is empty")
	}

	names(calls) <- ord
	if(!evaluate) return(calls[seq_len(k)])

	if(k < nrow(rval)) {
		i <- seq_len(k)
		rval <- rval[i, , drop = FALSE]
		ord <- ord[i]
		calls <- calls[i]
		coefTables <- coefTables[i]
	}

	if(nVarying) {
		varlev <- ord %% nVariants
		varlev[varlev == 0L] <- nVariants
		rval[, nVars + seq_len(nVarying)] <- variants[varlev, ]
	}

	rval <- as.data.frame(rval, stringsAsFactors = TRUE)
	row.names(rval) <- ord

	# Convert columns with presence/absence of terms to factors
	tfac <- which(!(allTerms %in% gmCoefNames))
	rval[tfac] <- lapply(rval[tfac], factor, levels = NaN, labels = "+")

	i <- seq_along(allTerms)
	v <- order(termsOrder)
	rval[, i] <- rval[, v]
	allTerms <- allTerms[v]
	colnames(rval) <- c(allTerms, varyingNames, extraNames, "df", lik$name, ICName)

	if(nVarying) {
		variant.names <- vapply(variantsFlat, asChar, "", width.cutoff = 20L)
		vnum <- split(seq_len(sum(vlen)), rep(seq_len(nVarying), vlen))
		names(vnum) <- varyingNames
		for (i in varyingNames) rval[, i] <-
			factor(rval[, i], levels = vnum[[i]], labels = variant.names[vnum[[i]]])

	}

	rval <- rval[o <- order(rval[, ICName], decreasing = FALSE), ]
	coefTables <- coefTables[o]

	rval$delta <- rval[, ICName] - min(rval[, ICName])
	rval$weight <- Weights(rval$delta)
    mode(rval$df) <- "integer"

	rval <- 
    structure(rval,
		model.calls = calls[o],
		global = global.model,
		global.call = gmCall,
		terms = structure(allTerms, interceptLabel = interceptLabel),
		rank = IC,
		beta = strbeta,
		call = {
            cl <- match.call(expand.dots = TRUE)
            cl[[1L]] <- as.symbol("dredge")
            cl
        },
		coefTables = coefTables,
		nobs = gmNobs,
		vCols = varyingNames, ## XXX: remove
        column.types = {
			colTypes <- c(terms = length(allTerms), varying = length(varyingNames),
				extra = length(extraNames), df = 1L, loglik = 1L, ic = 1L, delta = 1L,
				weight = 1L)
			column.types <- rep(1L:length(colTypes), colTypes)
			names(column.types) <- colnames(rval)
			lv <- 1L:length(colTypes)
			factor(column.types, levels = lv, labels = names(colTypes)[lv])
		},
        class = c("model.selection", "data.frame")
	)

	if(length(warningList)) {
		class(warningList) <- c("warnings", "list")
		attr(rval, "warnings") <- warningList
	}

	if (!is.null(attr(allTerms0, "random.terms")))
		attr(rval, "random.terms") <- attr(allTerms0, "random.terms")

	if(doParallel) clusterCall(cluster, "rm",
		list = c(".pdredge_process_model", "pdredge_props"), envir = .GlobalEnv)


	return(rval)
} ######


`pdredge_process_model` <-
function(modv, envir = get("pdredge_props", .GlobalEnv)) {
	### modv == list(call = clVariant, id = modelId)
	result <- tryCatchWE(eval(modv$call, get("gmEnv", envir)))
	if (inherits(result$value, "condition")) return(result)

	fit1 <- result$value
	if(get("nExtra", envir) != 0L) {
		extraResult1 <- get("applyExtras", envir)(fit1)
		nExtra <- get("nExtra", envir)
		if(length(extraResult1) < nExtra) {
			tmp <- rep(NA_real_, nExtra)
			tmp[match(names(extraResult1), get("extraResultNames", envir))] <-
				extraResult1
			extraResult1 <- tmp
		}
	} else extraResult1 <- NULL
	ll <- .getLik(fit1)$logLik(fit1)

	#mcoef <- matchCoef(fit1, all.terms = get("allTerms", envir),
	# beta = get("beta", envir), allCoef = TRUE)
	mcoef <- eval(get("matchCoefCall", envir))

	list(value = c(mcoef, extraResult1, df = attr(ll, "df"), ll = ll,
		ic = get("IC", envir)(fit1)),
		nobs = nobs(fit1),
		coefTable = attr(mcoef, "coefTable"),
		warnings = result$warnings)

}
