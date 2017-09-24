\encoding{utf-8}
\name{loo}
\alias{loo}
\alias{loo.default}
\alias{loo.lm}
\title{Leave-one-out cross-validation}
\usage{
loo(object, type = c("loglik", "rmse"), ...)
}
\arguments{
\item{object}{a fitted object model, currently only \code{lm}/\code{glm} is
accepted.}

\item{type}{the criterion to use, given as a character string,
either \code{"rmse"} for Root-Mean-Square Error or \code{"loglik"}
for log-likelihood.}

\item{...}{other arguments are currently ignored.}
}
\value{
\code{loo} returns a single numeric value of \acronym{RMSE} or
    mean log-likelihood.
}
\description{
Computes the \acronym{RMSE}/log-likelihood based on leave-one-out cross-validation.
}
\details{
Leave-one-out cross validation is \var{K}-fold cross validation, with \var{K}
equal to the number of data points in the set \var{N}. For all \var{i} from 1 
to \var{N}, the model is fitted to all the data except for point \var{i} and a 
prediction is made for that point. The average error is computed and used to 
evaluate the model.
}
\examples{
fm <- lm(y ~ X1 + X2 + X3 + X4, Cement)
loo(fm, type = "l")
loo(fm, type = "r")

## Compare LOO_RMSE and AIC/c
options(na.action = na.fail)
dd <- dredge(fm, rank = loo, extra = list(AIC, AICc), type = "rmse")
plot(loo ~ AIC, dd, ylab = expression(LOO[RMSE]), xlab = "AIC/c")
points(loo ~ AICc, data = dd, pch = 19)
legend("topleft", legend = c("AIC", "AICc"), pch = c(1, 19))
}
\author{
Kamil Barto\enc{ń}{n}, based on code by Carsten Dormann
}
\references{
Dormann, C. et al. (\emph{in prep.}) Model averaging in ecology: evidence,
approach examples.
}
\keyword{models}
