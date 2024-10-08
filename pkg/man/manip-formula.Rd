\name{Formula manipulation}
\alias{simplify.formula}
\alias{expand.formula}

\encoding{utf-8}

\title{Manipulate model formulas}
\description{
\code{simplify.formula} rewrites a \code{formula} into shorthand notation.
Currently only the factor crossing operator \code{*} is applied, so an
expanded expression such as \code{a+b+a:b} becomes \code{a*b}.
\code{expand.formula} does the opposite, additionally expanding other
expressions, i.e. all nesting (\code{/}), grouping and \code{^}. }

\usage{
simplify.formula(x)
expand.formula(x)
}

\arguments{
  \item{x}{a \code{formula} or an object from which it can be extracted (such as
  a fitted model object). }
}
\author{Kamil Barto\enc{ń}{n}}

\seealso{
\lxcode{formula}{stats}

\lcode{delete.response}, \lcode{drop.terms}, and 
\lcode{reformulate}

}


\examples{

simplify.formula(y ~ a + b + a:b + (c + b)^2)
simplify.formula(y ~ a + b + a:b + 0)

expand.formula(~ a * b)

}

\keyword{manip}
