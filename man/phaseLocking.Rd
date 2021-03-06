% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phase_locking.R
\name{phaseLocking}
\alias{phaseLocking}
\title{pahseLocking calculation for two vectors}
\usage{
phaseLocking(x, y, n = 1, m = 1)
}
\arguments{
\item{x}{- a vector, should be a time series vector}

\item{y}{- a vector, should be a time series vector}

\item{n}{- time period tatio}

\item{m}{- time period tatio}
}
\value{
list with 1. entropy rho 2. gamma  3. Strobo index lambda 4. Strobo index angle
}
\description{
pahseLocking is written for testing purposes.
if we are calculating phase locking for a large matrix,
running it per pair is slow
}
\examples{
x <- sin(1:200)
y <- sin(2*1:200)
n=1
m=2
phaseLocking(x, y, n=1, m=2)
}
\references{
[1]
}

