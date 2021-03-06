% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allGenerics.R, R/allMethods.R
\docType{methods}
\name{chmGetURL}
\alias{chmGetURL}
\alias{chmGetURL,character-method}
\alias{chmGetURL,ngchm-method}
\title{Get the URL for an installed NGCHM.}
\usage{
chmGetURL(chm, ...)

\S4method{chmGetURL}{character}(chm, server = NULL, ...)

\S4method{chmGetURL}{ngchm}(chm, server = NULL, ...)
}
\arguments{
\item{chm}{The NGCHM for which the URL is required.}

\item{...}{Ignored.}

\item{server}{The server on which to view the NGCHM}

\item{server}{The server on which to view the NGCHM}
}
\description{
Return the URL for accessing the specified Next Generation Clustered Heat Map (NGCHM) on the specified server.
}
\seealso{
ngchmServer-class

ngchm-class
}

