\name{rppa.load}
\alias{rppa.load}
\title{Load a RPPA slide from MIRACLE}
\usage{
rppa.load(connection = NULL, barcode = NA, slideIndex = NA,
  securityToken = NA, baseUrl = "http://localhost:8080/MIRACLE/spotExport/",
  filter.diameter = T, filter.neg.values = T, filter.flag = T,
  apply.shifts = T)
}
\arguments{
  \item{connection}{RCurl connection object obtained
  through \code{\link{rppa.authenticate}}}

  \item{barcode}{Barcode of the slide}

  \item{slideIndex}{Database ID of the slide in MIRACLE}

  \item{securityToken}{UUID of the slide used for secure
  access without \code{connection}}

  \item{baseUrl}{URL pointing to the spotExport controller
  of the used MIRACLE instance}

  \item{filter.diameter}{if true filter spots with diameter
  > 200 (set to NA)}

  \item{filter.negative}{if true filter negative values
  (set to NA)}

  \item{filter.flag}{if true filter flagged values, e.g.
  where flag property is not 0 (set to NA)}

  \item{apply.shifts}{if true apply horizontal and vertical
  shifts}
}
\description{
This function loads a RPPA slide from MIRACLE using one of
the following as input: a connection object obtained with
\code{\link{rppa.authenticate}} together with either
barcode or slideIndex or a securityToken where user
authentication is completly omitted. In addition we need to
specify MIRACLE's URL as baseUrl, specify signal filter
options and optionally correct for block shifts in the
process.
}
\examples{
baseUrl <- "http://192.168.0.1:8080/MIRACLE/spotExport/"
conn <- rppa.authenticate(user="mlist", password="xxx", baseUrl=baseUrl)
rppa.load(connection=conn, slideIndex=1, baseUrl=baseUrl)
rppa.load(connection=conn, barcode="ABCD")
rppa.load(securityToken="abcd-xzft-djfk-2345")
}
\keyword{io}
\keyword{json}
\keyword{load}
\keyword{read}

