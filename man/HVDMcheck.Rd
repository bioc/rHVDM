\name{HVDMcheck}

\alias{HVDMcheck}

\title{Checks that an object of class exprSet is compliant with HVDM}
\description{
  This method issues warnings about the compliance of the exprSet object to be used with rHVDM.
  It just issues warnings and does not attempt to correct the problems encountered. the \code{pdata} argument
  is optional; if missing, the method will check the phenodata that is inside
  the exprSet.
}
\usage{
HVDMcheck(eset,pdata)
}
\arguments{
  \item{eset}{ an exprSet object (Biobase) }
  \item{pdata}{ a dataframe (optional argument) }
}
\details{ The checks performed are: A) the phenoData has all the required fields. B) the time field is numeric, has
  a zero time point, does not have negative or repeated values. C) the row names in the phenoData 
  are consistent with the columns names found in the exprSet. D) the data is not log-transformed.
  E) Standard deviations have been supplied in the exprSet, and they are non-negative.
}

\value{
	no value returned.  
}
\references{ 
M. Barenco, D. Tomescu, D. Brewer, R. Callard, J. Stark, M. Hubank (2006) Ranked predictions
of p53 targets using Hidden Variable Dynamic Modelling. \emph{Genome Biology}, \bold{V7(3)}, R25.
}
\author{ Martino Barenco }

\note{ If the \code{pdata} argument is ignored (it can be), 
the method checks the phenoData that is inside the exprSet.}

\seealso{\code{\link{training}}}

\examples{
data(HVDMexample)
HVDMcheck(fiveGyMAS5)
}
\keyword{ IO }
\keyword{ ts }

