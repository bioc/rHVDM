\name{estimerrors}

\alias{estimerrors}

\title{computes the standard deviation of the measurement error using pre-calculated tables specifict to the plattform or user-defined table}
\description{
  This method takes as main input and ExpressionSet object, and using pre-calculated tables contained in rHVDM, it returns the same object with an extra slot that contains the standard deviation of the measurement error. The latter is indispensable to then run rHVDM proper. The table can be user-defined (using the \code{errtable} parameter). To check which tables are stored in rHVDM, just run the command with the \code{eset} parameter left empty, it will return a list of the supported rHVDM plattforms. "Plattform" should be understood in a generic way here and comprises both the actual physical plattform and the mathematical method (summary) used to handle the primary data (for example, GCRMA or MAS5 for affymetrix microarrays would use two different tables). All parameters of this function are optional, see below for details on the effect of omitting one or more of them. The method assumes the the amount of error is dependent on both on the particular microarray and the signal intensity. For more details on how it is computed, please refer to the textual description of rHVDM.
}
\usage{
estimerrors(eset,plattid,refchips,errtable)
}
\arguments{
  \item{eset}{ an object of class ExpressionSet (Biobase), every time this parameter is omitted the function returns a list of supported plattforms (identifier+plattform description).
	}
	\item{plattid}{ an optional argument (character format) allowing to specify the plattform identifiers. These identifiers are rHVDM specific. Any identifier will be overriden by a table given as input of the \code{errtable} argument. In future versions of rHVDM, if this argument and the \code{errtable} are missing, a search will be done among the supported plattforms to find a matching plattform, based on the individual genes identifiers. For now, omission of both arguments will be conductive to outputting a list of supported plattforms.
	}
	\item{refchips}{ a vector of names or column indexes compatible with the ExpressionSet useds (note that this compatibility is not verified by the function, it is up to the user to supply compatible array names). Some microarray might be known to be of not very good quality and although the measurement error will be computed for these chips it is better to leave them aside during early stages of the computation. The list to be input should contain all the "good" microarrays. This parameter is optional and if omitted, all microarrays are used.
	}
	\item{errtable}{ a nx2 table in matrix format. The first column of the matrix should contain the natural log of signal intensity in ascending order and the second the corresponding variance.
	}
}
\details{ 
	Local interpolation and individual array normalisation are used to estimate the standard deviation of the measurement error for each individual transcript on each microarray. The precise method is laid out in the paper cited below, with some minor modifications that are explained in the textual description of rHVDM.	
}

\value{
	Normally an updated eset is returned. In case some crucial element is missing, the original eset is returned. If the latter is missing too a list of supported plattform is returned (thus, when using this function, it is better to be careful of what lies at the LHS of the assignement arrow).
}
\references{ 
M. Barenco, D. Tomescu, D. Brewer, R. Callard, J. Stark, M. Hubank (2006) Ranked predictions
of p53 targets using Hidden Variable Dynamic Modelling. \emph{Genome Biology}, \bold{V7(3)}, R25.
}
\author{ Martino Barenco }

\note{ The HTML report is generated in the working directory.
}

\seealso{\code{\link{HVDMcheck}},\code{\link{training}},\code{\link{screening}},\code{\link{fitgene}}}

\examples{
data(HVDMexample)
tHVDMp53<-training(eset=fiveGyMAS5,genes=p53traingenes,degrate=0.8,actname="p53")
#leave that as an example for now
#HVDMreport(tHVDMp53)
#HVDMreport(tHVDMp53, name="p53training")
}
\keyword{ IO }
\keyword{ ts }
