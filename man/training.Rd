\name{training}

\alias{training}

\title{Performs the HVDM training step and returns a list containing the results}
\description{
  This method performs the training step of the HVDM algorithm. It returns a list that will then be used in the subsequent
  screening step. 
}
\usage{
training(eset,genes,transforms,degrate,actname,pdata,forcetransforms)
}
\arguments{
  \item{eset}{ an ExpressionSet object (Biobase) }
  \item{genes}{ a vector containg the gene identifiers of the training genes }
  \item{transforms}{ a vector containing the kintetic parameter identifiers that have to be transformed during optimisation (optional) }
  \item{degrate}{ value of the anchoring gene degradation rate, expressed in inverse unit time (optional) }
  \item{pdata}{ a data frame, phenoData to be used for the training (optional)}
  \item{actname}{ name of the transcription factor of interest (optional)}
  \item{forcetransforms}{ Boolean, whether the transformation in argument transforms have to be applied }
}
\details{ 
The first entry in the genes vector is the anchoring gene. This means that the sensitivity (Sj) for this genes is set
at 1.0 by default and that if a degradation rate is supplied it applies to that gene.

An exponential transform is set by default for both the basal (Bj) and degradation (Dj) rates (through the transforms argument). This
forces the values for both these parameters to be positive. It also helps to reach a better fit. To turn this off, set the \code{forcetransforms}
switch to FALSE. Even in this case the degradation rate will not be allowed to take non-positive values as it causes problems with the
differential operator used internally. The value in the vector indicates the parameter to be transformed: "Bj": basal rate of transcription, "Sj": sensitivity, "Dj": degrdation rate.
The entry label indicates the transform to be applied; presently, only log-tranforms are implemented (ie "exp").

The \code{degrate} argument is optional, but it is recommended to provide the algorithm with an externally measured degradation rate, as this
greatly improves the accuracy and robustness of the outcome.

The \code{pdata} argument is also optional. By default the method will use the phenoData contained in the expression set.
This argument can be used for excluding a time point, or an entire replicate. To extract the phenoData from the expression set, 
use \code{dataframe<-pData(eset)}. The \code{dataframe} object obtained can then be manipulated as desired.

The default name of the transcription factor is "trfact1". 
}

\value{
	a list containing the results (see documentation for more details).
}
\references{ 
M. Barenco, D. Tomescu, D. Brewer, R. Callard, J. Stark, M. Hubank (2006) Ranked predictions
of p53 targets using Hidden Variable Dynamic Modelling. \emph{Genome Biology}, \bold{V7(3)}, R25.
}
\author{ Martino Barenco }

\note{ It is recommended to run the HVDMcheck method before running this command.
}

\seealso{\code{\link{HVDMcheck}},\code{\link{screening}},\code{\link{fitgene}},\code{\link{HVDMreport}}}

\examples{
data(HVDMexample)
tHVDMp53<-training(eset=fiveGyMAS5,genes=p53traingenes,degrate=0.8,actname="p53")
}

\keyword{ ts }


