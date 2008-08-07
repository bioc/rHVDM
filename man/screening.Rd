\name{screening}

\alias{screening}

\title{Fits the optimal kinetic parameter values for several genes.}
\description{
  This method fits the three kinetic parameter values for each gene in a user-supplied vector. 
  It returns a list containing the results.
}
\usage{
screening(eset,genes,HVDM,transforms,cl1zscorelow,cl1modelscorehigh,cl1degraterange)
}

\arguments{
  \item{eset}{ an ExpressionSet object (Biobase) }
  \item{genes}{ a vector containing the genes identifiers to be screened (in character format) }
  \item{transforms}{ a vector containing the kintetic parameter identifiers that have to be transformed during optimisation (optional) }
  \item{HVDM}{ a HVDM object (see details) }
  \item{cl1zscorelow}{ the sensitivity Z-score cutoff value for a gene to be classified as a putative target }
  \item{cl1modelscorehigh}{ the model score cutoff value for a gene to be classified as a putative target }
  \item{cl1degraterange}{ the degradation rate bounds applied for a gene to be classified as a putative target }
}

\details{ 
This screening step can only be applied if a \code{training()} step has already been run. 
The output to the \code{training()} step can be given as the "HVDM" argument. A fit of each gene in the "genes" argument
is then performed. 

Alternatively an output to a previously run \code{screening()} step can be supplied as an
"HVDM" argument. In this case, the fit is not performed once again. Feeding a previous output
of \code{screening()} to that same function again is useful only if the various bounds altered in the secon run.
In the case this option is used, the "eset" and "genes" arguments do not have to be supplied (they will
be ignored anyway).

The output of that function is a list containing the results. The relevant data frame is in the "results" member of the output. 
Putative targets can be identified using the "class1" field of that data frame (see example).

Bounds determining whether a gene is a target of the transcription factor under review have to be supplied. They are:

- cl1zscorelow: lower bound for the sensitivity Z-score (default: 2.5)

- cl1modelscore: upper bound for the model score (default: 100.0). 
This default will have to be changed in most cases. As a rule of thumb,  
5x the model score for the genes in the training set can be used.

- cl1degraterange: lower and upper bounds for the degradation rate (default: c(0.05,5.0)). This is to 
exclude those genes with an absurd degradation rate, measured in (unit time)\^(-1). In our example the
unit time is an hour. In the case the unit time is different, these default bounds will have to be altered accordingly.

An exponential transform is set by default for both the basal (Bj) and degradation (Dj) rates (through the transforms argument). This
forces the values for both these parameters to be positive. It also helps to reach a better fit. To turn this off let transforms=c(). Even in this case the degradation rate will not be allowed to take non positive values as it causes problems with the
differential operator used internally. The value in the vector indicates the parameter to be transformed: "Bj": basal rate of transcription, "Sj": sensitivity, "Dj": degradation rate.
The entry label indicates the transform to be applied; presently, only log-tranforms are implemented (ie "exp").
}

\value{
	a list containing the results (see documentation for more details).
}
\references{ 
M. Barenco, D. Tomescu, D. Brewer, R. Callard, J. Stark, M. Hubank (2006) Ranked predictions
of p53 targets using Hidden Variable Dynamic Modelling. \emph{Genome Biology}, \bold{V7(3)}, R25.
}
\author{ Martino Barenco }

\note{ Obviously, the expression set given as a \code{eset} argument has to be the same as the one used
for the training set.
}

\seealso{\code{\link{training}},\code{\link{HVDMreport}},\code{\link{fitgene}}}

\examples{
data(HVDMexample)
tHVDMp53<-training(eset=fiveGyMAS5,genes=p53traingenes,degrate=0.8,actname="p53")
screenp53<-screening(eset=fiveGyMAS5,genes=genestoscreen[1:10],HVDM=tHVDMp53)

#extracting a list of putative p53 targets
p53targets<-screenp53$results[screenp53$results$class1,]

#shifting the goal posts
screenp53B<-screening(HVDM=screenp53,cl1zscorelow=3.5)
}
\keyword{ ts }

