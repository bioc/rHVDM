\name{fitgene.nl}

\alias{fitgene.nl}

\title{Fits the optimal kinetic parameter values for a particular gene and model.}
\description{
  This method fits the three kinetic parameter values for a particular gene. It returns a list containing the results.
}
\usage{
fitgene.nl(eset,gene,tHVDM,transforms,firstguess,model)
}

\arguments{
  \item{eset}{ an ExpressionSet object (Biobase) }
  \item{gene}{ the gene identifier in character format }
  \item{transforms}{ a vector containing the kinetic parameter identifiers that have to be transformed during optimisation (optional) }
  \item{tHVDM}{ the output of the training set }
  \item{firstguess}{ first guess for the fitting (optional, see details)}
  \item{model}{ model specification of the production function in character format" MM (Michelis-Menten, default) or hill (hill function) }
}
  
\details{ 
An exponential transform is set by default for both the basal (Bj), degradation (Dj) and Kj rates (through the transforms argument). This
forces the values for both these parameters to be positive. It also helps to reach a better fit. To turn this off let transforsm=c(). Even in this case the degradation rate will not be allowed to take non positive values as it causes problems with the
differential operator used internally. The value in the vector indicates the parameter to be transformed: "Bj": basal rate of transcription, "Vj": sensitivity, "Dj": degrdation rate.
The entry label indicates the transform to be applied.

This \code{fitgene.nl()} step can only be applied after a \code{training()} step. The output to the \code{training.nl()} step has to be fed through 
the \code{tHVDM} argument.

The \code{firstguess} argument is optional (a first guess is generated internally by default). 
However a first guess can be supplied by the user which can take several forms. 
It can either be a vector with three entries containing a first guess for the basal rate, 
the sensitivity, the degradation rate (in that order). 
Alternatively, another output from the \code{fitgene()} function (for example from a gene that
 has a similar expression profile) can be supplied as a \code{firstguess} argument.

The \code{model} argument is only used to choose the model (MM: Michelis-Menten, hill: Hill function).
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
for the training step.
}

\seealso{\code{\link{training}},\code{\link{screening}},\code{\link{HVDMreport}}}

\examples{
#load data and fit a linear model
data(HVDMexample)
rm(fiveGyMAS5)
data(HVDMexample2)
tp532<-training(eset=twodosesMAS5,genes=p53traingenes,degrate=0.8)
#formulate constraints
CONSTRAINTS<-c(350,35)
names(CONSTRAINTS)<-c("trfact1.5Gy.1.4","trfact1.hGy.1.4")
#specify individual gene models
GENEMODELS<-rep("MM",5)
GENEMODELS[c(1,3)]<-"hill"
names(GENEMODELS)<-p53traingenes
#run the model
tp53hyb<-training.nl(inputHVDM=tp532,constraints=CONSTRAINTS,genemodels=GENEMODELS)

#fit an individual gene
sHVDMex<-fitgene.nl(eset=twodosesMAS5,gene="213293_s_at",tHVDM=tp53hyb)
}

\keyword{ ts }


