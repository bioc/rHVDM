\name{training.nl}

\alias{training.nl}

\title{Performs the HVDM training step and returns a list containing the results}
\description{
  This method performs the nonlinear training step of the HVDM algorithm. It returns a list that will then be used in the subsequent
  screening step. 
}
\usage{
training.nl(inputHVDM,transforms,constraints,forcetransforms,genemodels,firstguess)
}
\arguments{
  \item{inputHVDM}{ an HVDM training object }
  \item{transforms}{ a vector containing the kinetic parameter identifiers that have to be transformed during optimisation (optional) }
  \item{constraints}{ "known" values for the activator signal }
  \item{forcetransforms}{ Boolean, whether the transformation in argument transforms have to be applied }
  \item{genemodels}{the type of model that has to be used for each gene}
  \item{firstguess}{a first guess for all the parameters}
}
\details{ 
Contrary to the linear training function (without .nl suffix), this function takes as main input another training object (either a linear training object or a non-linear one, where all the genes have been fitted witha MIchelis-Menten model ("MM").

An exponential transform is set by default for the basal (Bj) and degradation (Dj) rates, as well as for the kinetic parameters (Vj and Kj) in the production function (through the transforms argument). This forces the values for these parameters to be positive For the exponent Nj, the \code{expp1} function ensures it is greater than one, when the hill formulation is used. To turn this off, set the \code{forcetransforms} switch to FALSE. Even in this case the degradation rate will not be allowed to take non-positive values as it causes problems with the differential operator used internally. The value in the vector indicates the parameter to be transformed: "Bj": basal rate of transcription, "Sj": sensitivity, "Dj": degrdation rate, etc.. The entry label indicates the transform to be applied.

The \code{constraints} argument is used to specify values for know values of some of the parameters, typically the strength of the activator for some time points. See example below for the syntax.

The \code{genemodels} argument is compulsory and used to specify the model used for each individual gene. MM means Michaelis-Menten model whereas hill means a hill function is used. The general form of the production function is $B_j+V_j*f(t)^N_j/(K_j^N_j+f(t)^N_j)$ where Nj=1 for the MM model and Nj>1 in the hill case. See below for the syntax of the input.

The {firstguess} argument is not in use yet and will offer the possibility to enter a first guess for the fitting.
}

\value{
	a list containing the results.
}
\references{ 
M. Barenco, D. Tomescu, D. Brewer, R. Callard, J. Stark, M. Hubank (2006) Ranked predictions
of p53 targets using Hidden Variable Dynamic Modelling. \emph{Genome Biology}, \bold{V7(3)}, R25.
}
\author{ Martino Barenco }

\seealso{\code{\link{HVDMcheck}},\code{\link{screening.nl}},\code{\link{fitgene}},\code{\link{HVDMreport}},\code{\link{training}}}

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
}

\keyword{ ts }


