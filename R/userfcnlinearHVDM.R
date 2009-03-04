#user function(s) for rHVDM (linear model only)

training<-function(eset,genes,transforms=c("exp"="Dj","exp"="Bj"),degrate,actname="trfact1",pdata,forcetransforms=TRUE){
	if(missing(pdata)) pdata<-pData(eset)
	if(forcetransforms) tHVDM<-.initialisetrainingHVDM(eset=eset,trainingset=genes,degrate=degrate,pdata=pdata,trfname=actname,transforms=transforms)
	else tHVDM<-.initialisetrainingHVDM(eset=eset,trainingset=genes,degrate=degrate,pdata=pdata,trfname=actname,transforms=c())
	tHVDM<-.tcfirstguess(tHVDM=tHVDM)
	results<-.optim(HVDM=tHVDM)
	tHVDM<-.freeparsevaluate(HVDM=tHVDM,x=results$par)
	results$par<-NULL  #this is not needed any more in the HVDM object
	tHVDM$results<-results
	tHVDM$scores<-.scorout(tHVDM)
	tHVDM$itgenes<-.screenall(eset,genes,tHVDM)
	tHVDM$eset<-deparse(substitute(eset))
	tHVDM$type<-c("training")
	tHVDM
}

fitgene.lin<-function(eset,gene,tHVDM,transforms=c("exp"="Dj","exp"="Bj"),firstguess){
	sHVDM<-.fitgene(eset=eset,gene=gene,tHVDM=tHVDM,transforms=transforms,firstguess=firstguess)
	sHVDM$tset<-tHVDM$itgenes
	sHVDM$eset<-deparse(substitute(eset))
	sHVDM$tHVDMname<-deparse(substitute(tHVDM))
	sHVDM$type<-c("indgene")
	sHVDM
}

screening<-function(eset,genes,HVDM,transforms=c("exp"="Dj","exp"="Bj"),cl1zscorelow=2.5,cl1modelscorehigh=100.0,cl1degraterange=c(0.01,5.0)){
	if(HVDM$type=="training"){
		reslis<-.screenall(eset=eset,genes=genes,tHVDM=HVDM,transforms=transforms)
	}
	else if(HVDM$type=="screening"){
		reslis<-HVDM$results
		HVDM<-HVDM$tHVDM
	}
	class1<-((reslis$sens_z_score>=cl1zscorelow) & (reslis$model_score<=cl1modelscorehigh) & (reslis$degradation>=cl1degraterange[1]) & (reslis$degradation<=cl1degraterange[2]))
	reslis$class1<-class1
	ordered<-order(-reslis$sens_z_score)
	results<-reslis[ordered,]
	res<-vector("list")
	res$results<-results
	res$tHVDM<-HVDM
	res$transforms<-transforms
	bounds<-vector("list")
	bounds$zscore<-cl1zscorelow
	bounds$modelscore<-cl1modelscorehigh
	bounds$degrate<-cl1degraterange
	res$class1bounds<-bounds
	res$type<-"screening"
	res
}