#user functions for rHVDM (nonlinear model)


fitgene.nl<-function(eset,gene,tHVDM,transforms=c("exp"="Dj","exp"="Bj","exp"="Kj","expp1"="Nj"),firstguess,model="MM"){
	sHVDM<-.fitgene.nl(eset=eset,gene=gene,tHVDM=tHVDM,transforms=transforms,firstguess=firstguess,model=model)
	sHVDM$tset<-tHVDM$itgenes
	sHVDM$eset<-deparse(substitute(eset))
	sHVDM$tHVDMname<-deparse(substitute(tHVDM))
	sHVDM$type<-c("indgene")
	sHVDM
}

training.nl<-function(inputHVDM,transforms=c("exp"="Dj","exp"="Bj","exp"="Kj","expp1"="Nj"),constraints,
						forcetransforms=TRUE,genemodels,firstguess){
	if(inputHVDM$type=="training"){
		#In case the input HVDM object is a linear fit, a first fit with MM throughout is performed 
		#through use of a special function hidden from the end user
		iptHVDM<-.training.lintoMM(inputHVDM=inputHVDM,transforms=transforms,constraints=constraints,forcetransforms=forcetransforms)
	}
	else{
		iptHVDM<-inputHVDM
	}
	#here have to test whether the genemodels are all MM or absent, in which case, subsequent treatment are skipped as the main job has already been done
	if(missing(genemodels)){
		tHVDMnl<-iptHVDM
	}
	else{
		if(prod(genemodels=="MM")){
			tHVDMnl<-iptHVDM
		}
		else{
			if(forcetransforms) tHVDMnl<-.initialisetrainingHVDM.nl(HVDM=iptHVDM,constraints=constraints,
																transforms=transforms,genemodels=genemodels)
			else tHVDMnl<-.initialisetrainingHVDM.nl(HVDM=iptHVDM,constraints=constraints,transforms=c(),genemodels=genemodels)
			#here introduce the possibility to input ad hoc first guesses (input a parameter in training.nl)
			if (missing(firstguess)){
				tHVDMnl<-.tcfirstguess.nl(tHVDM=tHVDMnl,inputHVDM=iptHVDM)
			}
			else{
				bidon="INPUTFIRSTGUESSHERE"
			}
			tHVDMnl$type<-"training.nl"
			results<-.optim(HVDM=tHVDMnl)
			tHVDMnl<-.freeparsevaluate(HVDM=tHVDMnl,x=results$par)
			results$par<-NULL  #this is not needed any more in the HVDM object
			tHVDMnl$results<-results
			tHVDMnl$scores<-.scorout(tHVDMnl)
		}
	}
	ligenes<-row.names(inputHVDM$par$genemap)
	tHVDMnl$itgenes<-.screenall.nl(eset=eval(parse(text=tHVDMnl$eset)),ligenes,tHVDMnl)
	tHVDMnl
}

screening.nl<-function(eset,genes,HVDM,transforms=c("exp"="Dj","exp"="Bj",exp="Kj",expp1="Nj"),cl1zscorelow=2.5,cl1modelscorehigh=100.0,cl1degraterange=c(0.01,5.0),criterion="BIC"){
	if(HVDM$type=="training.nl"){
		reslis<-.screenall.nl(eset=eset,genes=genes,tHVDM=HVDM,transforms=transforms,criterion=criterion)
	}
	else if(HVDM$type=="screening.nl"){
		reslis<-HVDM$results
		HVDM<-HVDM$tHVDM
	}
	class1<-((reslis$Vj_z_score>=cl1zscorelow) & (reslis$model_score<=cl1modelscorehigh) & (reslis$degradation>=cl1degraterange[1]) & (reslis$degradation<=cl1degraterange[2]))
	reslis$class1<-class1
	ordered<-order(-reslis$Vj_z_score)
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
	res$type<-"screening.nl"
	res
}