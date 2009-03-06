#functions for nonlinear HVDM

.screenall.nl<-function(eset,genes,tHVDM,transforms=c("exp"="Dj","exp"="Bj","exp"="Kj","expp1"="Nj"),criterion="BIC"){
	n<-length(genes)
	res<-data.frame(row.names=genes,
					model_score=rep(0,n),
					lm_message=rep(0,n),
					npar=rep(0,n),
					Bj=rep(0,n),
					Vj=rep(0,n),
					Kj=rep(0,n),
					Nj=rep(0,n),
					Dj=rep(0,n),
					Vj_z_score=rep(0,n),
					hess_rank=rep(0,n),
					Bj_d=rep(0,n),
					Vj_d=rep(0,n),
					Kj_d=rep(0,n),
					Nj_d=rep(0,n),
					Dj_d=rep(0,n),
					Bj_u=rep(0,n),
					Vj_u=rep(0,n),
					Kj_u=rep(0,n),
					Nj_u=rep(0,n),
					Dj_u=rep(0,n),
					altscore=rep(0,n),
					linscore=rep(0,n)
					)
	for(gene in genes){
		print(gene)
		res[gene,]<-.genescreen.nl(eset=eset,gene=gene,tHVDM=tHVDM,template=res[gene,],transforms=transforms,criterion=criterion)
	}
	res
}


.genescreen.nl<-function(eset,gene,tHVDM,template,transforms=c("exp"="Dj","exp"="Bj","exp"="Kj","expp1"="Nj"),criterion="BIC"){
	#wrapper function for screening individual genes [tbc] works only for linear model
	res<-template
	sHVDM<-.fitgene.best(eset=eset,gene=gene,tHVDM=tHVDM,transforms=transforms,criterion="BIC")
	res[gene,]$linscore<-sHVDM$linscore
	res[gene,]$model_score<-sHVDM$scores$total
	res[gene,]$lm_message<-sHVDM$results$message
	hout<-.hessstuff.nl(sHVDM)
	centfit<-sHVDM$par$parameters[sHVDM$distribute$free]
	npar<-length(centfit)
	res[gene,]$npar<-npar
	res[gene,]$Bj<-centfit[1]
	res[gene,]$Vj<-centfit[2]
	res[gene,]$Kj<-centfit[3]
	res[gene,]$Dj<-centfit[4]
	res[gene,]$Nj<-1
	if(npar>4){
		res[gene,]$Nj<-centfit[4]
		res[gene,]$Dj<-centfit[5]
		res[gene,]$altscore<-sHVDM$score$info_MM$score
	}
	else{
		res[gene,]$altscore<-sHVDM$score$info_hill$score
	}
	if(class(hout)=="list"){
		res[gene,]$Vj_z_score<-(centfit[2]*2*1.96)/(hout[[2]][2]-hout[[1]][2])
		res[gene,]$hess_rank<-npar
		res[gene,]$Bj_d<-hout[[1]][1]
		res[gene,]$Vj_d<-hout[[1]][2]
		res[gene,]$Kj_d<-hout[[1]][3]
		res[gene,]$Nj_d<-1
		res[gene,]$Dj_d<-hout[[1]][4]
		res[gene,]$Bj_u<-hout[[2]][1]
		res[gene,]$Vj_u<-hout[[2]][2]
		res[gene,]$Kj_u<-hout[[2]][3]
		res[gene,]$Nj_u<-1
		res[gene,]$Dj_u<-hout[[2]][4]
		if(npar>4){
			res[gene,]$Nj_d<-hout[[1]][4]
			res[gene,]$Dj_d<-hout[[1]][5]
			res[gene,]$Nj_u<-hout[[2]][4]
			res[gene,]$Dj_u<-hout[[2]][5]
		}
	}
	else{
		res[gene,]$Vj_z_score<--1
		res[gene,]$hess_rank<-hout
		res[gene,]$Bj_d<-centfit[1]
		res[gene,]$Vj_d<-centfit[2]
		res[gene,]$Kj_d<-centfit[3]
		res[gene,]$Nj_d<-1
		res[gene,]$Dj_d<-centfit[4]
		res[gene,]$Bj_u<-centfit[1]
		res[gene,]$Vj_u<-centfit[2]
		res[gene,]$Kj_u<-centfit[3]
		res[gene,]$Nj_u<-1
		res[gene,]$Dj_u<-centfit[4]
		if(npar>4){
			res[gene,]$Nj_d<-centfit[4]
			res[gene,]$Dj_d<-centfit[5]
			res[gene,]$Nj_u<-centfit[4]
			res[gene,]$Dj_u<-centfit[5]
		}
	}
	res
}

.hessstuff.nl<-function(sHVDM){#NB works also for nonlinear model
	#returns the individual (raw) sigmas and correlations/ zero if hessian is singular
	H<-sHVDM$results$hessian
	dimen<-dim(H)[1]
	hr<-.hessrank(H)
	if(hr<dimen) res<-hr
	else{
		vcov<-solve(H)
		sdev<-1.96*diag(vcov)^0.5  #95% confidence interval
		work<-sHVDM
		central<-.exportfree(work)
		nams<-names(central)
		allupbnd<-(.importfree(HVDM=sHVDM,x=central[nams]+sdev[nams]))$par$parameters[nams]
		alllobnd<-(.importfree(HVDM=sHVDM,x=central[nams]-sdev[nams]))$par$parameters[nams]
		res<-list(alllobnd,allupbnd)
	}
	res
}

.hessrank<-function(H){
	#returns the rank of a matrix using singular values decomposition
	cns<-svd(H)$d
	cns<-cns/cns[1]
	sum(cns>1e-15)
}

#screening step

.fitgene.best<-function(eset,gene,tHVDM,transforms=c("exp"="Dj","exp"="Bj","exp"="Kj","expp1"="Nj"),criterion="BIC"){
	#run an individual gene fit with both models and returns the best according to the specified criterion
	MM<-fitgene.nl(eset=eset,gene=gene,transforms=transforms,tHVDM=tHVDM,model="MM")
	linscore<-MM$linscore
	fguess<-.exportfree(MM,raw=FALSE)
	Njguess<-4.0
	Djguessmult<-10.0
	Djguess<-fguess[paste(gene,"Dj",sep=".")]
	fguess[paste(gene,"Dj",sep=".")]<-max(Djguess*Djguessmult,0.5)
	fguess[1]<-max(fguess[1],0.1)
	fguess[3]<-min(fguess[3],1e10)
	fguess<-append(fguess,1.1)
	names(fguess)[5]<-paste(gene,"Nj",sep=".")
	hill<-fitgene.nl(eset=eset,gene=gene,tHVDM=tHVDM,transforms=transforms,firstguess=fguess,model="hill")
	cnt=1
	while(hill$score$total-MM$score$total>1e-5){
		message(".")
		if(cnt>20){
			Djguessmult<-(Djguessmult-1.)*0.9+1
		}
		Njguess<-(Njguess-1.)*0.9+1
		fguess[paste(gene,"Dj",sep=".")]<-max(Djguess*Djguessmult,1e-2)
		fguess[paste(gene,"Nj",sep=".")]<-Njguess
		hill<-fitgene.nl(eset=eset,gene=gene,tHVDM=tHVDM,transforms=transforms,firstguess=fguess,model="hill")
		cnt=cnt+1
	}
	if(criterion=="BIC"){
		best="MM"
		if(MM$score$info$BIC>hill$score$info$BIC) best="hill"
	}
	else if(criterion=="AIC"){
		best="MM"
		if(MM$score$info$AIC>hill$score$info$AIC) best="hill"
	}
	else if(criterion=="Q"){
		best="MM"
		if(MM$score$info$Qstat<hill$score$info$Qstat) best="hill"
	}
	else{
		warning("Unknown criterion for model selection (choose among BIC, AIC, Qstat). BIC will be used by default")
		best="MM"
		if(MM$score$info$BIC>hill$score$info$BIC) best="hill"
	}
	if(best=="MM"){
		res<-MM
		res$score$info_hill<-hill$score$info
	}
	else{
		res<-hill
		res$score$info_MM<-MM$score$info
	}
	res$linscore<-linscore
	res
}

.fitgene.nl<-function(eset,gene,tHVDM,transforms,firstguess,model){
	sHVDM<-.initialisescreeningHVDM.nl(eset=eset,gene=gene,tHVDM=tHVDM,transforms=transforms,model=model)
	linscore<-0
	if(missing(firstguess)){
		#estimate parameter values by doing a linear fit first
		lingHVDM<-fitgene.lin(eset=eset,gene=gene,tHVDM=tHVDM,transforms=c("exp"="Dj","exp"="Bj"),firstguess)
		linscore<-lingHVDM$scores$total
		fguess<-.kpfirstguess.nl(sHVDM=sHVDM,linHVDM=lingHVDM,gene=gene)
	}
	else if (class(firstguess)=="numeric" || class(firstguess)=="integer") fguess<-firstguess
	else if (!(is.null(firstguess$type)) && firstguess$type=="indgene"){
		fguess<-.exportfree(HVDM=firstguess,raw=FALSE)
	}
	else{
		lingHVDM<-fitgene.lin(eset=eset,gene=gene,tHVDM=tHVDM,transforms=c("exp"="Dj","exp"="Bj"),firstguess)
		fguess<-.kpfirstguess.nl(sHVDM=sHVDM,linHVDM=lingHVDM,gene=gene)
	}
	sHVDM<-.importfree(HVDM=sHVDM,x=fguess,raw=FALSE)
	results<-.optim(HVDM=sHVDM)
	sHVDM<-.freeparsevaluate(HVDM=sHVDM,x=results$par)
	results$par<-NULL  #this is not needed any more in the HVDM object
	sHVDM$results<-results
	sHVDM$scores<-.scorout(sHVDM)
	sHVDM$fguess<-fguess
	sHVDM$linscore<-linscore
	sHVDM
}

#training step

#a special instance of the training function whereby the input is a linear HVDM training object and used to generate a solution
#to the Michaelis Menten model
#no genemodel option obviously
.training.lintoMM<-function(inputHVDM,transforms=c("exp"="Dj","exp"="Bj","exp"="Kj","expp1"="Nj"),constraints,
																				forcetransforms=TRUE,firstguess){
	if(forcetransforms) tHVDMnl<-.initialisetrainingHVDM.nl(HVDM=inputHVDM,constraints=constraints,
																transforms=transforms)
	else tHVDMnl<-.initialisetrainingHVDM.nl(HVDM=inputHVDM,constraints=constraints,transforms=c())
	#here introduce the possibility to input ad hoc first guesses (input a parameter in training.nl)
	if (missing(firstguess)){
		tHVDMnl<-.tcfirstguess.nl(tHVDM=tHVDMnl,inputHVDM=inputHVDM)
	}
	else{
		bidon="INPUTFIRSTGUESSHERE"
	}
	tHVDMnl$type<-c("training.nl")
	results<-.optim(HVDM=tHVDMnl)
	tHVDMnl<-.freeparsevaluate(HVDM=tHVDMnl,x=results$par)
	results$par<-NULL  #this is not needed any more in the HVDM object
	tHVDMnl$results<-results
	tHVDMnl$scores<-.scorout(tHVDMnl)
	#tHVDMnl$eset<-deparse(substitute(eset))
	tHVDMnl
}

#individual genes

.initialisescreeningHVDM.nl<-function(eset,gene,tHVDM,transforms,model){
	#this only works for the linear model [tbc]
	#the transforms input parameter can be omitted, in which case the transforms applied to the training genes
	#will be set as defaults
	sHVDM<-.createHVDMobj(eset=eset,genenames=gene,genemodels=model,activatornames=tHVDM$par$activators,pdata=tHVDM$pdata)
	sHVDM$pdata<-tHVDM$pdata
	#input known activation profile
	knownnames<-vector("character")
	for (i in 1:length(tHVDM$tc)){
		exprep<-paste(sHVDM$par$activators,sHVDM$tc[[i]]$experiment,sHVDM$tc[[i]]$replicate,sep=".")
		knownnames<-append(knownnames,paste(exprep,sHVDM$tc[[i]]$time,sep="."))
	}
	known<-tHVDM$par$parameters[knownnames]
	sHVDM<-.importknown(HVDM=sHVDM,known=known)
	#import tranformations of input data and start creating the mappings
	if (missing(transforms)) transforms<-tHVDM$distribute$transforms
	transf<-vector("list")
	if(length(transforms)>0){
		fctns<-unique(names(transforms))
		transf<-vector("list",length=length(fctns))
		for(i in 1:length(fctns)){
			transf[[i]]<-vector("list")
			transf[[i]]$pars<-c(outer(gene,transforms[names(transforms)==fctns[i]],FUN=paste,sep="."))
			if (fctns[i]=="exp"){
				transf[[i]]$fctn<-exp
				transf[[i]]$inverse<-log
			}
			else if (fctns[i]=="mexp"){ 
				transf[[i]]$fctn<-.mexp
				transf[[i]]$inverse<-.invmexp
			}
			else if (fctns[i]=="expp1"){ 
				transf[[i]]$fctn<-.expp1
				transf[[i]]$inverse<-.invexpp1
			}

		}
	}
	sHVDM$distribute<-list("known"=known,"transf"=transf,"transforms"=transforms)
	sHVDM<-.cleanupdistributeobj(sHVDM)
	sHVDM
}

#firstguess for individual gene

.kpfirstguess.nl<-function(sHVDM,linHVDM,gene,Kjmult=1.,Djmult=10.){
	#generate the highest values in the linear model
	#names of the activator variables
	ntc<-length(sHVDM$tc)
	tfconcnames<-c()
	actnam<-sHVDM$par$activators
	for (i in c(1:ntc)){
		tfconcnames<-append(tfconcnames,paste(actnam,sHVDM$tc[[i]]$experiment,sHVDM$tc[[i]]$replicate,sHVDM$tc[[i]]$time,sep="."))
	}
	maxact<-max(sHVDM$par$parameters[tfconcnames])
	Kj<-maxact
	sHVDM$par$parameters[paste(gene,"Dj",sep=".")]<-linHVDM$par$parameters[paste(gene,"Dj",sep=".")]*Djmult
	sHVDM$par$parameters[paste(gene,"Bj",sep=".")]<-linHVDM$par$parameters[paste(gene,"Bj",sep=".")]
	#Vj and Kj guess such that Vj/Kj=Sj, first choose Kj big (same for all genes), then Vj accordingly
	sHVDM$par$parameters[paste(gene,"Kj",sep=".")]<-Kj*Kjmult
	sHVDM$par$parameters[paste(gene,"Vj",sep=".")]<-Kj*linHVDM$par$parameters[paste(gene,"Sj",sep=".")]
	#for those genes who have a hill model, simply choose a value of 2.5 for Nj
	if(as.character(sHVDM$par$genemap[gene,"model"])=="hill"){
		sHVDM$par$parameters[paste(gene,"Nj",sep=".")]<-2.5
	}
	res<-.exportfree(sHVDM,raw=FALSE)
	res
}


#initialisation

	#first guess for training set (for nonlinear model)
	
.tcfirstguess.nl<-function(tHVDM,inputHVDM){
	#uses the results of the nonlinear modelling
	#keep the same degradation and basal rates for the genes
	#evaluate Kj and Vj from Sj in the linear model
	
	#generate the highest values in the linear model
	#names of the activator variables
	ntc<-length(tHVDM$tc)
	tfconcnames<-c()
	actnam<-tHVDM$par$activators
	for (i in c(1:ntc)){
		tfconcnames<-append(tfconcnames,paste(actnam,tHVDM$tc[[i]]$experiment,tHVDM$tc[[i]]$replicate,tHVDM$tc[[i]]$time,sep="."))
	}
	if (inputHVDM$type=="training"){
		maxact<-max(inputHVDM$par$parameters[tfconcnames])
		Kj<-maxact
	}
	genes<-rownames(tHVDM$par$genemap)
	for (gene in genes){
		tHVDM$par$parameters[paste(gene,"Dj",sep=".")]<-inputHVDM$par$parameters[paste(gene,"Dj",sep=".")]*10
		tHVDM$par$parameters[paste(gene,"Bj",sep=".")]<-inputHVDM$par$parameters[paste(gene,"Bj",sep=".")]
		#Vj and Kj guess such that Vj/Kj=Sj, first choose Kj big (same for all genes), then Vj accordingly
		if (inputHVDM$type=="training"){
			tHVDM$par$parameters[paste(gene,"Kj",sep=".")]<-Kj
			tHVDM$par$parameters[paste(gene,"Vj",sep=".")]<-Kj*inputHVDM$par$parameters[paste(gene,"Sj",sep=".")]
		}
		else{
			tHVDM$par$parameters[paste(gene,"Kj",sep=".")]<-inputHVDM$par$parameters[paste(gene,"Kj",sep=".")]
			tHVDM$par$parameters[paste(gene,"Vj",sep=".")]<-inputHVDM$par$parameters[paste(gene,"Vj",sep=".")]
			#for those genes who have a hill model, simply choose a value of 1.1 for Nj
			if(as.character(tHVDM$par$genemap[gene,"model"])=="hill"){
				tHVDM$par$parameters[paste(gene,"Nj",sep=".")]<-2.5
			}
		}
	}
	#for the time courses,simply choose the values stored in the linear output, 
	#keeping the constraints that have been specified (and should therefore be stored in the tHVDM object)
	#generate names of the free time points
	n<-length(tHVDM$tc)
	actname<-tHVDM$par$activators
	tfconcnames<-setdiff(tfconcnames,names(tHVDM$distribute$known))
	tHVDM$par$parameters[tfconcnames]<-inputHVDM$par$parameters[tfconcnames]
	namesknown<-names(tHVDM$distribute$known)
	tHVDM$par$parameters[namesknown]<-tHVDM$distribute$known[namesknown]
	tHVDM
}


#structures creation

#training object

.initialisetrainingHVDM.nl<-function(HVDM,constraints,transforms=c(exp="Dj",exp="Bj",exp="Kj",expp1="Nj"),genemodels){
#initialises a nonlinear training object for the nonlinear model, leverages on the linear fit
	genenames<-row.names(HVDM$par$genemap)
	if (missing(genemodels)){
		genemodels<-rep("MM",length(genenames))
		names(genemodels)<-genenames
	}
	pdata<-HVDM$pdata
	tHVDM<-.createHVDMobj(eset=eval(parse(text=HVDM$eset)),genenames=genenames,
						genemodels=genemodels,activatornames=HVDM$par$activators,pdata=pdata)
	tHVDM$pdata<-pdata
	#input of known data
	known<-c()
	knownnames<-c()
	#zero hours time points=known
	for (i in 1:length(tHVDM$tc)){
		known<-append(known,0.0)
		knownnames<-append(knownnames,paste(HVDM$par$activators,tHVDM$tc[[i]]$experiment,tHVDM$tc[[i]]$replicate,0,sep="."))
	}
	#check whether the degradation rate was known and if yes, insert it
	putativename<-paste(genenames[1],"Dj",sep=".")
	if (putativename %in% names(HVDM$distribute$known)){
		known<-append(known,HVDM$distribute$known[putativename])
		knownnames<-append(knownnames,putativename)
	}
	#finally, import the constraints passed to the function
	if (!missing(constraints)){#!!!!!!! the consistency of the input constraints is not checked !!!!
		for(i in c(1:length(constraints))){
			known<-append(known,constraints[i])
			knownnames<-append(knownnames,names(constraints[i]))
		}
	}
	names(known)<-knownnames
	tHVDM<-.importknown(HVDM=tHVDM,known=known)
	#import transformations of input data
	transf<-vector("list")
	if (length(transforms)>0){
		fctns<-unique(names(transforms))
		transf<-vector("list",length=length(fctns))
		for(i in 1:length(fctns)){
			transf[[i]]<-vector("list")
			transf[[i]]$pars<-c(outer(genenames,transforms[names(transforms)==fctns[i]],FUN=paste,sep="."))
			if(fctns[i]=="exp"){
				transf[[i]]$fctn<-exp
				transf[[i]]$inverse<-log
			}
			else if (fctns[i]=="mexp"){
				transf[[i]]$fctn<-.mexp
				transf[[i]]$inverse<-.invmexp
			}
			else if (fctns[i]=="expp1"){
				transf[[i]]$fctn<-.expp1
				transf[[i]]$inverse<-.invexpp1
			}
		}
	}
	tHVDM$distribute<-list("known"=known,"transf"=transf,"transforms"=transforms)
	tHVDM<-.cleanupdistributeobj(tHVDM)
	tHVDM$eset<-HVDM$eset
	tHVDM
}