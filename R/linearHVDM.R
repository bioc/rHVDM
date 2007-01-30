#functions for individual gene in linear situation (ie linear model) [tbc]
#these functions will have to be rewritten completely for other situations

.screenall<-function(eset,genes,tHVDM,transforms=c("exp"="Dj","exp"="Bj")){
	n<-length(genes)
	res<-data.frame(row.names=genes,
					model_score=rep(0,n),
					lm_message=rep(0,n),
					basal=rep(0,n),
					sensitivity=rep(0,n),
					degradation=rep(0,n),
					sens_z_score=rep(0,n),
					hess_rank=rep(0,n),
					basal_d=rep(0,n),
					sensitivity_d=rep(0,n),
					degradation_d=rep(0,n),
					basal_u=rep(0,n),
					sensitivity_u=rep(0,n),
					degradation_u=rep(0,n)
					)
	for(gene in genes){
		res[gene,]<-.genescreen(eset=eset,gene=gene,tHVDM=tHVDM,template=res[gene,],transforms=transforms)
	}
	res
}

.genescreen<-function(eset,gene,tHVDM,template,transforms,firstguess){
	#wrapper function for screening individual genes [tbc] works only for linear model
	res<-template
	sHVDM<-.fitgene(eset=eset,gene=gene,tHVDM=tHVDM,transforms=transforms,firstguess=firstguess)
	res[gene,]$model_score<-sHVDM$scores$total
	res[gene,]$lm_message<-sHVDM$results$message
	hout<-.hessstuff(sHVDM)
	centfit<-sHVDM$par$parameters[sHVDM$distribute$free]
	res[gene,]$basal<-centfit[1]
	res[gene,]$sensitivity<-centfit[2]
	res[gene,]$degradation<-centfit[3]
	if(class(hout)=="list"){
		res[gene,]$sens_z_score<-(centfit[2]*2*1.96)/(hout[[2]][2]-hout[[1]][2])
		res[gene,]$hess_rank<-3
		res[gene,]$basal_d<-hout[[1]][1]
		res[gene,]$sensitivity_d<-hout[[1]][2]
		res[gene,]$degradation_d<-hout[[1]][3]
		res[gene,]$basal_u<-hout[[2]][1]
		res[gene,]$sensitivity_u<-hout[[2]][2]
		res[gene,]$degradation_u<-hout[[2]][3]
	}
	else{
		res[gene,]$sens_z_score<--1
		res[gene,]$hess_rank<-hout
		res[gene,]$basal_d<-centfit[1]
		res[gene,]$sensitivity_d<-centfit[2]
		res[gene,]$degradation_d<-centfit[3]
		res[gene,]$basal_u<-centfit[1]
		res[gene,]$sensitivity_u<-centfit[2]
		res[gene,]$degradation_u<-centfit[3]
	}
	res
}

.hessstuff<-function(sHVDM){#NB works only for 3x3 matrix, ie linear model [tbc]
	#returns the individual (raw) sigmas and correlations/ zero if hessian is singular
	H<-sHVDM$results$hessian
	hr<-.hessrank(H)
	if(hr<3) res<-hr
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

.fitgene<-function(eset,gene,tHVDM,transforms,firstguess){
	sHVDM<-.initialisescreeningHVDM(eset=eset,gene=gene,tHVDM=tHVDM,transforms=transforms)
	if (missing(firstguess)) fguess<-.kpfirstguess(sHVDM,gene)
	else if (class(firstguess)=="numeric" || class(firstguess)=="integer") fguess<-firstguess
	else if (!(is.null(firstguess$type)) && firstguess$type=="indgene"){
		fguess<-.exportfree(HVDM=firstguess,raw=FALSE)
	}
	else fguess<-.kpfirstguess(sHVDM,gene)
	names(fguess)<-paste(gene,c("Bj","Sj","Dj"),sep=".")
	sHVDM<-.importfree(HVDM=sHVDM,x=fguess,raw=FALSE)
	results<-.optim(HVDM=sHVDM)
	sHVDM<-.freeparsevaluate(HVDM=sHVDM,x=results$par)
	results$par<-NULL  #this is not needed any more in the HVDM object
	sHVDM$results<-results
	sHVDM$scores<-.scorout(sHVDM)
	names(fguess)<-c("Bj","Sj","Dj")
	sHVDM$fguess<-fguess
	sHVDM
}

# ****** initialisation ***************
   
   #firstguess for training set (N.B. this is specific for the training set and for the linear model)
   
.tcfirstguess<-function(tHVDM){
	genes<-rownames(tHVDM$par$genemap)
	agene<-genes[1] #for this gene, sensitivity and degradation rate are known
	ntc<-length(tHVDM$tc)
	#determine basal rate for this gene
	tp0<-vector("character",ntc)
	for(i in 1:ntc) tp0[i]<-paste(tHVDM$tc[[i]]$experiment,tHVDM$tc[[i]]$replicate,"0",sep=".")
	tsnam<-paste(agene,tp0,sep=".")
	w<-1/tHVDM$dm$sdev[tsnam]
	w<-w/sum(w)
	Dj<-tHVDM$par$parameters[paste(agene,"Dj",sep=".")]
	Bj<-sum(Dj*tHVDM$dm$signal[tsnam]*w)
	tHVDM$par$parameters[paste(agene,"Bj",sep=".")]<-Bj
	#determine transcription factor profiles using anchoring gene only
	actpro<-tHVDM$par$activators
	for(i in 1:ntc){
		exprept<-paste(tHVDM$tc[[i]]$experiment,tHVDM$tc[[i]]$replicate,tHVDM$tc[[i]]$time,sep=".")
		outnames<-paste(actpro,exprept,sep=".")
		x<-tHVDM$dm$signal[paste(agene,exprept,sep=".")]
		f<-tHVDM$tc[[i]]$A %*% x + Dj*x -Bj
		tHVDM$par$parameters[outnames]<-f
	}
	known<-tHVDM$distribute$known
	tHVDM$par$parameters[names(known)]<-known[names(known)]
	for (gene in genes[2:length(genes)]){
		tHVDM$par$parameters[paste(gene,tHVDM$par$linear,sep=".")]<-.kpfirstguess(HVDM=tHVDM,gene=gene)
	}
	tHVDM
}

	#firstguess for kinetic parameters of an individual gene (specific for normal gene should work in
	#screening as well as in training situation)
	#returns a vector containing basal,sensitivity,degradation

.kpfirstguess<-function(HVDM,gene){
	w<-vector("numeric")
	p<-vector("numeric")
	actpro<-HVDM$par$activators
	ntc<-length(HVDM$tc)
	for (i in 1:ntc){
		tc<-HVDM$tc[[i]]
		exprept<-paste(tc$experiment,tc$replicate,tc$time,sep=".")
		xj<-HVDM$dm$signal[paste(gene,exprept,sep=".")]
		y<-tc$A %*% xj
		f<-HVDM$par$parameters[paste(actpro,exprept,sep=".")]
		cons<-rep(1.0,length(xj))
		Xp<-rbind(cons,f,-xj)
		XpX<-Xp %*% t(Xp)
		Xpy<-Xp %*% y
		p<-cbind(p,solve(XpX,Xpy))
		w<-append(w,sum(1/HVDM$dm$sdev[paste(gene,exprept,sep=".")]))
	}
	w<-w/sum(w)
	res<-c(p %*% w)
	#do a sanity check for those parameters that go through a transform
	trans<-HVDM$distribute$transforms
	if(.isthisinto("Bj",trans) && res[1]<=0.0) res[1]<-1.0
	if(.isthisinto("Sj",trans) && res[2]<=0) res[2]<-1.0
	if(.isthisinto("Dj",trans) && res[3]<=0) res[3]<-1.0
	res
}

   #structures creation
   
.initialisescreeningHVDM<-function(eset,gene,tHVDM,transforms){
	#this only works for the linear model [tbc]
	#the transforms input parameter can be omitted, in which case the transforms applied to the training genes
	#will be set as defaults
	sHVDM<-.createHVDMobj(eset=eset,genenames=gene,activatornames=tHVDM$par$activators,pdata=tHVDM$pdata)
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
		}
	}
	sHVDM$distribute<-list("known"=known,"transf"=transf,"transforms"=transforms)
	sHVDM<-.cleanupdistributeobj(sHVDM)
	sHVDM
}
   
   
.initialisetrainingHVDM<-function(eset,trainingset,degrate,pdata,trfname="trfact1",transforms){
	#this only works for the linear model (GB2006) [tbc]
	#anchoring gene is the one in the first position of the training set list
	if (missing(pdata)) pdata<-pData(eset)
	tHVDM<-.createHVDMobj(eset=eset,genenames=trainingset,activatornames=c(trfname),pdata=pdata)
	tHVDM$pdata<-pdata
	#input known degradation rate and also sensitivity=1.0 for the same gene
	bidon<-paste(trainingset[1],"Sj",sep=".")
	known<-1.0
	knownnames<-bidon
	if(!missing(degrate)){
		bidon<-paste(trainingset[1],"Dj",sep=".")
		known<-append(known,degrate)
		knownnames<-append(knownnames,bidon)
	}
	#input known zero hours time points for activation profile
	for (i in 1:length(tHVDM$tc)){
		known<-append(known,0.0)
		knownnames<-append(knownnames,paste(trfname,tHVDM$tc[[i]]$experiment,tHVDM$tc[[i]]$replicate,0,sep="."))
	}
	names(known)<-knownnames
	tHVDM<-.importknown(HVDM=tHVDM,known=known)
	#import trandformations of input data and start creating the mappings
	transf<-vector("list")
	if(length(transforms)>0){
		fctns<-unique(names(transforms))
		transf<-vector("list",length=length(fctns))
		for(i in 1:length(fctns)){
			transf[[i]]<-vector("list")
			transf[[i]]$pars<-c(outer(trainingset,transforms[names(transforms)==fctns[i]],FUN=paste,sep="."))
			if (fctns[i]=="exp"){
				transf[[i]]$fctn<-exp
				transf[[i]]$inverse<-log
			}
			else if (fctns[i]=="mexp"){ 
				transf[[i]]$fctn<-.mexp
				transf[[i]]$inverse<-.invmexp
			}
		}
	}
	tHVDM$distribute<-list("known"=known,"transf"=transf,"transforms"=transforms)
	tHVDM<-.cleanupdistributeobj(tHVDM)
	tHVDM
}