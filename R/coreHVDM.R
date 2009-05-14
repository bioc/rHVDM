#core functions for HVDM, inclusion of more complicated models should require only minimal modifications
#the tag [tbc] indicates where new models have to be added

# ****** model fitting ****************

.optim<-function(HVDM){
	firstguess<-.exportfree(HVDM)
	outnames<-names(HVDM$dm$estimate)
	func<-function(p,HVDMobject){
		HVDMi<-.freeparsevaluate(HVDM=HVDMobject,x=p)
		dm<-HVDMi$dm
		(dm$signal[outnames]-dm$estimate[outnames])/dm$sdev[outnames]
	}
	out<-nls.lm(par=firstguess,fn=func,HVDMobject=HVDM)
	out$diag<-NULL
	out$fvec<-NULL
	out
}


# ****** evaluation *******************

.freeparsevaluate<-function(HVDM,x){
	HVDM<-.importfree(HVDM=HVDM,x=x)
	HVDM<-.evaluate(HVDM)
}

.evaluate<-function(HVDM){
#evaluates model output for parameter values stored in HVDM list
	ntc<-length(HVDM$tc)
	genes<-rownames(HVDM$par$genemap)
	acts<-HVDM$par$activators
	for (gene in genes){
		model<-as.character(HVDM$par$genemap[gene,"model"])
		gparnames<-paste(gene,HVDM$par[model][[1]],sep=".")
		#correct degradation rates that would render the system singular [tbc] maybe...
		degaddress<-paste(gene,"Dj",sep=".")
		x<-HVDM$par$parameters[degaddress]
		HVDM$par$parameters[degaddress]<-min(max(x,1e-10),1e+10)
		gpar<-HVDM$par$parameters[gparnames]
		#print(gpar)
		for(i in 1:ntc){
			tc<-HVDM$tc[[i]]
			exprep<-paste(tc$experiment,tc$replicate,sep=".")
			outnames<-paste(gene,exprep,tc$time,sep=".")
			if(model=="linear"){
				activatorname<-acts[1]
				actprof<-HVDM$par$parameters[paste(activatorname,exprep,tc$time,sep=".")]
				prod<-gpar[1]+gpar[2]*actprof
				mat<-tc$A+diag(gpar[3],length(tc$time))
				HVDM$dm$estimate[outnames]<-solve(mat,prod)
			} #[tbc] other model types will have to be inserted here as "else if"s
			else if(model=="MM"){#MMstands for Michaelis-Menten
				activatorname<-acts[1]
				trfconc<-HVDM$par$parameters[paste(activatorname,exprep,tc$time,sep=".")] #transcription factor concentration
				prod<-gpar[1]+gpar[2]*trfconc/(trfconc+gpar[3]) #Bj+Vj*f/(f+Kj)
				mat<-tc$A+diag(gpar[4],length(tc$time)) #the degradation rate is now the fourth parameter
				HVDM$dm$estimate[outnames]<-solve(mat,prod)
			}
			else if(model=="hill"){#hill is not an accronym
				activatorname<-acts[1]
				trfconc<-HVDM$par$parameters[paste(activatorname,exprep,tc$time,sep=".")] #transcription factor concentration
				
				prod<-.prodhill.patch(trfconc,gpar)
				mat<-tc$A+diag(gpar[5],length(tc$time)) #the degradation rate is now the fourth parameter
				HVDM$dm$estimate[outnames]<-solve(mat,prod)
			} 
		}
	}
	HVDM
}


.prodhill.patch<-function(trfconc,gpar){
	#this patch is added to deal with the cases where the activator is negative (it is set to zero where the input is negative)
	#if the exponent is equal to -1, then no change
	if (gpar[4]>1.0){
		trfconc[trfconc<0]<-0
	}
	#this patch had to be added to deal with the cases where infinity is met in the production term when computing the hill function
	fn<-trfconc^gpar[4]
	#print(fn)
	infp<-is.infinite(fn)
	kn<-gpar[3]^gpar[4]*trfconc^0
	infk<-is.infinite(kn)
	#now, deal with 4 cases
	#none is infinite:normal formula
	#both are infinite: fraction is=0.5
	#kn only is infinite: fraction is zero
	#fn only is infinite: fraction is one
	frac<-fn/(kn+fn)
	frac[infp]<-1.
	frac[infk]<-0.
	frac[infp & infk]<-0.5
	frac[fn==0]<-0.0
	gpar[1]+frac*gpar[2]
}


    #computation of model score distribution

.scorout<-function(HVDM){
	res<-vector("list")
	dm<-HVDM$dm
	devs<-((-dm$estimate+dm$signal)/dm$sdev)
	s2<-devs^2
	#create names
	genes<-rownames(HVDM$par$genemap)
	ntc<-length(HVDM$tc)
	exprep<-rep(0,ntc)
	timevecs<-vector("list",ntc)
	for (i in 1:ntc){
		exprep[i]<-paste(HVDM$tc[[i]]$experiment,HVDM$tc[[i]]$replicate,sep=".")
		timevecs[[i]]<-HVDM$tc[[i]]$time
	}
	#score by gene
	genescore<-rep(0,length(genes))
	names(genescore)<-genes
	for(gene in genes){
		s2names<-vector("character")
		for(i in 1:ntc) s2names<-append(s2names,paste(gene,exprep[i],timevecs[[i]],sep="."))
		genescore[gene]<-sum(s2[s2names])
	}
	#score by timecourse
	tcscore<-rep(0,ntc)
	ntp<-rep(0,ntc)
	names(tcscore)<-exprep
	withintc<-vector("list",ntc)
	for(i in 1:ntc){
		vroot<-paste(exprep[i],timevecs[[i]],sep=".")
		s2names<-outer(genes,vroot,FUN=paste,sep=".")
		tcscore[exprep[i]]<-sum(s2[s2names])
		ntp[i]<-length(timevecs[[i]])
		#score within time courses
		withintc[[i]]<-vector("list")
		withintc[[i]]$name<-exprep[i]
		withintc[[i]]$score<-rep(0,length(timevecs[[i]]))
		names(withintc[[i]]$score)<-as.character(timevecs[[i]])
		for(j in names(withintc[[i]]$score)){
			suffix<-paste(exprep[i],j,sep=".")
			s2names<-paste(genes,suffix,sep=".")
			withintc[[i]]$score[j]<-sum(s2[s2names])
		}
	}
	res$total<-sum(s2)
	res$bygene<-genescore
	res$bytc<-tcscore
	res$bytcpertp<-tcscore/ntp
	res$withintc<-withintc
	res$s2<-s2
	res$devs<-devs
	#add in various indicators: BIC, AIC, Q-stat
	info<-vector("list")
	info$score<-res$total
	info$paracount<-length(HVDM$distribute$free)
	info$datacount<-length(s2)
	datacount<-info$datacount
	paracount<-info$paracount
	info$BIC<-datacount*log(res$total/datacount)+paracount*log(datacount)
	info$AIC<-2*paracount+datacount*(log(-4*asin(-1)*res$total/datacount)+1)
	info$Qstat<-pchisq(res$total,datacount-paracount,lower.tail=FALSE)
	res$info<-info
	res
}



   
   
   #DATA SHIFTING EVALUATIONS ETC.
   
   #data input and output (these two functions can be used with screening object as well)
   
.importfree<-function(HVDM,x,raw=TRUE){
#import free parameters in HVDM data structure
#raw=T means the transform is effected upon import
#the vector x must be properly labelled
	if(raw){
		trans<-HVDM$distribute$transf
		if(length(trans)>0){
			for (i in 1:length(trans)){
				x[trans[[i]]$pars]<-trans[[i]]$fctn(x[trans[[i]]$pars])
			}
		}
	}
	HVDM$par$parameters[HVDM$distribute$free]<-x[HVDM$distribute$free]
	HVDM
}
   
.exportfree<-function(HVDM,raw=TRUE){
#export vector of free parameters (raw values)
#raw=TRUE means the inverse transform is effected
	fnames<-HVDM$distribute$free
	res<-HVDM$par$parameters[fnames]
	#transform using inverse functions
	if (raw){
		trans<-HVDM$distribute$transf
		if(length(trans)>0){
			for (i in 1:length(trans)){
				res[trans[[i]]$pars]<-trans[[i]]$inverse(res[trans[[i]]$pars])
			}
		}
	}
	res
}
   
  
.cleanupdistributeobj<-function(HVDM){
	#list of free parameters
	freeparameterslist<-setdiff(names(HVDM$par$parameters),names(HVDM$distribute$known))
	#remove known parameters from parameters that have to go through a transformation
	ltr<-length(HVDM$distribute$transf)
	if(ltr>0){
		for (i in 1:ltr){
			HVDM$distribute$transf[[i]]$pars<-setdiff(HVDM$distribute$transf[[i]]$pars,names(HVDM$distribute$known))
			HVDM$distribute$transf[[i]]$pars<-intersect(HVDM$distribute$transf[[i]]$pars,names(HVDM$par$parameters))
		}
	}
	HVDM$distribute$free<-freeparameterslist
	HVDM
}	
   
.importknown<-function(HVDM,known){
	#known is a straight vector with known parameter values names must correspond to the names 
	#in the HVDM object. This can also also be used as an import method, e.g. after firstguess has been
	#generated
	namesk<-names(known)
	HVDM$par$parameters[namesk]<-known[namesk]
	HVDM
}
    
    
   #CREATE DATA STRUCTURE 
    
    #create HVDM object
    
.createHVDMobj<-function(eset,genenames,genemodels,activatornames=c("trfac1"),pdata){
#this is fairly general and should not require an update for different typres of models
	ngenes<-length(genenames)
	if (missing(genemodels)) genemodels<-rep("linear",ngenes)  #[tbc] allow several types of models
	res<-vector("list")
	res$tc<-.timecourses(pdata)
	res$dm<-.datamodel(Eset=eset,timecourses=res$tc,genes=genenames)
	res$par<-.paramlist(timecourses=res$tc,genenames=genenames,genemodels=genemodels,activatornames=activatornames)
	res
}


	#create "parameters" object
	
.paramlist<-function(timecourses,genenames,genemodels,activatornames){
	res<-list("linear"=c("Bj","Sj","Dj"),"MM"=c("Bj","Vj","Kj","Dj"),"hill"=c("Bj","Vj","Kj","Nj","Dj"),"activators"=activatornames) # [tbc] add other model types here
	paramvec<-vector("numeric")
	dfgnames<-genenames
	dfgmodel<-genemodels
	nexprep<-length(timecourses)
	ptvec<-1
	#transcription factor(s) activity profile(s)
	for (i in 1:nexprep){
		for (activ in activatornames){
			dfname<-paste(activ,timecourses[[i]]$experiment,timecourses[[i]]$replicate,sep=".")
			actilen<-length(timecourses[[i]]$time)
			toadd<-rep(0.0,actilen)
			names(toadd)<-paste(dfname,timecourses[[i]]$time,sep=".")
			paramvec<-append(paramvec,toadd)
		}
	}
	ptgene<-1
	for (gene in genenames){
		if (genemodels[ptgene]=="linear"){
			paramnames=c("Bj","Sj","Dj")
			paramvals<-rep(1.0,3)
			names(paramvals)<-paste(gene,paramnames,sep=".")
			dfgnames[ptgene]<-gene
			paramvec<-append(paramvec,paramvals)
		}
		else if(genemodels[ptgene]=="MM"){
			paramnames=c("Bj","Vj","Kj","Dj")
			paramvals<-rep(1.0,4)
			names(paramvals)<-paste(gene,paramnames,sep=".")
			dfgnames[ptgene]<-gene
			paramvec<-append(paramvec,paramvals)
		}
		else if(genemodels[ptgene]=="hill"){
			paramnames=c("Bj","Vj","Kj","Nj","Dj")
			paramvals<-rep(1.0,5)
			names(paramvals)<-paste(gene,paramnames,sep=".")
			dfgnames[ptgene]<-gene
			paramvec<-append(paramvec,paramvals)
		}
		#else .... other types of models
		ptgene<-ptgene+1
	}
	res$genemap<-data.frame(model=dfgmodel,row.names=dfgnames)
	res$parameters<-paramvec
	res
}	

	#create "data and model prediction" object (i.e. a list)

.datamodel<-function(Eset,timecourses,genes){
	ngenes<-length(genes)
	nexperepl<-length(timecourses) #this is the time courses count (not length of individual time courses)
	datavec<-vector("numeric")
	sevec<-vector("numeric")
	dfnames<-vector("character")
	for (i in 1:nexperepl){
		for(gene in genes){
			dfname<-paste(gene,timecourses[[i]]$experiment,timecourses[[i]]$replicate,sep=".")
			dfnames<-append(dfnames,dfname)
			colabels<-timecourses[[i]]$explabel
			toadd<-exprs(Eset)[gene,colabels]
			names(toadd)<-paste(dfname,timecourses[[i]]$time,sep=".")
			datavec<-append(datavec,toadd)
			toadd<-assayData(Eset)$se.exprs[gene,colabels]
			names(toadd)<-paste(dfname,timecourses[[i]]$time,sep=".")
			sevec<-append(sevec,toadd)
		}
	}
	modelvec<-datavec*0
	list(signal=datavec,sdev=sevec,estimate=modelvec)
}


    #create time course list from phenodata
    
.timecourses<-function(pdata){
	namesoftc<-.listoftimecourses(pdata)
	ntcs<-length(namesoftc)
	for (i in 1:ntcs){
		bidon<-.timevectorandexplabels(pdata,namesoftc,i)
		namesoftc[[i]]$time<-bidon[[1]]
		namesoftc[[i]]$explabel<-bidon[[2]]
		namesoftc[[i]]$A<-.diffopw(namesoftc[[i]]$time)   #[tbc] allow more generalised bounds
	}
	namesoftc
}

.timevectorandexplabels<-function(pdata,tclist,i){
#returns time vector and corresponding experiment labels, ordered
	repl<-tclist[[i]]$replicate
	exper<-tclist[[i]]$experiment
	partpd<-pdata[pdata$experiment==exper & pdata$replicate==repl,]
	exlab<-rownames(partpd)
	time<-partpd[exlab,"time"]
	orderindex<-order(time)
	list(time[orderindex],exlab[orderindex])
}

.listoftimecourses<-function(pdata){
	res<-vector("list")
	bidon<-pdata$experiment
	listoftc<-c()
	for (i in 1:length(bidon)){
		experiment<-pdata[i,"experiment"]
		replicate<-pdata[i,"replicate"]
		thistc<-paste(experiment,replicate)
		if (!(.isthisinto(thistc,listoftc))){
			listoftc<-append(listoftc,thistc)
			tcnumb<-length(listoftc)
			res[[tcnumb]]<-vector("list")
			res[[tcnumb]]$experiment<-experiment
			res[[tcnumb]]$replicate<-replicate
		}
	}
	res[order(listoftc)]
}

.isthisinto<-function(this,intothat){
	if(length(intothat)==0) res<-FALSE
	else res<-sum(intothat==this)
	res
}


# ****** differential operator ********

.diffop<-function(timevec,interprange,slopeatzero=TRUE){
	#testing that the timevector is appropriately formatted
	if (timevec[1]!=0.0) {
		print("the first time point is not equal to zero, exiting");
		break;
	}
	else{	
		n<-length(timevec)
		res<-array(rep(0,n*n),dim=c(n,n))
		for(i in 1:n){
			p=interprange[i,1];
			r=interprange[i,2];
			for (k in 1:n){
				if( (k>r) || (k<p) || ((slopeatzero) && (i==1)) ){
					res[i,k]=0.0;
				}
				else{
					if(k==i){
						for(j in p:r){
							if(j!=i) res[i,k]=res[i,k]+1.0/(timevec[i]-timevec[j])
						}
						if( (p==1) && (slopeatzero) ) res[i,k]=res[i,k]+1/timevec[i];
					}
					else{
						res[i,k]=1/(timevec[k]-timevec[i])
						for(j in p:r){
							if( (j!=k) && (j!=i) ) res[i,k]=res[i,k]*(timevec[i]-timevec[j])/(timevec[k]-timevec[j]);
						}
						if( (p==1) && (slopeatzero) ){ 
							if(k!=1) res[i,k]=res[i,k]*timevec[i]/timevec[k]
							else{
								s<-0;
								for (j in 2:r) s<-s-1/timevec[j];
								s<-1-s*timevec[i];
								res[i,k]<-res[i,k]*s;
							}
						}
					}
				}
			}
		}
	}
	res;
}

.diffopw<-function(timevec,slopeatzero=TRUE){
	#wrapper for the .diffop function with default parameters
	#use 5 interpolation points where possible
	n<-length(timevec);
	ranges<-array(rep(0,2*n),dim=c(n,2))
	for(i in 1:n){
		ranges[i,1]<-max(i-2,1);
		ranges[i,2]<-min(i+2,n);
	}
	.diffop(timevec,interprange=ranges,slopeatzero=slopeatzero);
}

#useful functions

.mexp<-function(x){
	sm0<-(x<0)
	res<-x
	res[sm0]<-exp(x[sm0])
	res[!sm0]<-x[!sm0]+1
	res
}

.invmexp<-function(x){
	gt1<-(x>1)
	res<-x
	res[gt1]<-x[gt1]-1
	res[!gt1]<-log(x[!gt1])
	res
}

#The function below (expp1, with inverse invexpp1), makes sure its output always exceeds 1

.expp1<-function(x){
	res<-exp(x)+1
	res
}

.invexpp1<-function(x){
	leqq1<-(x<=1)
	res<-x
	res[leqq1]<-1e-15
	res[!leqq1]<-log(x[!leqq1]-1.0)
	res
}

#check time vectors for HVDMcheck method

.checktimevectors<-function(pdata){
	res<-TRUE
	message("*** checking time vectors:")
	if(!((class(pdata$time)=="numeric") | (class(pdata$time)=="integer"))){
		message("   B) the time field in pdata is not numeric")
		res<-FALSE
	}
	#inventory of the number of time courses
	n<-length(pdata$time)
	exper<-vector("numeric")
	repli<-vector("numeric")
	expre<-vector("numeric")
	for(i in 1:n){
		thise<-as.character(pdata[i,"experiment"])
		thisr<-as.character(pdata[i,"replicate"])
		thisone<-paste(thise,thisr,sep=".")
		if (!(thisone %in% expre)){
			expre<-append(expre,thisone)
			exper<-append(exper,thise)
			repli<-append(repli,thisr)
		}
	}
	m<-length(expre)
	for	(i in 1:m){
		ttime<-pdata[(pdata$experiment==exper[i] & pdata$replicate==repli[i]),"time"]
		#check there are no repeats
		if (length(unique(ttime))<length(ttime)){
			message(paste("   B) some time values in",expre[i],"seem to be duplicated"))
			res<-FALSE
		}
		#check there is a zero time point
		if(sum(ttime==0)<1){
			message(paste("   B) there is no zero time point in",expre[i]))
			res<-FALSE
		}
		#check all elements are nonnegative
		if(sum(ttime<0)>0){
			message(paste("   B) some time values in",expre[i],"are negative"))
			res<-FALSE
		}
	}
	if(res) message("					OK")
	res
}