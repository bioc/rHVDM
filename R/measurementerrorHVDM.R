#set of hidden functions used to compute the measurement error in an eset

.returnsigma<-function(x,tabl){
#looks up value log(x) in a table and return interpolated value
#the table should have two colums where the first is log transformed of the signal
#and the second the square of the variance of the measurement error
#the function returns the standard devation
	if(x<=0){
		lx<-(-15)
	}
	else{
		lx<-log(x)
	}
	n<-dim(tabl)[1]
	if(lx<tabl[1,1]){
		res<-tabl[1,2]^0.5
	}
	else{
		if(lx>=tabl[n,1]){
			rel<-tabl[n,2]^0.5/exp(tabl[n,1])
			res<-rel*x
		}
		else{#do the normal search with interpolation and all that
			lowbnd<-1
			upebnd<-n
			while(upebnd-lowbnd>1){
				tesbnd<-(upebnd+lowbnd)/2
				if(lx<tabl[tesbnd,1]){
					upebnd<-tesbnd
				}
				else{
					lowbnd<-tesbnd
				}
			}
			#the difference between upper and lower bound is now one, so have to interpolate a value
			totrendj<-tabl[upebnd,1]-tabl[lowbnd,1]
			parrendj<-lx-tabl[lowbnd,1]
			p<-parrendj/totrendj
			res<-((1-p)*tabl[lowbnd,2]+p*tabl[upebnd,2])^0.5
			res
		}
	}
}

.returnmulfac<-function(x,tabl){
#looks up value log(x) in a table and return interpolated value
#the table should have two colums where the first is log transformed of the signal
#and the second the square of the variance of the measurement error
#the function returns the standard devation
	if(x<=0){
		lx<-(-15)
	}
	else{
		lx<-log(x)
	}
	n<-dim(tabl)[1]
	if(lx<tabl[1,1]){
		res<-tabl[1,2]
	}
	else{
		if(lx>=tabl[n,1]){
			rel<-tabl[n,2]
			res<-rel
		}
		else{#do the normal search with interpolation and all that
			lowbnd<-1
			upebnd<-n
			while(upebnd-lowbnd>1){
				tesbnd<-(upebnd+lowbnd)/2
				if(lx<tabl[tesbnd,1]){
					upebnd<-tesbnd
				}
				else{
					lowbnd<-tesbnd
				}
			}
			#the difference between upper and lower bound is now one, so have to interpolate a value
			totrendj<-tabl[upebnd,1]-tabl[lowbnd,1]
			parrendj<-lx-tabl[lowbnd,1]
			p<-parrendj/totrendj
			res<-((1-p)*tabl[lowbnd,2]+p*tabl[upebnd,2])
			res
		}
	}
}


.calculatecorrectionfactors<-function(X,errbnd,refs){
	bnds<-c(-20,errbnd)
	n<-length(errbnd)
	for(i in c(1:(n-1))){
		bnds[i+1]<-errbnd[i]/2+errbnd[i+1]/2
	}
	bnds[n+1]<-log(max(X)+1)
	#compute average per observation
	dims<-dim(X)
	mulavg<-rep(1/dims[2],dims[2])
	avre<-as.vector(X %*% mulavg)
	deviates<-matrix(rep(0,dims[2]*n),n,dims[2])
	avgates<-0*c(1:n)
	Xmju2<-((X-avre)^2)^0.5
	for(i in c(1:n)){
		sel<-((log(avre)>bnds[i]) & (log(avre)<=bnds[i+1]))
		sub<-t(Xmju2[sel,])
		dimz<-dim(sub)
		ultz<-rep(1/dimz[2],dimz[2])
		deviates[i,]<-as.vector(sub %*% ultz)
		avgates[i]<-sum(avre[sel]*ultz)
	}
	print(avgates)
	print(deviates)
	deviates<-cbind(avgates,deviates)
	rnams<-paste("r",c(1:n))
	dimnames(deviates)<-list(rnams,c("clcenterz",dimnames(X)[[2]]))
	#finally, normalise to array references
	refarr<-deviates[,refs]
	dims<-dim(refarr)
	mulref<-rep(1/dims[2],dims[2])
	normaliz<-as.vector(refarr %*% mulref)
	res<-deviates/normaliz
	res[,1]<-log(avgates)
	print(res)
	res
}

.computerrs<-function(eset,errtable,refs){
	message("hi there")
	sanitycheck<-TRUE  #overidden for now, maybe eventulally check that the overall signal corresponds to the
					   #data that was used for calibrating the plattform
	if(sanitycheck){
		X<-exprs(eset)
		mftab<-.calculatecorrectionfactors(X=X,errbnd=errtable[,1], refs=refs)
		sig<-X
		dnams<-dimnames(X)
		for(aarray in dnams[[2]]){
			for(gene in dnams[[1]]){
				res<-.returnsigma(X[gene,aarray],errtable)
				multfac<-.returnmulfac(X[gene,aarray],mftab[,c("clcenterz",aarray)])
				#message(paste(aarray,gene,X[gene,aarray],res,multfac))
				sig[gene,aarray]<-res*multfac
			}
			message(paste(aarray,"done"))
		}
		res<-assayDataElementReplace(eset,"se.exprs",sig)
	}
	else{
		message("sanity check not passed, errors not computed")
		res<-eset
	}
	res
}

#stucture of .plattform list:
#indexed list under each index, there is a list whose constituents are in text format unless specified
#.plattform[[i]]$id: plattform identification
#.plattform[[i]]$description: short description
#.plattform[[i]]$source: source of calibration data
#.plattform[[i]]$errtable: nx2 matrix containing the log of the signal in the first column and the variance in the second
#.plattform[[i]]$geneids: vector of gene ids (in text format) that could be used to identify the plattform

.plattformmatch<-function(eset){#must return the ***index*** of the plattform or a "notfound" string
	res<-"notfound"
	res
}

.returnplattidnumber<-function(plattid,repository=rHVDMplattforms){#index of plattid returned
	res<-"notfound"
	n<-length(repository)
	for (i in c(1:n)){
		if(repository[[i]]$id==plattid){
			res<-i
		}
	}
	res
}

.returnsupportedplattforms<-function(){
	n<-length(rHVDMplattforms)
	if (n==0){
		message("no microarray plattforms are supported in rHVDM yet")
		return()
	}
	else{
		plattform_id<-c(1:n)
		description<-c(1:n)
		for(i in c(1:n)){
			plattform_id[i]<-rHVDMplattforms[[i]]$id
			description[i]<-rHVDMplattforms[[i]]$description
		}
		res<-cbind(plattform_id,description)
		res
	}
}

#some utilities that will be helpful not directly in rHVDM but to compute and populate the known plattforms repository
#add suffix .ut to make things less confusing (hopefully)
#calculation of the table

.fillplattrepository.ut<-function(repository,X,id,description,source,overwrite=FALSE){
	#all inputs are compulsory
	#to rewrite an existing entry, change the overwrite switch to TRUE
	#checking that everything has been input (X and repository are not verified)
	if(missing(id)){
		message("plattorform *id* has to be given")
		res<-repository
	}
	else{
		if (missing(description)){
			message("plattform *description* has to be given")
			res<-repository
		}
		else{
			if (missing(source)){
				message("some indication of where the data comes from (the *source*) has to be given")
				res<-repository
			}
			else{
				#check that the id is not already used
				newdex<-.returnplattidnumber(plattid=id,repository=repository)
				if (newdex=="notfound"){
					newdex<-length(repository+1)
					res<-.repositoryupdate(repository=repository,X=X,id=id,description=description,source=source,slot=newdex)
				}
				else{
					message(paste("plattform",id,"allready in repository"))
					if(overwrite){
						message("overwriting the old entry for that id")
						res<-.repositoryupdate(repository=repository,X=X,id=id,description=description,source=source,slot=newdex)
					}
					else{
						message("this entry cannot be overwritten, change the 'overwrite' switch")
						res<-repository
					}
				}
			}
		}
	}
	res
}

#.plattform[[i]]$id: plattform identification
#.plattform[[i]]$description: short description
#.plattform[[i]]$source: source of calibration data
#.plattform[[i]]$errtable: nx2 matrix containing the log of the signal in the first column and the variance in the second
#.plattform[[i]]$geneids: vector of gene ids (in text format) that could be used to identify the plattform

.repositoryupdate<-function(repository,X,id,description,source,slot){
	resplatt<-list()
	resplatt$id<-id
	resplatt$description<-description
	resplatt$source<-source
	errtable<-.bin.ut(X)
	resplatt$errtable<-errtable
	geneids<-dimnames(X)[[1]]
	resplatt$geneids<-geneids
	res<-repository
	res[[slot]]<-resplatt
	res
}

.calcvar.ut<-function(X){
	#returns a vector with the variances
	dimz<-dim(X)
	X2<-X^2
	sumv<-rep(1/dimz[2],dimz[2])
	mju<-(X %*% sumv)
	mju2<-(X2 %*% sumv)
	var<-dimz[2]*(mju2-mju^2)/(dimz[2]-1)
	var
}

.bin.ut<-function(X,ngroup=30){
	probbs<-c(0,c(1:ngroup))/ngroup
	qs<-quantile(X,probs=probbs)
	qs[ngroup+1]<-qs[ngroup+1]+1
	obsvar<-as.vector(.calcvar.ut(X))
	res<-matrix(rep(0,ngroup*2),ngroup,2)
	for(i in c(1:ngroup)){
		boolettes<-((X>=qs[i]) & (X<qs[i+1]))
		nobs<-sum(boolettes*1)
		meann<-sum((boolettes*1)*X)/nobs
		varr<-sum((boolettes*1)*obsvar)/nobs
		res[i,1]<-log(meann)
		res[i,2]<-varr
	}
	res
}

.binb.ut<-function(X,quant,ngroup=30){
	probbs<-c(0,c(1:ngroup))/ngroup
	obsvar<-as.vector(.calcvar.ut(X))
	obsq<-obsvar*0
	for(i in length(obsq)){
		obsq[[i]]<-quantile(X[i,],prob=quant)
	}
	qs<-quantile(obsq,probs=probbs)
	qs[ngroup+1]<-qs[ngroup+1]+1
	res<-matrix(rep(0,ngroup*2),ngroup,2)
	for(i in c(1:ngroup)){
		boolettes<-((obsq>=qs[i]) & (obsq<qs[i+1]))
		nobs<-sum(boolettes*1)
		meann<-sum(obsq[boolettes])/nobs
		varr<-sum(obsvar[boolettes])/nobs
		res[i,1]<-log(meann)
		res[i,2]<-varr
	}
	res
}

