#set of hidden functions used to compute the measurement error in an eset

.computeinderrtable<-function(X,errtable,refchips){
	#method as described in GB paper
	#prepare the receptacle
	N<-dim(X)[2]
	if(missing(refchips)){
		refchips<-c(1:N)
	}
	m<-dim(errtable)[1]
	res<-matrix(rep(0,m*N),m,N)
	dimnames(res)[[2]]<-dimnames(X)[[2]]
	#compute average
	coeffmul<-rep(1./N,N)
	aver<-as.vector(X %*% coeffmul)
	difsq<-(X-aver)^2
	#determine bounds for classification
	bnds<-c(-20,errtable[,1])
	for(i in c(1:(m-1))){
		bnds[i+1]<-errtable[i,1]/2+errtable[i+1,1]/2
	}
	bnds[m+1]<-log(max(aver)+1)
	if(refchips==0){
		wei<-rep(0,m)
	}
	for(i in c(1:m)){
		keep<-((log(aver)>bnds[i]) & (log(aver)<=bnds[i+1]))
		ngenes<-sum(keep*1)
		coe<-rep(1./ngenes,ngenes)
		res[i,]<-as.vector(t(difsq[keep,]) %*% coe)
		if(refchips==0){
			wei[i]<-ngenes
		}
	}
	if(refchips==0){
		refchips<-.retlessnoisy(X=res,weight=wei,N=3)
		print(refchips)
	}
	#compute average for reference group
	nref<-length(refchips)
	coe<-rep(1./nref,nref)
	avrefchi<-res[,refchips] %*% coe
	ratio<-as.vector(errtable[,2]/avrefchi)
	res<-cbind(errtable[,1],res*ratio)
	dimnames(res)[[2]][1]<-"avg"
	res
}

.retlessnoisy<-function(X,weight,N){
	dimz<-dim(X)
	cweff<-rep(1./dimz[2],dimz[2])
	avgs<-as.vector(X %*% cweff)
	X<-X/avgs
	resis<-as.vector((t(X) %*% weight)/sum(weight))
	names(resis)<-dimnames(X)[[2]]
	resis<-sort(resis)
	res<-names(resis[c(1:N)])
	print(res)
	res
}

.computerrs<-function(eset,errtable,refs){
	#vectorized version
	sanitycheck<-TRUE  #overidden for now, maybe eventulally check that the overall signal corresponds to the
					   #data that was used for calibrating the plattform
	if(sanitycheck){
		X<-exprs(eset)
		sigtab<-.computeinderrtable(X=X,errtable=errtable, refchips=refs)
		sig<-X*0
		N<-dim(X)[1]
		M<-dim(X)[2]
		dimz<-dim(sigtab)
		#sift through error classes
		nclass<-dimz[1]
		message("pre-computations done")
		#do first class first (asolute error constant accross the class)
		mask<-(X<=exp(sigtab[1,1]))*1
		absolu<-(sigtab[1,1+c(1:M)]^0.5)
		sig<-sig+t(t(mask)*absolu)
		message("first class done")
		#loop through the intermediate classes
		for(i in c(1:(nclass-1))){
			mask<-( (X>exp(sigtab[i,1])) & (X<=exp(sigtab[i+1,1])) )*1.
			Ps<-(log(X+1e-10)-sigtab[i,1])/(sigtab[i+1,1]-sigtab[i,1])
			sig <- sig + ((t(t(Ps)*sigtab[i+1,1+c(1:M)]+t(1.-Ps)*sigtab[i,1+c(1:M)]))*mask)^0.5
			message(paste("class",i+1,"/",nclass+1,"done"))
		}
		#do highest class (relative error is constant)
		mask<-(X>exp(sigtab[nclass,1]))*1
		rel<-(sigtab[nclass,1+c(1:M)]^0.5)/exp(sigtab[nclass,1])
		sig<-sig+t(t(mask)*rel)*X
		message("highest class done")
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

