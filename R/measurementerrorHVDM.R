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

.returnplattidnumber<-function(plattid,repository=.rHVDMplattforms){#index of plattid returned
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
	n<-length(.rHVDMplattforms)
	if (n==0){
		message("no microarray plattforms are supported in rHVDM yet")
		return()
	}
	else{
		plattform_id<-c(1:n)
		description<-c(1:n)
		for(i in c(1:n)){
			plattform_id[i]<-.rHVDMplattforms[[i]]$id
			description[i]<-.rHVDMplattforms[[i]]$description
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
#this last member will be left out from (1.9.5, march 2009)

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

#data for the plattform tables

.rHVDMplattforms<-list()

#affy U133

affy133<-list()
affy133$id<-"affy_HGU133A"
affy133$description<-"Affymetrix Human genome 133A 3-prime array, MAS5 summary"
affy133$source<-"Affymetrix spike-in experiment, data accessible on the Affimetrix website"
affy133$errtable<-t(matrix(c(
-0.2106170,     3.030952,
0.4637223  ,   4.488031,
0.8341946  ,   6.786352,
1.1527715  ,   8.866159,
1.4474100  ,  13.244326,
1.7220945  ,  15.747942,
1.9901272  ,  17.893599,
2.2386661  ,  22.929129,
2.4586614  ,  24.337833,
2.6499757  ,  27.111229,
2.8247291  ,  29.096279,
2.9791003  ,  31.097968,
3.1231124  ,  34.673634,
3.2639854  ,  40.087314,
3.4078585  ,  41.540187,
3.5490655  ,  48.236199,
3.6917189  ,  55.368312,
3.8422731  ,  59.359361,
3.9964226  ,  66.993106,
4.1579789  ,  75.847599,
4.3225134  , 103.653371,
4.4947673  , 119.709417,
4.6758389  , 148.446633,
4.8570693  , 170.795388,
5.0548622  , 227.148177,
5.2864374  ,313.226515,
5.5693897  , 468.947419,
5.9312727  , 770.988764,
6.4451406 , 1836.468715,
7.6395270, 13000.404178
),2,30))

.rHVDMplattforms[[1]]<-affy133
rm(affy133)

#affy gene expression array

affyge<-list()
affyge$id<-"affy_ID10ST"
affyge$description<-"Affy Human genome arrays (1.0 ST), summarized using PLIER with quantile normalization"
affyge$source<-"computed on technical tripicates of a T-cell cell line (MOLT4), at the Institute of Child Health, London, UK. m.hubank@ich.ucl.ac.uk"
affyge$errtable<-t(matrix(c(
0.3966758 ,    7.357677,
1.5387863 ,    7.390087,
2.0228899 ,    9.665248,
2.3300754 ,    9.426613,
2.5602427 ,   11.262398,
2.7582411 ,   14.556250,
2.9275274 ,   16.102466,
3.0902940 ,   18.270034,
3.2452735 ,   23.117356,
3.3987511 ,   28.259277,
3.5488812 ,   33.033751,
3.6961924 ,   39.944499,
3.8452394 ,   46.567401,
3.9972747 ,   55.044254,
4.1510491 ,   68.557068,
4.3093431 ,   83.363928,
4.4709532 ,  104.748336,
4.6328518 ,  127.738549,
4.7909148 ,  145.132212,
4.9529966 ,  183.702869,
5.1177288 ,  226.670256,
5.2795822 ,  277.571413,
5.4389206 ,  339.290293,
5.6017597 ,  387.109394,
5.7791409 ,  550.898651,
5.9753123 ,  849.800103,
6.2182958 , 1126.856340,
6.5159188 , 1812.932056,
6.9818574 , 4645.247199,
7.9740574 ,27974.882980
),2,30))

.rHVDMplattforms[[2]]<-affyge
rm(affyge)
