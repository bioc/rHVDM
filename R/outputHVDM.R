# ****** output functions *************

.SUFFX.rHVDM<-"_files"

.sreport<-function(lHVDM,name){
	
	#file name and directories management
	if (missing(name)){
		rname<-as.character(floor(runif(1,100000,999999)))
		prefix<-paste("HVDM_",rname,"_screening",sep="")
	}
	else{
		prefix<-name
	}
	dirname<-paste(prefix,.SUFFX.rHVDM,sep="")
	dir.create(dirname)
	HTMLfile<-paste(prefix,".HTML",sep="")
	
	#header section
	tHVDM<-lHVDM$tHVDM
	HTML(x="<a name=Top></a>",file=HTMLfile,append=FALSE)
	HTML.title("Report for HVDM screening() step",HR=1,file=HTMLfile)
	
	#links
	HTML("<a href=#Input>1) Inputs-training</a>  <a href=#Inputsc>2) Inputs-screening</a>   <a href=#results>3) Results</a>  <a href=#class1>4) Putative targets</a>",file=HTMLfile)
	
	HTML("<hr  width=100% size=20 noshade>",file=HTMLfile)
	
	HTML(x="<a name=Input></a>",file=HTMLfile)
	HTML.title(x="1) Parameters input to the training() command",HR=2,file=HTMLfile)
	HTML("<a href=#Top>Top</a>  <a href=#Inputsc>2) Inputs-screening</a>   <a href=#results>3) Results</a>  <a href=#class1>4) Putative targets</a>",file=HTMLfile)
	#name of expression set used
	HTML.title(x=paste("Expression set:",tHVDM$eset),HR=4,file=HTMLfile)
	
	#name of transcription factor
	def<-""
	if(tHVDM$par$activators=="trfact1") def<-"(default name)"
	HTML.title(paste("Transcription factor name:",tHVDM$par$activators,def),HR=4,file=HTMLfile)
	
	#transformations used for the fitting
	if(length(tHVDM$distribute$transf)==0) HTML.title("Transformations applied for the fitting: None",HR=4,file=HTMLfile)
	else{
		HTML.title("Transformations applied for the fitting (training):",HR=4,file=HTMLfile)
		HTML(tHVDM$distribute$transforms,file=HTMLfile)
		HTML(names(tHVDM$distribute$transforms),file=HTMLfile)
	}
	
	#set of training genes
	HTML.title(x="Training genes:",HR=4,file=HTMLfile)
	HTML(x=rownames(tHVDM$par$genemap),file=HTMLfile,align="left")
	
	#name of anchoring gene
	HTML.title(x="Anchoring gene:",HR=4,file=HTMLfile)
	anchgene<-rownames(tHVDM$par$genemap)[1]
	testname<-paste(anchgene,"Dj",sep=".")
	if (testname %in% tHVDM$distribute$free){
		HTML(x=paste(anchgene,": not anchored"),file=HTMLfile)
		HTML(x="<i> degradation rate has not been given, try parameter degrate in function training() </i>",file=HTMLfile)
	}
	else HTML(x=paste(anchgene,", degradation rate: ",tHVDM$par$parameters[testname]," (unit time)^(-1)",sep=""),file=HTMLfile)

	#Pheno data
	HTML.title(x="Pheno data:",HR=4,file=HTMLfile)
	HTML(x=tHVDM$pdata,file=HTMLfile,align="left")
	HTML("<a href=#Top>Top</a>",file=HTMLfile)
	HTML("<hr  width=100% size=20 noshade>",file=HTMLfile)

	HTML(x="<a name=Inputsc></a>",file=HTMLfile)
	HTML.title(x="1) Parameters input to the screening() command",HR=2,file=HTMLfile)
	
	HTML("<a href=#Top>2) Top</a>  <a href=#Input>1) Inputs-training</a>  <a href=#results>3) Results</a>  <a href=#class1>4) Putative targets</a>",file=HTMLfile)

	#classification bounds
	HTML.title("Classification Bounds",HR=3,file=HTMLfile)
	bounds<-lHVDM$class1bounds
	HTML.title("Sensitivity Z-score (lower bound)",HR=4,file=HTMLfile)
	HTML.title(bounds$zscore,HR=4,file=HTMLfile)
	HTML.title("Model Score (upper bound)",HR=4,file=HTMLfile)
	HTML.title(bounds$modelscore,HR=4,file=HTMLfile)
	HTML.title("degradation rate (range), in (unit time)^(-1)",HR=4,file=HTMLfile)
	HTML.title(bounds$degrate,HR=4,file=HTMLfile)
	
	#transformations used for the fitting
	if(length(lHVDM$transforms)==0) HTML.title("Transformations applied for the fitting (screening step): None",HR=4,file=HTMLfile)
	else{
		HTML.title("Transformations applied for the fitting (screening step):",HR=4,file=HTMLfile)
		HTML(lHVDM$transforms,file=HTMLfile)
		HTML(names(lHVDM$transforms),file=HTMLfile)
	}
	
	HTML("<a href=#Top>Top</a>",file=HTMLfile)
	HTML("<hr  width=100% size=20 noshade>",file=HTMLfile)

	#results section
	HTML(x="<a name=results></a>",file=HTMLfile)
	HTML.title(x="3) Results of the screening",HR=2,file=HTMLfile)
	
	HTML("<a href=#Top>Top</a> <a href=#Input>1) Inputs-training</a>  <a href=#Inputsc>2) Inputs-screening</a>  <a href=#class1>4) Putative targets</a>",file=HTMLfile)
	
	class1<-lHVDM$results
	class1<-class1[class1$class1,c("model_score","sens_z_score","basal","sensitivity","degradation")]
	n<-length(class1$model_score)
    bounds<-lHVDM$class1bounds	
	nzscore<-sum(lHVDM$results$sens_z_score>=bounds$zscore)
	nmodsco<-sum(lHVDM$results$model_score<=bounds$modelscore)
	ndegrat<-sum((lHVDM$results$degradation>=bounds$degrate[1]) & (lHVDM$results$degradation<=bounds$degrate[2]))
	nhessia<-sum(lHVDM$results$hess_rank>=3)
	HTML.title(paste(nzscore,"genes pass the sensitivity Z-score criteria"),HR=4,file=HTMLfile)
	HTML.title(paste(nmodsco,"genes pass model score criteria"),HR=4,file=HTMLfile)
	HTML.title(paste(ndegrat,"genes pass the degradation rate criteria"),HR=4,file=HTMLfile)
	HTML.title(paste("For",nhessia,"genes, the Hessian has full rank"),HR=4,file=HTMLfile)
	HTML.title(paste(n,"genes pass all the criteria above / are putative",lHVDM$tHVDM$par$activators,"targets (see next section)"),HR=4,file=HTMLfile)
	
	#Include the Figure here
	.screenoutput(HVDM=lHVDM,file=HTMLfile,rname=prefix)
	
	HTML("<a href=#Top>Top</a>",file=HTMLfile)
	HTML("<hr  width=100% size=20 noshade>",file=HTMLfile)

	#list of putative genes section
	HTML(x="<a name=class1></a>",file=HTMLfile)
	HTML.title(x=paste("4) List of putative",lHVDM$tHVDM$par$activators,"targets"),HR=2,file=HTMLfile)

	HTML("<a href=#Top>Top</a> <a href=#Input>1) Inputs-training</a>  <a href=#Inputsc>2) Inputs-screening</a>   <a href=#results>3) Results</a>",file=HTMLfile)
	HTML.title("the genes are ranked by decreasing sensitivity Z-score",HR=4,file=HTMLfile)
	
	class1<-lHVDM$results
	class1<-class1[class1$class1,c("model_score","sens_z_score","basal","sensitivity","degradation")]
	n<-length(class1$model_score)
	class1$rank<-c(1:n)
	HTML(class1[,c("rank","model_score","sens_z_score","basal","sensitivity","degradation")],file=HTMLfile,align="left")
	
	HTML("<a href=#Top>Top</a>",file=HTMLfile)
	
	#after creating the report, this function returns the name of the HTML file and its directory location
	res<-vector("list")
	res$file<-HTMLfile
	res$directory<-getwd()
	res
}

.screenoutput<-function(HVDM,file,rname){
		
	dirname<-paste(rname,.SUFFX.rHVDM,sep="")

	scrlist<-HVDM$results
	bounds<-HVDM$class1bounds
	degbnds<-bounds$degrate
	scrlist<-scrlist[(scrlist$degradation>=degbnds[1]) & (scrlist$degradation<=degbnds[2]) & !(scrlist$sens_z_score==-1),c("model_score","sens_z_score")]
	trlist<-HVDM$tHVDM$itgenes[,c("model_score","sens_z_score")]
	
	pngname<-"Z_vs_model_score.png"
	png(pngname,width=700,height=450)
	plot(scrlist$model_score,scrlist$sens_z_score,log="x",xlab="Model score",ylab="Sensitivity Z-score")
	lines(trlist$model_score,trlist$sens_z_score,type="p",col="green")
	xrng<-range(trlist$model_score,scrlist$model_score)
	yrng<-range(scrlist$model_score,scrlist$sens_z_score)
	lines(rep(bounds$modelscore,2),yrng,col="red")
	lines(xrng,rep(bounds$zscore,2),col="red")
	dev.off()
	pngname<-.moveindirectory(filename=pngname,dirname=dirname)
	HTMLInsertGraph(GraphFileName=pngname,Align="left",file=file,WidthHTML=700,
					Caption="Training genes are in green. Only those genes that pass both the degradation rate criteria and whose Hessian has full rank are plotted above.")
}

.greport<-function(sHVDM,name){
	#filenames and directory management
	genename<-rownames(sHVDM$par$genemap)
	if (missing(name)){
		prefix<-paste("HVDM_",as.character(floor(runif(1,100000,999999))),"__",genename,sep="")	}
	else{
		prefix<-name
	}
	HTMLfile<-paste(prefix,".HTML",sep="")
	dirname<-paste(prefix,.SUFFX.rHVDM,sep="")
	dir.create(dirname)

	#header section
	HTML(x="<a name=Top></a>",file=HTMLfile,append=FALSE)
	HTML.title(paste("Report for HVDM fitgene() step (",genename,")",sep=""),HR=1,file=HTMLfile)
	
	#links
	HTML("<a href=#Input>1) Inputs</a>  <a href=#Modelsco>2) Model score</a>   <a href=#Parfit>3) Parameter fit</a>  <a href=#Genefit>4) Comparison of model and data, gene by gene</a>",file=HTMLfile)
	
	HTML("<hr  width=100% size=20 noshade>",file=HTMLfile)
	
	HTML(x="<a name=Input></a>",file=HTMLfile)
	HTML.title(x="1) Parameters input to the fitgene() command",HR=2,file=HTMLfile)
	HTML("<a href=#Top>Top</a>   <a href=#Modelsco>2) Model score</a>   <a href=#Parfit>3) Parameter fit</a>  <a href=#Genefit>4) Comparison of model and data, gene by gene</a>",file=HTMLfile)
	#name of expression set used
	HTML.title(x=paste("Expression set:",sHVDM$eset),HR=4,file=HTMLfile)
	
	#name of training object
	HTML.title(paste("training set object used:",sHVDM$tHVDMname),HR=4,file=HTMLfile)
	
	#name of transcription factor
	def<-""
	if(sHVDM$par$activators=="trfact1") def<-"(default name)"
	HTML.title(paste("Transcription factor name:",sHVDM$par$activators,def),HR=4,file=HTMLfile)
	
	#set of training genes
	HTML.title(x="Training genes:",HR=4,file=HTMLfile)
	HTML(x=rownames(sHVDM$tset),file=HTMLfile,align="left")
	
	#transformations used for the fitting
	if(length(sHVDM$distribute$transf)==0) HTML.title("Transformations applied for the fitting: None",HR=4,file=HTMLfile)
	else{
		HTML.title("Transformations applied for the fitting:",HR=4,file=HTMLfile)
		HTML(sHVDM$distribute$transforms,file=HTMLfile)
		HTML(names(sHVDM$distribute$transforms),file=HTMLfile)
	}
	
	#first guess:
	HTML.title(x="First guess:",HR=4,file=HTMLfile)
	HTML(x=sHVDM$fguess,file=HTMLfile,align="left")
	
	HTML("<hr  width=100% size=20 noshade>",file=HTMLfile)
	
	#chisquared, qqplot and score distribution
	HTML(x="<a name=Modelsco></a>",file=HTMLfile)
	HTML.title(x="2) Model score",HR=2,file=HTMLfile)
	HTML("<a href=#Top>Top</a> <a href=#Input>1) Inputs</a>    <a href=#Parfit>3) Parameter fit</a>  <a href=#Genefit>4) Comparison of model and data, gene by gene</a>",file=HTMLfile)
	.scoreoutput(HVDM=sHVDM,file=HTMLfile,rname=prefix)
	HTML("<a href=#Top>Top</a>",file=HTMLfile)
	HTML("<hr  width=100% size=20 noshade>",file=HTMLfile)
	
	#parameter fit and hessian stuff
	HTML(x="<a name=Parfit></a>",file=HTMLfile)
	HTML.title(x="3) Parameter fit",HR=2,file=HTMLfile)
	HTML("<a href=#Top>Top</a> <a href=#Input>1) Inputs</a>  <a href=#Modelsco>2) Model score</a>     <a href=#Genefit>4) Comparison of model and data, gene by gene</a>",file=HTMLfile)
	if(length(sHVDM$distribute$free)==3){
		.parameteroutputg(tHVDM=sHVDM,file=HTMLfile,prefix=prefix)
	}
	else{
		.parameteroutputg.nl(tHVDM=sHVDM,file=HTMLfile,prefix=prefix)
	}
	HTML("<a href=#Top>Top</a>",file=HTMLfile)
	HTML("<hr  width=100% size=20 noshade>",file=HTMLfile)
	
	#comparison of model prediction and data, gene by gene
	HTML(x="<a name=Genefit></a>",file=HTMLfile)
	HTML.title(x=paste("4) Comparison of model fit and data for ",genename,sep=""),HR=2,file=HTMLfile)
	HTML("<a href=#Top>Top</a> <a href=#Input>1) Inputs</a>  <a href=#Modelsco>2) Model score</a>   <a href=#Parfit>3) Parameter fit</a>  ",file=HTMLfile)
	HTML.title(x="black: data + 95% confidence intervals / red: model fit",HR=4,file=HTMLfile)
	.allgenesoutput(tHVDM=sHVDM,file=HTMLfile,rname=prefix)
	HTML("<a href=#Top>Top</a>",file=HTMLfile)
	
	#after creating the report, this function returns the name of the HTML file and its directory location
	res<-vector("list")
	res$file<-HTMLfile
	res$directory<-getwd()
	res
}

.treport<-function(tHVDM,name){
	
	#filenames and directory management
	if (missing(name)){
		prefix<-paste("HVDM_",as.character(floor(runif(1,100000,999999))),"_training",sep="")	}
	else{
		prefix<-name
	}
	HTMLfile<-paste(prefix,".HTML",sep="")
	dirname<-paste(prefix,.SUFFX.rHVDM,sep="")
	dir.create(dirname)
	
	#header section
	HTML(x="<a name=Top></a>",file=HTMLfile,append=FALSE)
	HTML.title("Report for HVDM training() step",HR=1,file=HTMLfile)
	
	#links
	HTML("<a href=#Input>1) Inputs</a>  <a href=#Modelsco>2) Model score</a>   <a href=#Parfit>3) Parameter fit</a>  <a href=#Genefit>4) Comparison of model and data, gene by gene</a>",file=HTMLfile)
	
	HTML("<hr  width=100% size=20 noshade>",file=HTMLfile)
	
	HTML(x="<a name=Input></a>",file=HTMLfile)
	HTML.title(x="1) Parameters input to the training() command",HR=2,file=HTMLfile)
	HTML("<a href=#Top>Top</a>   <a href=#Modelsco>2) Model score</a>   <a href=#Parfit>3) Parameter fit</a>  <a href=#Genefit>4) Comparison of model and data, gene by gene</a>",file=HTMLfile)
	#name of expression set used
	HTML.title(x=paste("Expression set:",tHVDM$eset),HR=4,file=HTMLfile)
	
	#name of transcription factor
	def<-""
	if(tHVDM$par$activators=="trfact1") def<-"(default name)"
	HTML.title(paste("Transcription factor name:",tHVDM$par$activators,def),HR=4,file=HTMLfile)
	
	#transformations used for the fitting
	if(length(tHVDM$distribute$transf)==0) HTML.title("Transformations applied for the fitting: None",HR=4,file=HTMLfile)
	else{
		HTML.title("Transformations applied for the fitting:",HR=4,file=HTMLfile)
		HTML(tHVDM$distribute$transforms,file=HTMLfile)
		HTML(names(tHVDM$distribute$transforms),file=HTMLfile)
	}
	
	#set of training genes
	HTML.title(x="Training genes:",HR=4,file=HTMLfile)
	HTML(x=rownames(tHVDM$par$genemap),file=HTMLfile,align="left")
	
	#name of anchoring gene
	HTML.title(x="Anchoring gene:",HR=4,file=HTMLfile)
	anchgene<-rownames(tHVDM$par$genemap)[1]
	testname<-paste(anchgene,"Dj",sep=".")
	if (testname %in% tHVDM$distribute$free){
		HTML(x=paste(anchgene,": not anchored"),file=HTMLfile)
		HTML(x="<i> degradation rate has not been given, try parameter degrate in function training() </i>",file=HTMLfile)
	}
	else HTML(x=paste(anchgene,", degradation rate: ",tHVDM$par$parameters[testname]," (unit time)^(-1)",sep=""),file=HTMLfile)
	
	#Pheno data
	HTML.title(x="Pheno data:",HR=4,file=HTMLfile)
	HTML(x=tHVDM$pdata,file=HTMLfile,align="left")
	HTML("<a href=#Top>Top</a>",file=HTMLfile)
	HTML("<hr  width=100% size=20 noshade>",file=HTMLfile)
	
	
	#chisquared, qqplot and score distribution
	HTML(x="<a name=Modelsco></a>",file=HTMLfile)
	HTML.title(x="2) Model score",HR=2,file=HTMLfile)
	HTML("<a href=#Top>Top</a> <a href=#Input>1) Inputs</a>    <a href=#Parfit>3) Parameter fit</a>  <a href=#Genefit>4) Comparison of model and data, gene by gene</a>",file=HTMLfile)
	.scoreoutput(HVDM=tHVDM,file=HTMLfile,rname=prefix)
	HTML("<a href=#Top>Top</a>",file=HTMLfile)
	HTML("<hr  width=100% size=20 noshade>",file=HTMLfile)
	
	#parameter fit and hessian stuff
	HTML(x="<a name=Parfit></a>",file=HTMLfile)
	HTML.title(x="3) Parameter fit",HR=2,file=HTMLfile)
	HTML("<a href=#Top>Top</a> <a href=#Input>1) Inputs</a>  <a href=#Modelsco>2) Model score</a>     <a href=#Genefit>4) Comparison of model and data, gene by gene</a>",file=HTMLfile)
	.parameteroutput(tHVDM=tHVDM,file=HTMLfile,rname=prefix)
	HTML("<a href=#Top>Top</a>",file=HTMLfile)
	HTML("<hr  width=100% size=20 noshade>",file=HTMLfile)
	
	#comparison of model prediction and data, gene by gene
	HTML(x="<a name=Genefit></a>",file=HTMLfile)
	HTML.title(x="4) Comparison of model fit and data, gene by gene",HR=2,file=HTMLfile)
	HTML("<a href=#Top>Top</a> <a href=#Input>1) Inputs</a>  <a href=#Modelsco>2) Model score</a>   <a href=#Parfit>3) Parameter fit</a>  ",file=HTMLfile)
	HTML.title(x="black: data + 95% confidence intervals / red: model fit",HR=4,file=HTMLfile)
	.allgenesoutput(tHVDM=tHVDM,file=HTMLfile,rname=prefix)
	HTML("<a href=#Top>Top</a>",file=HTMLfile)
	
	#after creating the report, this function returns the name of the HTML file and its directory location
	res<-vector("list")
	res$file<-HTMLfile
	res$directory<-getwd()
	res
}

.parameteroutputg<-function(tHVDM,file,prefix){
	#name of directory
	dirname<-paste(prefix,.SUFFX.rHVDM,sep="")

	#Hessian stuff
	H<-tHVDM$results$hessian
	nfp<-length(tHVDM$distribute$free)
	rnk<-.hessrank(H)
	if (rnk==nfp){
		HTML.title(x="NB: The hessian has full rank",HR=5,file=file)
		ciS=TRUE
	}
	else{
		HTML.title(x="NB: The hessian is singular",HR=5,file=file)
		HTML.title(x="Confidence intervals will not be plotted",HR=5,file=file)
		HTML.title(x="The covariance matrix can not be computed",HR=5,file=file)
		ciS=FALSE
	}
	
	#extract transformed bounds (if applicable)
	if (ciS){
		vcov<-solve(H)
		sdev<-diag(vcov)^0.5
		work<-tHVDM
		central<-.exportfree(work)
		nams<-names(central)
		allupbnd<-(.importfree(HVDM=tHVDM,x=central[nams]+1.96*sdev[nams]))$par$parameters
		alllobnd<-(.importfree(HVDM=tHVDM,x=central[nams]-1.96*sdev[nams]))$par$parameters
		#computation of zscore sensitivity
		nameex<-paste(names(tHVDM$scores$bygene),"Sj",sep=".")#replace Sj with Vj in case the nuber of parameters exceeds 3
		upbnd<-(.importfree(HVDM=tHVDM,x=central[nams]+sdev[nams]))$par$parameters
		lobnd<-(.importfree(HVDM=tHVDM,x=central[nams]-sdev[nams]))$par$parameters
		zscore<-2*central[nameex]/(upbnd[nameex]-lobnd[nameex])
		rm(work)
	}
	else zscore<--1
	
	
	#kinetic parameters
	HTML.title(x="Kinetic parameters",HR=3,file=file)
	genes<-rownames(tHVDM$par$genemap)
	kpar<-tHVDM$par$linear  
	pngname<-"kinetic_parms.png"
	png(pngname,width=800,height=300)
	old=par(mfrow=c(1,3),las=3,mar=c(7,4,4,2)+0.1)
	for(k in kpar){  #[tbc} this is for a normal linear model only
		if (k=="Bj"){
			cen<-tHVDM$tset$basal
			upb<-tHVDM$tset$basal_u
			lob<-tHVDM$tset$basal_d
		}
		else if (k=="Sj"){
			cen<-tHVDM$tset$sensitivity
			upb<-tHVDM$tset$sensitivity_u
			lob<-tHVDM$tset$sensitivity_d
		}
		else if (k=="Dj"){
			cen<-tHVDM$tset$degradation
			upb<-tHVDM$tset$degradation_u
			lob<-tHVDM$tset$degradation_d
		}
		tgenes<-rownames(tHVDM$tset)
		parnames<-paste(genes,k,sep=".")
		col<-c(rep("green",length(parnames)),rep("grey",length(tgenes)))
		vals<-c(tHVDM$par$parameters[parnames],cen)
		#determine bar colours
		#prepare plot
		names(vals)<-c(genes,tgenes)
		if (ciS){
			upbnd<-c(allupbnd[parnames],upb)
			lobnd<-c(alllobnd[parnames],lob)
			rng<-range(0,vals,upbnd,lobnd)
			tbp<-(lobnd!=upbnd)
			x<-barplot(vals,col=col,main=paste("parameter",k),ylim=rng,ylab="parameter value")
			arrows(x[tbp],upbnd[tbp],x[tbp],lobnd[tbp],angle=90,code=3,length=0.05,col="red")
		}
		else barplot(vals,col=col,main=paste("parameter",k),ylab="parameter value")
	}
	par(old)
	dev.off()
	pngname<-.moveindirectory(filename=pngname,dirname=dirname)
	HTMLInsertGraph(GraphFileName=pngname,Align="left",file=file,WidthHTML=800,
					Caption="Green:gene under review, Grey:training genes, fitted individually. If the hessian is singular, corresponding confidence intervals cannot be determined.")
	HTMLbr(x=1,file=file)
	
	#zscore plot
	zscores<-c(zscore,tHVDM$tset$sens_z_score)
	names(zscores)<-c(names(tHVDM$scores$bygene),rownames(tHVDM$tset))
	colours<-rep("grey",length(zscores))
	colours[1]<-"green"
	HTML.title(x="Sensitivity Z-score",HR=3,file=file)
	pngname<-"sens_zscores.png"
	png(pngname,width=400,height=350)
	old<-par(las=3,mar=c(7,4,4,2)+0.1)
	barplot(zscores,col=colours,ylab="sensitivity Z-score")
	par(old)
	dev.off()
	pngname<-.moveindirectory(filename=pngname,dirname=dirname)
	HTMLInsertGraph(GraphFileName=pngname,Align="left",file=file,WidthHTML=400,
					Caption="Z-scores of the training set are in grey (the gene under review is in green). A Z-score equal to -1 indicates that it could not be computed because of singularity problems.")
	HTMLbr(x=1,file=file)
		
	#eigenvalues of variance-covariance matrix (if applicable)
	if (ciS){
		HTML.title(x="Eigenvalues of the covariance matrix",HR=3,file=file)
		ev<-eigen(vcov)$values
		pngname<-"hessian_evalues.png"
		png(pngname,width=300,height=200)
		old=par(mar=c(2,2,2,1)+0.1,omi=c(0,0,0,0))
		barplot(ev,xlab="eigenvalue rank",ylab="eigenvalue")
		par(old)
		dev.off()
		pngname<-.moveindirectory(filename=pngname,dirname=dirname)
		HTMLInsertGraph(GraphFileName=pngname,Align="left",file=file,WidthHTML=300,
					Caption="The covariance matrix is estimated by inverting the Hessian.")
		HTMLbr(x=1,file=file)
	}
	
	#output correlation matrix
	if (ciS){ #print out correlation matrix (if applicable)
		corr<-round(cov2cor(vcov),3)
		HTML.title(x="Correlation matrix of the free parameters",HR=3,file=file)
		dim<-length(tHVDM$distribute$free)
		#subdimension 6 or 7, make sure remainder is not 1
		subd<-6
		if(dim%%subd==1) subd<-7
		nstep<-ceiling(dim/subd)
		for(i in 1:nstep){
			colrng<-c(1,subd)+(i-1)*subd
			colrng[2]<-min(dim,colrng[2])
			for (j in 1:nstep){
				linrng<-c(1,subd)+(j-1)*subd
				linrng[2]<-min(dim,linrng[2])
				HTML.title(x=paste("lines",linrng[1],"-",linrng[2]," / ","columns",colrng[1],"-",colrng[2]),HR=5,file=file)
				HTML(corr[c(linrng[1]:linrng[2]),c(colrng[1]:colrng[2])],file=file,align="left")
			}
		}
		HTMLbr(x=1,file=file)
	}
}

.parameteroutputg.nl<-function(tHVDM,file,prefix){
	#name of directory
	dirname<-paste(prefix,.SUFFX.rHVDM,sep="")

	#Hessian stuff
	H<-tHVDM$results$hessian
	nfp<-length(tHVDM$distribute$free)
	rnk<-.hessrank(H)
	if (rnk==nfp){
		HTML.title(x="NB: The hessian has full rank",HR=5,file=file)
		ciS=TRUE
	}
	else{
		HTML.title(x="NB: The hessian is singular",HR=5,file=file)
		HTML.title(x="Confidence intervals will not be plotted",HR=5,file=file)
		HTML.title(x="The covariance matrix can not be computed",HR=5,file=file)
		ciS=FALSE
	}
	
	#extract transformed bounds (if applicable)
	if (ciS){
		vcov<-solve(H)
		sdev<-diag(vcov)^0.5
		work<-tHVDM
		central<-.exportfree(work)
		nams<-names(central)
		allupbnd<-(.importfree(HVDM=tHVDM,x=central[nams]+1.96*sdev[nams]))$par$parameters
		alllobnd<-(.importfree(HVDM=tHVDM,x=central[nams]-1.96*sdev[nams]))$par$parameters
		#computation of zscore of Vj (no more sensitivity with nonlinear model)
		nameex<-paste(names(tHVDM$scores$bygene),"Vj",sep=".")
		upbnd<-(.importfree(HVDM=tHVDM,x=central[nams]+sdev[nams]))$par$parameters
		lobnd<-(.importfree(HVDM=tHVDM,x=central[nams]-sdev[nams]))$par$parameters
		zscore<-2*central[nameex]/(upbnd[nameex]-lobnd[nameex])
		rm(work)
	}
	else zscore<--1
	
	
	#kinetic parameters
	HTML.title(x="Kinetic parameters",HR=3,file=file)
	genes<-rownames(tHVDM$par$genemap)
	kpar<-tHVDM$par$hill
	if (length(tHVDM$distribute$free)==4){
		kpar<-tHVDM$par$MM
	}
	pngname<-"kinetic_parms.png"
	png(pngname,width=800,height=600)
	old=par(mfrow=c(2,3),las=3,mar=c(7,4,4,2)+0.1)
	for(k in kpar){
		cen<-tHVDM$tset[,k]
		upb<-tHVDM$tset[,paste(k,"u",sep="_")]
		lob<-tHVDM$tset[,paste(k,"d",sep="_")]
		tgenes<-rownames(tHVDM$tset)
		parnames<-paste(genes,k,sep=".")
		col<-c(rep("green",length(parnames)),rep("grey",length(tgenes)))
		vals<-c(tHVDM$par$parameters[parnames],cen)
		#determine bar colours
		#prepare plot
		names(vals)<-c(genes,tgenes)
		if (ciS){
			upbnd<-c(allupbnd[parnames],upb)
			lobnd<-c(alllobnd[parnames],lob)
			rng<-range(0,vals,upbnd,lobnd)
			tbp<-(lobnd!=upbnd)
			x<-barplot(vals,col=col,main=paste("parameter",k),ylim=rng,ylab="parameter value")
			arrows(x[tbp],upbnd[tbp],x[tbp],lobnd[tbp],angle=90,code=3,length=0.05,col="red")
		}
		else barplot(vals,col=col,main=paste("parameter",k),ylab="parameter value")
	}
	par(old)
	dev.off()
	pngname<-.moveindirectory(filename=pngname,dirname=dirname)
	HTMLInsertGraph(GraphFileName=pngname,Align="left",file=file,WidthHTML=800,
					Caption="Green:gene under review, Grey:training genes, fitted individually. If the hessian is singular, corresponding confidence intervals cannot be determined.")
	HTMLbr(x=1,file=file)
	
	#zscore plot
	zscores<-c(zscore,tHVDM$tset$Vj_z_score)
	names(zscores)<-c(names(tHVDM$scores$bygene),rownames(tHVDM$tset))
	colours<-rep("grey",length(zscores))
	colours[1]<-"green"
	HTML.title(x="Vj Z-score",HR=3,file=file)
	pngname<-"sens_zscores.png"
	png(pngname,width=400,height=350)
	old<-par(las=3,mar=c(7,4,4,2)+0.1)
	barplot(zscores,col=colours,ylab="sensitivity Z-score")
	par(old)
	dev.off()
	pngname<-.moveindirectory(filename=pngname,dirname=dirname)
	HTMLInsertGraph(GraphFileName=pngname,Align="left",file=file,WidthHTML=400,
					Caption="Z-scores of the training set are in grey (the gene under review is in green). A Z-score equal to -1 indicates that it could not be computed because of singularity problems.")
	HTMLbr(x=1,file=file)
		
	#eigenvalues of variance-covariance matrix (if applicable)
	if (ciS){
		HTML.title(x="Eigenvalues of the covariance matrix",HR=3,file=file)
		ev<-eigen(vcov)$values
		pngname<-"hessian_evalues.png"
		png(pngname,width=300,height=200)
		old=par(mar=c(2,2,2,1)+0.1,omi=c(0,0,0,0))
		barplot(ev,xlab="eigenvalue rank",ylab="eigenvalue")
		par(old)
		dev.off()
		pngname<-.moveindirectory(filename=pngname,dirname=dirname)
		HTMLInsertGraph(GraphFileName=pngname,Align="left",file=file,WidthHTML=300,
					Caption="The covariance matrix is estimated by inverting the Hessian.")
		HTMLbr(x=1,file=file)
	}
	
	#output correlation matrix
	if (ciS){ #print out correlation matrix (if applicable)
		corr<-round(cov2cor(vcov),3)
		HTML.title(x="Correlation matrix of the free parameters",HR=3,file=file)
		dim<-length(tHVDM$distribute$free)
		#subdimension 6 or 7, make sure remainder is not 1
		subd<-6
		if(dim%%subd==1) subd<-7
		nstep<-ceiling(dim/subd)
		for(i in 1:nstep){
			colrng<-c(1,subd)+(i-1)*subd
			colrng[2]<-min(dim,colrng[2])
			for (j in 1:nstep){
				linrng<-c(1,subd)+(j-1)*subd
				linrng[2]<-min(dim,linrng[2])
				HTML.title(x=paste("lines",linrng[1],"-",linrng[2]," / ","columns",colrng[1],"-",colrng[2]),HR=5,file=file)
				HTML(corr[c(linrng[1]:linrng[2]),c(colrng[1]:colrng[2])],file=file,align="left")
			}
		}
		HTMLbr(x=1,file=file)
	}
}

.parameteroutput<-function(tHVDM,file,rname){
	#directory name
	dirname<-paste(rname,.SUFFX.rHVDM,sep="")
	
	#Hessian stuff
	H<-tHVDM$results$hessian
	nfp<-length(tHVDM$distribute$free)
	rnk<-.hessrank(H)
	if (rnk==nfp){
		HTML.title(x="NB: The hessian has full rank",HR=5,file=file)
		ciS=TRUE
	}
	else{
		HTML.title(x="NB: The hessian is singular",HR=5,file=file)
		HTML.title(x="Confidence intervals will not be plotted",HR=5,file=file)
		HTML.title(x="The covariance matrix can not be computed",HR=5,file=file)
		ciS=FALSE
	}
	
	#extract transformed bounds (if applicable)
	if (ciS){
		vcov<-solve(H)
		sdev<-1.96*diag(vcov)^0.5
		work<-tHVDM
		central<-.exportfree(work)
		nams<-names(central)
		allupbnd<-(.importfree(HVDM=tHVDM,x=central[nams]+sdev[nams]))$par$parameters
		alllobnd<-(.importfree(HVDM=tHVDM,x=central[nams]-sdev[nams]))$par$parameters
		rm(work)
	}
	
	#kinetic parameters
	HTML.title(x="Kinetic parameters",HR=3,file=file)
	genes<-rownames(tHVDM$par$genemap)
	kpar<-c("bidon")
	for(gene in genes){
		modell<-as.character(tHVDM$par$genemap[gene,"model"])
		params<-tHVDM$par[modell][[1]]
		kpar<-union(kpar,params)
	}
	kpar<-setdiff(kpar,c("bidon"))
	nkpar<-length(kpar)
	grid<-c(ceiling(nkpar/3),min(nkpar,3))  
	pngname<-"kinetic_pars.png"
	png(pngname,width=800/3*grid[2],height=300*grid[1])
	old=par(mfrow=grid,las=3,mar=c(7,4,4,2)+0.1)
	for(k in kpar){  #[tbc] this should work for any model
		parnames<-paste(genes,k,sep=".")
		genez<-genes
		names(genez)<-parnames
		parnames<-intersect(parnames,names(tHVDM$par$parameters))
		genez<-genez[parnames]
		col<-rep("green",length(parnames))
		names(col)<-parnames
		vals<-tHVDM$par$parameters[parnames]
		#determine bar colours
		free<-vals[tHVDM$distribute$free]
		free<-names(free[!is.na(free)])
		col[free]<-"grey"
		if (ciS){ #plot with error bars
			#prepare plot
			names(vals)<-genez[parnames]
			upbnd<-allupbnd[parnames]
			lobnd<-alllobnd[parnames]
			rng<-range(0,vals,upbnd,lobnd)
			tbp<-(lobnd!=upbnd)
			x<-barplot(vals,col=col,main=paste("parameter",k),ylim=rng,ylab="parameter value")
			arrows(x[tbp],upbnd[tbp],x[tbp],lobnd[tbp],angle=90,code=3,length=0.05,col="red")
		}
		else{ #plot without error bars
			names(vals)<-genez[parnames]
			barplot(vals,col=col,main=paste("parameter",k),ylab="parameter value")
		}
	}
	par(old)
	dev.off()
	pngname<-.moveindirectory(filename=pngname,dirname=dirname)
	if (ciS) HTMLInsertGraph(GraphFileName=pngname,Align="left",file=file,WidthHTML=800,
					Caption="Fixed parameters are indicated in green and have no error bars associated.\n
							Error bars indicate best fit +/- 1.96*(standard deviation), as extracted from the Hessian.\n
							If a parameter transform has been applied for the fit, the corresponding error bar is transformed accordingly.")
	else HTMLInsertGraph(GraphFileName=pngname,Align="left",file=file,WidthHTML=800,
					Caption="No error bars are given, as the Hessian is Singular.\n
							Fixed parameters are indicated in green.")
	HTMLbr(x=1,file=file)
	
	#transcription factor activity
	HTML.title(x=paste("Transcription factor activity (or concentration if non-linear model) of",tHVDM$par$activators),HR=3,file=file)
	pngname<-"trfact_activity.png"
	tc<-tHVDM$tc
	ntc<-length(tc)
	grid<-c(ceiling(ntc/3),min(ntc,3))
	png(pngname,width=800/3*grid[2],height=300*grid[1])
	old<-par(mfrow=grid)
	trfac<-tHVDM$par$activators
	for(i in 1:ntc){
		exprep<-paste(tc[[i]]$experiment,tc[[i]]$replicate,sep=".")
		time<-tc[[i]]$time
		pnames<-paste(trfac,exprep,time,sep=".")
		vals<-tHVDM$par$parameters[pnames]
		col<-rep("green",length(time))
		names(col)<-pnames
		free<-vals[tHVDM$distribute$free]
		free<-names(free[!is.na(free)])
		col[free]<-"grey"
		if (ciS){#plot with error bars
			upbnd<-allupbnd[pnames]
			lobnd<-alllobnd[pnames]
			tbp<-(upbnd!=lobnd)
			rng<-range(0,vals,upbnd,lobnd)
			plot(x=time,y=vals,col=col, main=exprep,ylim=rng,pch="x",ylab="activity",xlab="time")
			lines(time,vals,col="red")
			if(sum(tbp)>0) arrows(time[tbp],upbnd[tbp],time[tbp],lobnd[tbp],angle=90,code=3,length=0.05)
		}
		else{ #plot without error bars
			plot(x=time,y=vals,col=col, main=exprep,pch="x",ylab="activity",xlab="time")
			lines(time,vals,col="red")
		}
	}
	rm(tc)
	par(old)
	dev.off()
	pngname<-.moveindirectory(filename=pngname,dirname=dirname)
	if (ciS) HTMLInsertGraph(GraphFileName=pngname,Align="left",file=file,WidthHTML=800/3*grid[2],
					Caption="Fixed time points are indicated with a green cross and have no error bars associated.\n
							Error bars indicate best fit +/- 1.96*(standard deviation), as extracted from the Hessian.\n
							If a transform has been applied for the fit, the corresponding error bar is transformed accordingly.")
	else HTMLInsertGraph(GraphFileName=pngname,Align="left",file=file,WidthHTML=800/3*grid[2],
					Caption="No error bars are given, as the Hessian is Singular.\n
							Fixed time points are indicated with a green cross.")
	HTMLbr(x=1,file=file)
	
	#zscores
	
	#eigenvalues of variance-covariance matrix (if applicable)
	if (ciS){
		HTML.title(x="Eigenvalues of the covariance matrix",HR=3,file=file)
		ev<-eigen(vcov)$values
		pngname<-"hessian_evalues.png"
		png(pngname,width=500,height=250)
		old=par(mar=c(2,2,2,1)+0.1,omi=c(0,0,0,0))
		barplot(ev,xlab="eigenvalue rank",ylab="eigenvalue")
		par(old)
		dev.off()
		pngname<-.moveindirectory(filename=pngname,dirname=dirname)
		HTMLInsertGraph(GraphFileName=pngname,Align="left",file=file,WidthHTML=500,
					Caption="The covariance matrix is estimated by inverting the Hessian. A single dominant eigenvalue indicates that the system is close to singularity: if this has not been done yet, try anchoring")
		HTMLbr(x=1,file=file)
	}
	
	#output correlation matrix
	if (ciS){ #print out correlation matrix (if applicable)
		corr<-round(cov2cor(vcov),3)
		HTML.title(x="Correlation matrix of the free parameters",HR=3,file=file)
		dim<-length(tHVDM$distribute$free)
		#subdimension 6 or 7, make sure remainder is not 1
		subd<-6
		if(dim%%subd==1) subd<-7
		nstep<-ceiling(dim/subd)
		for(i in 1:nstep){
			colrng<-c(1,subd)+(i-1)*subd
			colrng[2]<-min(dim,colrng[2])
			for (j in 1:nstep){
				linrng<-c(1,subd)+(j-1)*subd
				linrng[2]<-min(dim,linrng[2])
				HTML.title(x=paste("lines",linrng[1],"-",linrng[2]," / ","columns",colrng[1],"-",colrng[2]),HR=5,file=file)
				HTML(corr[c(linrng[1]:linrng[2]),c(colrng[1]:colrng[2])],file=file,align="left")
			}
		}
		HTMLbr(x=1,file=file)
	}
}

.allgenesoutput<-function(tHVDM,file,rname){
	genes<-rownames(tHVDM$par$genemap)
	for(gene in genes) .geneoutput(HVDM=tHVDM,gene=gene,file=file,rname=rname)
}

.geneoutput<-function(HVDM,gene,file,rname){
	#directory name
	dirname<-paste(rname,.SUFFX.rHVDM,sep="")
	
	#the rest...
	HTML.title(x=gene,HR=3,file=file)
	ntc<-length(HVDM$tc)
	grid<-c(ceiling(ntc/3),min(3,ntc))
	if (HVDM$type=="indgene") pngname<-paste("data_and_model___",gene,".png",sep="")
	else if (HVDM$type=="training") pngname<-paste("training_data_and_model___",gene,".png",sep="")
	png(pngname,width=800/3*grid[2],height=300*grid[1])
	old<-par(mfrow=grid)
	for(i in 1:ntc){
		time<-HVDM$tc[[i]]$time
		exprep<-paste(HVDM$tc[[i]]$experiment,HVDM$tc[[i]]$replicate,sep=".")
		dnames<-paste(gene,exprep,time,sep=".")
		signal<-HVDM$dm$signal[dnames]
		ci<-1.96*HVDM$dm$sdev[dnames]
		model<-HVDM$dm$estimate[dnames]
		high<-max(model,signal+ci)
		plot(time,signal,main=exprep,ylab="concentration",ylim=c(0,high),pch="x",xlab="time")
		arrows(time,signal+ci,time,signal-ci,angle=90,code=3,length=0.05)
		lines(time,model,col="red")
	}
	par(old)
	dev.off()
	pngname<-.moveindirectory(filename=pngname,dirname=dirname)
	HTMLInsertGraph(GraphFileName=pngname,Align="left",file=file,WidthHTML=800/3*grid[2])
	HTMLbr(x=1,file=file)
}

.scoreoutput<-function(HVDM,file,rname){
	dirname<-paste(rname,.SUFFX.rHVDM,sep="")

	scores<-HVDM$scores
	if (HVDM$type=="training"){#if it is a training HVDM object, the output is richer
		
		#model score total, number of parameters/observations/degrees of freedom
		npars<-length(HVDM$distribute$free)
		nobs<-length(HVDM$dm$estimate)
		HTML.title(x=paste("number of observations:",nobs),HR=4,file=file)
		HTML.title(x=paste("number of parameters to fit:",npars),HR=4,file=file)
		HTML.title(x=paste("degrees of freedom:",nobs-npars),HR=4,file=file)
		HTML.title(x=paste("Model score:",round(scores$total,2)),HR=4,file=file)
		
		scores$total
		#model score: chi quare and qqplot
		df<-nobs-npars
		ps<-c(1:1000)/1001
		qs<-sort(c(qchisq(ps,df=df),scores$total*c(0.95,1.05)))
		pngname<-"training_s0.png"
		png(pngname,width=700,height=400)
		old<-par(mfrow=c(1,2),las=3,mar=c(7,4,4,2)+0.1)
		plot(qs,dchisq(qs,df=df),type="l",xlab="scores",ylab="frequency",main=paste("chi2, df =",df))
		points(rep(scores$total,2),c(0,dchisq(scores$total,df=df)),col="red",pch="X")
		devs<-scores$devs
		qnorm01<-qnorm(p=c(1:nobs)/(nobs+1))
		plot(qnorm01,sort(devs),pch="x",main="quantile-quantile plot",xlab="N(0,1)",ylab="observed deviations",col="red")
		lines(qnorm01,qnorm01)
		par(old)
		dev.off()
		pngname<-.moveindirectory(filename=pngname,dirname=dirname)
		HTMLInsertGraph(GraphFileName=pngname,Align="left",file=file,WidthHTML=700)
		HTMLbr(x=1,file=file)
	
		#Distribution of model score: general
		HTML.title(x="Distribution of model score",HR=4,file=file)
		pngname<-"training_s1.png"
		png(pngname,width=800,height=400)
		old<-par(mfrow=c(1,3),las=3,mar=c(7,4,4,2)+0.1)
		barplot(scores$bygene,main="by gene",crt=45,ylab="score")
		barplot(scores$bytc,main="by time course",ylab="score")
		barplot(scores$bytcpertp,main="by time course\n(average per time point)",ylab="score")
		par(old)
		dev.off()
		pngname<-.moveindirectory(filename=pngname,dirname=dirname)
		HTMLInsertGraph(GraphFileName=pngname,Align="left",file=file,WidthHTML=800)
		HTMLbr(x=1,file=file)
		
		#Distribution of model score per time course
		HTML.title(x="\nDistribution of model score within time course by time point",HR=4,file=file)
		pngname<-"training_s2.png"
		ntc<-length(scores$withintc)
		grid<-c(ceiling(ntc/3),min(3,ntc))
		png(pngname,width=800/3*grid[2],height=300*grid[1])
		old<-par(mfrow=grid)
		#determine plot range (same scale for each time course)
		rng<-c(0,0)
		for(tc in scores$withintc){
			rng<-range(rng,tc$score)
		}
		for(tc in scores$withintc){
			barplot(tc$score,main=tc$name,ylim=rng,ylab="score")
		}
		par(old)
		dev.off()
		pngname<-.moveindirectory(filename=pngname,dirname=dirname)
		HTMLInsertGraph(GraphFileName=pngname,Align="left",file=file,WidthHTML=800/3*grid[2])
	}
	else if (HVDM$type=="indgene"){#Note that here rname is no longer a number but the prefix for the file output
		
		genename<-names(scores$bygene)
		
		#model score total, number of parameters/observations/degrees of freedom
		npars<-length(HVDM$distribute$free)
		nobs<-length(HVDM$dm$estimate)
		HTML.title(x=paste("number of observations:",nobs),HR=4,file=file)
		HTML.title(x=paste("number of parameters to fit:",npars),HR=4,file=file)
		HTML.title(x=paste("degrees of freedom:",nobs-npars),HR=4,file=file)
		HTML.title(x=paste("Model score:",round(scores$total,2)),HR=4,file=file)
		
		scores$total
		#model score: chi quare and qqplot
		df<-nobs-npars
		ps<-c(1:1000)/1001
		qs<-sort(c(qchisq(ps,df=df),scores$total*c(0.95,1.05)))
		pngname<-"indivgene_s0.png"
		png(pngname,width=700,height=400)
		old<-par(mfrow=c(1,2),las=3,mar=c(7,4,4,2)+0.1)
		plot(qs,dchisq(qs,df=df),type="l",xlab="scores",ylab="frequency",main=paste("chi2, df =",df))
		points(rep(scores$total,2),c(0,dchisq(scores$total,df=df)),col="red",pch="X")
		devs<-scores$devs
		qnorm01<-qnorm(p=c(1:nobs)/(nobs+1))
		plot(qnorm01,sort(devs),pch="x",main="quantile-quantile plot",xlab="N(0,1)",ylab="observed deviations",col="red")
		lines(qnorm01,qnorm01)
		par(old)
		dev.off()
		pngname<-.moveindirectory(filename=pngname,dirname=dirname)
		HTMLInsertGraph(GraphFileName=pngname,Align="left",file=file,WidthHTML=700)
		HTMLbr(x=1,file=file)
		
		#Distribution of model score: general
		HTML.title(x=paste("Model score: comparison of ",genename," with training genes and time course distribution",sep=""),HR=4,file=file)
		pngname<-"indivgene_s1.png"
		png(pngname,width=800,height=400)
		old<-par(mfrow=c(1,3),las=3,mar=c(7,4,4,2)+0.1)
		tsetscore<-HVDM$tset$model_score
		names(tsetscore)<-rownames(HVDM$tset)
		genescore<-c(scores$bygene,tsetscore)
		colours<-rep("grey",length(genescore))
		colours[1]<-"green"
		barplot(genescore,main="by gene\n(grey: genes in the training set)",crt=45,col=colours,ylab="score")
		barplot(scores$bytc,main=paste("by time course for ",genename,sep=""),ylab="score")
		barplot(scores$bytcpertp,main=paste("by time course for ",genename,"\n(average per time point)",sep=""),ylab="score")
		par(old)
		dev.off()
		pngname<-.moveindirectory(filename=pngname,dirname=dirname)
		HTMLInsertGraph(GraphFileName=pngname,Align="left",file=file,WidthHTML=800)
		HTMLbr(x=1,file=file)
		
		#Distribution of model score per time course
		HTML.title(x=paste("\nDistribution of model score within time course by time point for ",genename,sep=""),HR=4,file=file)
		pngname<-"indivgene_s2.png"
		ntc<-length(scores$withintc)
		grid<-c(ceiling(ntc/3),min(ntc,3))
		png(pngname,width=800/3*grid[2],height=300*grid[1])
		old<-par(mfrow=grid)
		#determine plot range (same scale for each time course)
		rng<-c(0,0)
		for(tc in scores$withintc){
			rng<-range(rng,tc$score)
		}
		for(tc in scores$withintc){
			barplot(tc$score,main=tc$name,ylim=rng,ylab="score")
		}
		par(old)
		dev.off()
		pngname<-.moveindirectory(filename=pngname,dirname=dirname)
		HTMLInsertGraph(GraphFileName=pngname,Align="left",file=file,WidthHTML=800/3*grid[2])
		
	}
}

.moveindirectory<- function(filename,dirname){
	file.copy(from=filename,to=dirname,overwrite=TRUE)
	file.remove(filename)
	newfname<-paste(dirname,filename,sep=.Platform$file.sep)
	newfname
}