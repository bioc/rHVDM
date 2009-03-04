#user function(s) for rHVDM (should be OK for all types of models)

fitgene<-function(eset,gene,tHVDM,transforms,firstguess,criterion="BIC"){
	if (tHVDM$type=="training"){
		if(missing(transforms)) sHVDM<-fitgene.lin(eset=eset,gene=gene,tHVDM=tHVDM,firstguess=firstguess)
		else sHVDM<-fitgene.lin(eset=eset,gene=gene,tHVDM=tHVDM,transforms=transforms,firstguess=firstguess)
	}
	else{
		sHVDM<-.fitgene.best(eset=eset,gene=gene,tHVDM=tHVDM,criterion=criterion)
	}
	sHVDM$eset<-deparse(substitute(eset))
	sHVDM$tHVDMname<-deparse(substitute(tHVDM))
	sHVDM$type<-c("indgene")
	sHVDM
}

HVDMreport<-function(HVDM,name){
	if(HVDM$type=="training") res<-.treport(tHVDM=HVDM,name=name)
	else if(HVDM$type=="indgene") res<-.greport(sHVDM=HVDM,name=name)
	else if(HVDM$type=="screening") res<-.sreport(lHVDM=HVDM,name=name)
	else if(HVDM$type=="training.nl"){
		HVDM$type<-"training"
		res<-.treport(tHVDM=HVDM,name=name)
		HVDM$type<-"training.nl"
	}
	message(paste("the report",res$file,"was generated"))
	message(paste("in directory",res$directory))
}

HVDMcheck<-function(eset,pdata){
	if (missing(pdata)) pdata<-pData(eset)
	#check that the pheno data has all the required fields
	fieldnames<-names(pdata)
	nopdwarnings<-TRUE
	message("*** checking that pheno data has all the required fields:")
	if (!("experiment" %in% fieldnames)){
		message("   A) an experiment field is missing in the phenodata")
		nopdwarnings<-FALSE
	}
	if (!("replicate" %in% fieldnames)){
		message("   A) a replicate field is missing in the phenodata")
		nopdwarnings<-FALSE
	}
	if (!("time" %in% fieldnames)){
		message("   A) a time field is missing in the phenodata")
		nopdwarnings<-FALSE
	}
	if (nopdwarnings) message("					OK")
	else message("pdata NOT OK")
	timevecsOK<-TRUE
	if (nopdwarnings) timevecsOK<-.checktimevectors(pdata)
	else message("the time entries in pdata could not be checked")
	
	#check the names of pData can be found in the column names of expresssion
	expnamestest<-TRUE
	colnames<-dimnames(exprs(eset))[[2]]
	pdatanames<-rownames(pdata)
	message("*** checking that pdata and expression set coincide:")
	for (rn in pdatanames){
		if (!(rn %in% colnames)){
			message(paste("   C) could not find",rn,"data in the expression set"))
			expnamestest<-FALSE 
		}
	}
	if(expnamestest) message("					OK")
	else message("   one or more of the data points in the pdata could not be found in the expression set")
	#checking dynamic range
	message("*** checking dynamic range of expression values:")
	rng<-range(exprs(eset))
	rngOK<-TRUE
	if(rng[2]-rng[1]<20){
		message("   D) your data seems to be log-transformed, make sure that is it not")
		rngOK<-FALSE
	}
	else message("					OK")
	#checking the standard deviation
	sdOK<-TRUE
	sd<-assayData(eset)$se.exprs
	nas<-sum(is.na(sd))
	message("HAVE TO CHANGE THE MESSAGE WITH THE CORRECT ACCESSION MODIFICATION LABELS")
	message("*** checking standard deviation in assayData([eset])$se.exprs:")
	if (nas>0){
		message("   E) some standard deviation are NA in assayData([eset])$se.exprs")
		sdOK<-FALSE
	}
	else{
		nonpos<-sum(sd<0)
		if (nonpos>0){ 
			message("   E) some standard deviation you have given in [eset]@se.exprs are negative.")
			sdOK<-FALSE
		}
	}
	if(sdOK) message("					OK")
	else{ 
		message("      There seems to be a problem with the standard deviation")
		message("      of the measurement errors you have given in assayData([eset])$se.exprs.")
		message("      Read the package vignette (section 2.1)")
	}
	everythingOK<-nopdwarnings*timevecsOK*expnamestest*rngOK*sdOK
	message("  ")
	if(everythingOK) message("everything seems to be OK")
	else{
		 message("*** Summary ***")
		 message("There seems to be a problem with:")
		 if (!nopdwarnings) message("   A) your phenoData")
		 if (!timevecsOK) message("   B) the time vectors")
		 if (!expnamestest) message("   C) some columns names in the phenoData that do not correspond to column names in the expression set")
		 if (!rngOK) message("   D) the dynamic range")
		 if (!sdOK) message("   E) the standard deviation that is given in [eset]@se.exprs")
	} 
}

estimerrors<-function(eset,plattid,refchips,errtable){
#this function takes an eset as input and returns 
#a completed eset (with standard deviation of measurement error estimates) in one of the slots
	if(missing(eset)){#if the expression set is not given as an input the function just returns the list of supportd plattorms
		res<-.returnsupportedplattforms()
		res
	}
	else{
		if(class(eset)=="ExpressionSet"){#sanity check
			if(missing(errtable)){
				if(missing(plattid)){#a plattform identity has not been given 
					foundplattid<-.plattformmatch(eset=eset)
					if(foundplattid=="notfound"){#the search using genes ids did not yield anything
						message("could not find a suitable plattform, measurement error not estimated")
						res<-eset
						res
					}
					else{
						errtable<-rHVDMplattforms[[foundplattid]]$errtable
						res<-.computerrs(eset=eset,errtable=errtable,refs=refchips)
					}
				}
				else{#a plattform id has been entered by the user
					foundplattid<-.returnplattidnumber(plattid=plattid)
					if(foundplattid=="notfound"){#plattform identity not found
						message(paste("the plattform",plattid,"is not supported by rHVDM"))
						message("list of plattform supported by rHVDM:")
						res<-.returnsupportedplattforms()
						print(res)
						res<-eset
						res
					}
					else{#all is well
						errtable<-rHVDMplattforms[[foundplattid]]$errtable
						res<-.computerrs(eset=eset,errtable=errtable,refs=refchips)
						res
					}
				}
			}
			else{#an error table has been input by the user and this overrides the rest
				res<-.computerrs(eset=eset,errtable=errtable,refs=refchips)
			}
		}
		else{
			message("the object input as eset is not an ExpressionSet")
			res<-eset
			res
		}
	}
}