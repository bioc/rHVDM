#user function(s) for rHVDM (should be OK for all types of models)

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
	colnames<-dimnames(eset@exprs)[[2]]
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
	rng<-range(eset@exprs)
	rngOK<-TRUE
	if(rng[2]-rng[1]<20){
		message("   D) your data seems to be log-transformed, make sure that is it not")
		rngOK<-FALSE
	}
	else message("					OK")
	#checking the standard deviation
	sdOK<-TRUE
	sd<-eset@se.exprs
	nas<-sum(is.na(sd))
	message("*** checking standard deviation in [eset]@se.exprs:")
	if (nas>0){
		message("   E) some standard deviation are NA in [eset]@se.exprs")
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
		message("      of the measurement errors you have given in eset@se.exprs.")
		message("      Try (for example) [eset]@se.exprs<-[eset]@exprs*0.1+5.0")
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