
# Equation from Hussey et al. (2013, pg 241) (with correction published May2nd 2014)



#Corrected trophic positions


 Nlim<--5.92/-0.27
 
 k<--log((5.92 - Nlim)/-Nlim)

 TPbase<-2.5

hussey_tp<-function(Nbase=9, NTP){



 TP<-(log(Nlim - Nbase) - log(Nlim - NTP))/k + TPbase

return(TP)

}

# Original trophic position - additive fractionation

orig_tp<-function(Nbase=9, NTP) {


	TP_orig<-TPbase + (NTP - Nbase)/3.4

	return(TP_orig)
}


## PPMR estimate from slope of regression between trophic position ~ body mass class

ppmr<-function(logbase=2, b){

	ppmr<-logbase^(1/b)
	return(ppmr)

}



#### Simulation data. How does PPMR change with the new method?

## Need range of Nbase values, and randomly generated N15 values for set body mass classes
## Using log2 body mass class, as in Jennings et al. 2002


# function that creates a dataframe for all N15 and MASS estimates


nitrogen_data_func<-function(Nbase=c(2:14), mass=c(2:12),  TPbase=2, n=1){

## parameters for Husseys TP correction equation

nitrogen_dat<-matrix(nrow=length(Nbase)*length(mass), ncol=6)
nitrogen_dat[,1]<-rep(mass, times=length(Nbase))

nitrogen_dat[,3]<-rep(Nbase, each=length(mass))

i<-1

repeat{

j<-1


	start<-j + i*length(mass) - j*length(mass)
	finish<-j*length(mass)  + (i*length(mass) - 1) - (j*length(mass)-1)

	data<-nitrogen_dat[start:finish,]

#	data[j,2]<-runif(1,0,n) + data[j,3]

data[j,2]<-n + data[j,3] +6.5

repeat{   

	
	
	j<-j+1
	
	#data[j,2]<-runif(1, 0, n)+data[j-1,2]

	data[j,2]<-n+data[j-1,2]
	
	if(j==length(mass)) {break}
}

	data[,4]<-orig_tp(Nbase=data[1,3], NTP=data[,2])
	data[,5]<-hussey_tp(Nbase=data[1,3], NTP=data[,2])
	data[,6]<-i

	nitrogen_dat[start:finish,]<-data
	
	i<-i+1

	if(finish==dim(nitrogen_dat)[1]) {break}

}

colnames(nitrogen_dat)<-c("MASS", "N15", "Nbase", "orig_tp", "cor_tp", "community")
return(nitrogen_dat)

}
