
### Script accompanying Hertz et al. (2014) Estimation of predator-prey mass ratios using stable isotopes: sources of errors.
##  https://github.com/baumlab/ppmr-isotopes
##  
## Script written by JPW Robinson (June 2014)
## Script tested by E Hertz (October 2014)


#### Function details

### hussey_tp : scaled trophic position estimates. From Hussey et al. (2013, Ecology Letters, pg. 241. Correction published May2nd 2014)

### Arguments: NTP = d15N estimates (default = 9); Nbase = baseline 15N; TPbase = baseline trophic position; Nlim and k from Hussey et al. (2013)
### Returns: scaled trophic positions

 Nlim<--5.92/-0.27
 
 k<--log((5.92 - Nlim)/-Nlim)


hussey_tp<-function(Nbase=9, NTP, TPbase=2.5){


 TP<-(log(Nlim - Nbase) - log(Nlim - NTP))/k + TPbase

return(TP)

}

# orig_tp - additive fractionation - trophic position estimates.

## Arguments: Nbase = baseline 15N (default = 9); NTP = d15N estimates; TPbase = baseline trophic position
## Returns: additive trophic positions

orig_tp<-function(Nbase=9, NTP, TPbase=2) {


	TP_orig<-TPbase + (NTP - Nbase)/3.4

	return(TP_orig)
}


## ppmr - predator-prey mass ratio estimate from slope of regression between trophic position ~ body mass class
## From Jennings et al. (2002, Marine Ecology Progress Series, 240:11-20)
## Arguments: logbase = log base for mass class bins in trophic level ~ body mass regression; b = slope of regression between trophic position ~ body mass class

ppmr<-function(logbase=2, b){

	ppmr<-logbase^(1/b)
	return(ppmr)

}


#########################################################

############# PART II ###################################


#### Simulate d15N values for range of body mass classes.

## Need range of Nbase values, and randomly generated N15 values for set body mass classes
## Using log2 body mass class, as in Jennings et al. 2002


# nitrogen_data_func - function that creates a dataframe for all N15 and MASS estimates

## Arguments - Nbase = baseline N values; mass = range of body masses; TPbase = baseline trophic position
##           - n = change in d15N between mass classes (default = 1)

nitrogen_data_func<-function(Nbase=c(2:10), mass=c(2:12),TPbase=2, n=1){

## Set up Nbase and mass values.

nitrogen_dat<-matrix(nrow=length(Nbase)*length(mass), ncol=6)
nitrogen_dat[,1]<-rep(mass, times=length(Nbase))

nitrogen_dat[,3]<-rep(Nbase, each=length(mass))

i<-1

repeat{

j<-1


	start<-j + i*length(mass) - j*length(mass)
	finish<-j*length(mass)  + (i*length(mass) - 1) - (j*length(mass)-1)

	data<-nitrogen_dat[start:finish,]

## add n (default = 1) to each 15N mass class estimate

data[j,2]<-n + data[j,3] 

repeat{   

	
	j<-j+1
	
	data[j,2]<-n+data[j-1,2]
	
	if(j==length(mass)) {break}
}

### Estimate trophic position (additive and scaled)

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




