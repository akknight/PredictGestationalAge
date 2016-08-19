#GAPredictor.R
#This file is a slightly modified version of additional file 20 from
#Horvath Genome Biology 2013 14:R115  doi:10.1186/gb-2013-14-10-r115
#Changes made:
#1) Removed the anti.trafo function to transform adult age
#2) Supplied a new set of predictive CpGs and coefficients
#3) Minor changes to comments
#4) Minor changes to "NormalizeAndPredictGA.R", formerly Horvath's additional file 25

#Wrapper program for predicting gestational age from cord blood and blood spots
#Code modified from DNAm age calculator by Steve Horvath.
#Please see: Horvath Genome Biology 2013 14:R115   doi:10.1186/gb-2013-14-10-r115
# Note: some supplementary files must be downloaded from
#Horvath Genome Biology 2013 14:R115   doi:10.1186/gb-2013-14-10-r115
#See InstructionsForGAPredictor.docx (Additional file X) for instructions

#call required libraries

library(WGCNA)
library(sqldf)
library(impute)
library(flashClust)
library(dynamicTreeCut)
library(RPMM)

#Note: to install bioconductor packages, use:
#source("https://bioconductor.org/biocLite.R")
#biocLite("impute")


#load in supplemental file 22
probeAnnotation21kdatMethUsed=read.csv("suppl22.csv")

#load supplemental file 21

probeAnnotation27k=read.csv("suppl21.csv")

#Input the set of CpG coefficients here.

datClock=read.csv("cgprobesGApredictor.csv",as.is=T)

#Input your test dataset

dat0=read.csv("TestDataset.csv",as.is=TRUE)

#removes probes from the annotation file that were QC'ed out of your merged dataset

#The CpG identifiers in dat0 should be named CpGName

tf<-probeAnnotation21kdatMethUsed$Name%in%dat0$CpGName
length(tf)
newtf<-probeAnnotation21kdatMethUsed$Name[tf==T]
length(newtf)
komit<-probeAnnotation21kdatMethUsed[newtf,]
probeAnnotation21kdatMethUsed<-komit



#insert supplemental file 24

source("suppl24.txt")



nSamples=dim(dat0)[[2]]-1
nProbes= dim(dat0)[[1]]
# the following command may not be needed. But it is sometimes useful when you use read.csv.sql
dat0[,1]= gsub(x=dat0 [,1],pattern="\"",replacement="") 
#Create a log file which will be output into your directory
# The code looks a bit complicated because it serves to create a log file (for error checks etc). # It will automatically create a log file.
file.remove("LogFile.txt")

file.create("LogFile.txt")

DoNotProceed=FALSE  
cat(paste( "The methylation data set contains", nSamples, "samples (e.g. arrays)" ,nProbes, " probes."),file="LogFile.txt")
if (nSamples==0) {DoNotProceed=TRUE; cat(paste( "\n ERROR: There must be a data input error since there seem to be no samples.\n Make sure that you input a comma delimited file (.csv file)\n that can be read using the R command read.csv.sql . Samples correspond to columns in that file  ."), file="LogFile.txt",append=TRUE) } 
if (nProbes==0) {DoNotProceed=TRUE; cat(paste( "\n ERROR: There must be a data input error since there seem to be zero probes.\n Make sure that you input a comma delimited file (.csv file)\n that can be read using the R command read.csv.sql  CpGs correspond to rows.")   , file="LogFile.txt",append=TRUE) } 
if (  nSamples > nProbes  ) { cat(paste( "\n MAJOR WARNING: It worries me a lot that there are more samples than CpG probes.\n Make sure that probes correspond to rows and samples to columns.\n I wonder whether you want to first transpose the data and then resubmit them? In any event, I will proceed with the analysis."),file="LogFile.txt",append=TRUE) }
if (  is.numeric(dat0[,1]) ) { DoNotProceed=TRUE; cat(paste( "\n Error: The first column does not seem to contain probe identifiers (cg numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file contains probe identifiers such as cg00000292. Instead it contains ", dat0[1:3,1]  ),file="LogFile.txt",append=TRUE)  } 
if (  !is.character(dat0[,1]) ) {  cat(paste( "\n Major Warning: The first column does not seem to contain probe identifiers (cg numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file contains CpG probe identifiers such as cg00000292. Instead it contains ", dat0[1:3,1]  ),file="LogFile.txt",append=TRUE)  } 
datout=data.frame(Error=c("Input error. Please check the log file for details","Please read the instructions carefully."), Comment=c("", "email Anna Knight. anna.knight@emory.edu"))
if ( ! DoNotProceed ) 
  nonNumericColumn=rep(FALSE, dim(dat0)[[2]]-1)
for (i in 2:dim(dat0)[[2]] ){ nonNumericColumn[i-1]=! is.numeric(dat0[,i]) }
if (  sum(nonNumericColumn) >0 ) { cat(paste( "\n MAJOR WARNING: Possible input error. The following samples contain non-numeric beta values: ", colnames(dat0)[-1][ nonNumericColumn], "\n Hint: Maybe you use the wrong symbols for missing data. Make sure to code missing values as NA in the Excel file. To proceed, I will force the entries into numeric values but make sure this makes sense.\n" ),file="LogFile.txt",append=TRUE)  } 
XchromosomalCpGs=as.character(probeAnnotation27k$Name[probeAnnotation27k$Chr=="X"])
selectXchromosome=is.element(dat0[,1], XchromosomalCpGs )
selectXchromosome[is.na(selectXchromosome)]=FALSE
meanXchromosome=rep(NA, dim(dat0)[[2]]-1)
if (   sum(selectXchromosome) >=500 )  {
  meanXchromosome= as.numeric(apply( as.matrix(dat0[selectXchromosome,-1]),2,mean,na.rm=TRUE)) }
if (  sum(is.na(meanXchromosome)) >0 ) { cat(paste( "\n \n Comment: There are lots of missing values for X chromosomal probes for some of the samples. This is not a problem when it comes to estimating age but I cannot predict the gender of these samples.\n " ),file="LogFile.txt",append=TRUE)  }  
match1=match(probeAnnotation21kdatMethUsed$Name , dat0[,1])
if  ( sum( is.na(match1))>0 ) { 
  missingProbes= probeAnnotation21kdatMethUsed$Name[!is.element( probeAnnotation21kdatMethUsed$Name , dat0[,1])]    
  DoNotProceed=TRUE; cat(paste( "\n \n Input error: You forgot to include the following ", length(missingProbes), " CpG probes (or probe names):\n ", paste( missingProbes, sep="",collapse=", ")),file="LogFile.txt",append=TRUE)  } 


match1=match(probeAnnotation21kdatMethUsed$Name , dat0[,1])
if  ( sum( is.na(match1))>0 ) stop(paste(sum( is.na(match1)), "CpG probes cannot be matched"))
dat1= dat0[match1,]
asnumeric1=function(x) {as.numeric(as.character(x))}
dat1[,-1]=apply(as.matrix(dat1[,-1]),2,asnumeric1)
set.seed(1)



#Normalize the Test Dataset
normalizeData=T

#load in "Normalization.R"
source("NormalizeAndPredictGA.R")

#This Output file will contain the age prediction
write.table(datout,"Outputfile.csv", row.names=F, sep="," )
dat0UsedNormalized=data.frame(CpGName=colnames(datMethUsedNormalized), data.frame(t(datMethUsedNormalized) ))
write.table(dat0UsedNormalized,file="dat0UsedNormalized.csv",sep=",",row.names=F)

