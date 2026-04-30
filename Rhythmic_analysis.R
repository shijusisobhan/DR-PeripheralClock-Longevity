rm(list=ls())

# required library
library(rain)

#setwd('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/GitHub_folder/DR-PeripheralClock-Longevity/Data')

# Load your batch corrected TPM values
data_HC_LC<-read.csv(".../Data/Batch_corrected_TPM.csv")

# data_HC_LC<-read.csv("Batch_corrected_TPM.csv")

#-------------------------------------------------------------------------
# function for calculate fold change
# log2FC based on TPM data
#-----------------------------------------------------------------------

Max_min_amp<-function(x){
  
  if(min(x)==0){
    amp=max(x)
  }else
    
  { amp<- max(x)/min(x) }
  
  
}

#---------------------------------------------------------------------------

# subset control data only
data_HC<-data_HC_LC[,c(1,2:25)]
#--------------------------------------------------------------------------------

data=data_HC[-1] # data only
Gene<-data_HC[1] # Gene name

#--------------------------------------------------------------------
# Rearrange data which is in suitable for Rain 
# ZT0_1, ZT0_2, ZT2_1, ZT2_2........ZT22_1, ZT22_2
#------------------------------------------------------------------

# Original column indices
first_half <- 1:12
second_half <- 13:24
# Interleave the two halves
new_order <- as.vector(rbind(first_half, second_half))
data<-data[,new_order]

#-----------------------------------------------------------------------------------------------------------
# Rain analysis

menet.ossc <- rain(t(data), deltat = 2, period = 24,
                   nr.series = 2, peak.border = c(0.3, 0.7), verbose=TRUE)

##  Correct p-values for multiple testing using Benjamini-Hochberg and view the top genes
# the `transform` function lets you add new columns to the data.frame from the already present columns.
menet.ossc <- transform(menet.ossc, 
                        Q.Val = p.adjust(pVal, method = "BH"))

amplitude_data = apply(data , 1, Max_min_amp) # Find the max/min
log2_FC<-log2(amplitude_data)

Rain_statics_HC<-cbind(Gene,menet.ossc, log2_FC)

Rain_statics_sig_HC<-Rain_statics_HC[which(Rain_statics_HC$Q.Val<0.1 & Rain_statics_HC$log2_FC>0.6),]

# save the Rhythmic analysis results of control data
#write.csv(Rain_statics_HC,'Rain_results_HC.csv')

#----------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------

# subset DR data only
data_LC<-data_HC_LC[,c(1,26:49)]
#--------------------------------------------------------------------------------

data=data_LC[-1] # data only
Gene<-data_LC[1] # Gene name

#--------------------------------------------------------------------
# Rearrange data which is in suitable for Rain 
# ZT0_1, ZT0_2, ZT2_1, ZT2_2........ZT22_1, ZT22_2
#------------------------------------------------------------------

# Original column indices
first_half <- 1:12
second_half <- 13:24
# Interleave the two halves
new_order <- as.vector(rbind(first_half, second_half))
data<-data[,new_order]

#-----------------------------------------------------------------------------------------------------------
# Rain analysis

menet.ossc <- rain(t(data), deltat = 2, period = 24,
                   nr.series = 2, peak.border = c(0.3, 0.7), verbose=TRUE)

##  Correct p-values for multiple testing using Benjamini-Hochberg and view the top genes
# the `transform` function lets you add new columns to the data.frame from the already present columns.
menet.ossc <- transform(menet.ossc, 
                        Q.Val = p.adjust(pVal, method = "BH"))

amplitude_data = apply(data , 1, Max_min_amp) # Find the max/min
log2_FC<-log2(amplitude_data)

Rain_statics_LC<-cbind(Gene,menet.ossc, log2_FC)

Rain_statics_sig_LC<-Rain_statics_LC[which(Rain_statics_LC$Q.Val<0.1 & Rain_statics_LC$log2_FC>0.6),]

# save the Rhythmic analysis results of control data
#write.csv(Rain_statics_LC,'Rain_results_LC.csv')

#----------------------------------------------------------------------------------------------------





