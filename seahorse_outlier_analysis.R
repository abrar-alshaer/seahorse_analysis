#MAKE SURE TO PROVIDE THE FINAL TRANSPOSED DATA FILE - same input for python seahorse program (seahorse_analysis.py)
setwd("C:/Users/Abrar/Dropbox/UNC_OneDrive/Shaikh Lab/Seahorse/")
rm(list=ls())

#RUN THIS PROGRAM USING EITHER OCR OR ECAR FILES - AN EXAMPLE INPUT FILE "SeahorseECAR_Data_COMBO_DAY24_lean_HF_only" IS PLACED IN THE seahorse_analysis GITHUB FOLDER
#Combo file of all biological replicates (all experiments for one timepoint)
#ocr <- read.csv("Seahorse_Day10_ComboData_Waves1-3.csv", header = TRUE) #load in OCR data 

#SeahorseOCR_Data_COMBO_DAY21_lean_HF_only
#SeahorseOCR_Data_COMBO_DAY10_lean_HF_only & SeahorseECAR_Data_COMBO_DAY10_lean_HF_only
ocr <- read.csv("SeahorseECAR_Data_COMBO_DAY24_lean_HF_only.csv", header = TRUE) #load in OCR data 
file_name <- "SeahorseECAR_Data_COMBO_DAY24_lean_HF_only"

#####OUTLIER ANALYSIS

groups = levels(ocr$Group) #create a variable called groups for all levels (row groups) of the dataframe
columns = colnames(ocr)[3:length(colnames(ocr))] #store columns in variable

datalist_ocr = list()
outliers_data1 = list()
c = 1

#loop through each group (rows) for COMBO data
for(i in 1:length(groups)) {
  print(groups[i]) #print group
  sub_ocr <- subset(ocr, ocr$Group %in% groups[i]) #subset dataframe by the biological group
  #loop through each column within each group
  for(j in 1:length(columns)){
    print(paste("COLUMN:",columns[j]))
    out1 = quantile(sub_ocr[,columns[j]], c(1)/4, na.rm = TRUE)[[1]]-1.5*IQR(sub_ocr[,columns[j]], na.rm = TRUE)
    out1 = floor(out1*10)/10 #needed due to rounding error (this statement prevents rounding of the outlier cutoff so we don't detect very small deviations from the IQR)
    out2 = quantile(sub_ocr[,columns[j]], c(3)/4, na.rm = TRUE)[[1]]+1.5*IQR(sub_ocr[,columns[j]], na.rm = TRUE)
    #out2 = floor(out2*10)/10 #not necessary because it reduces the threshold for a high outlier (makes it so that we are too sensitive and detect very small deviations from the IQR)
    n = 1
    outliers_data1 = list()
    for(k in sub_ocr[,columns[j]]){
      if(is.nan(k) == FALSE && (k < out1 | k > out2)){
        outliers_data1[[n]] = k
        print(columns[j])
        print(k)
        n = n + 1
      }
    }
    #if there are outliers (list isn't empty)
    if(length(outliers_data1)>0)
    {
      datalist_ocr[[c]] <- c(groups[i],columns[j],outliers_data1)
      c = c+1
      #outliers_data1 = list()
    }
  }
}

#big_data_ocr = do.call(rbind, datalist_ocr) #combines all previous dataframes from for loop
#list_data = data.frame(datalist = unlist(datalist_ocr))

#If big_data_data1 from rbind function fails
matrix_data = as.matrix(datalist_ocr)
write.csv(matrix_data, paste(file_name,"Outliers_7.10.19.csv"))

#Check to see if the list is empty (no outliers) - if not, add column headers & send dataframe to a file
# if(length(big_data_ocr) < 1){
#   print("There are no OCR outliers")
# } else {
#   colnames(big_data_ocr) <- c("Sample", "Measurement", "Value")
#   write.csv(big_data_ocr, "seahorseOCR_outliers.csv")
# }

#write.csv(matrix_data, "seahorseOCR_outliers_11.10.18.csv")

#######replacing outliers with NaN

#extracting all data1 values from the outliers list and placing them in a new list (removing first 2 elements for group & timepoint)
#for the second list it's appending only element 2 - which is the column the outlier belongs to - useful for the next for loop when replacing outliers with NaNs
datalist_ocr_new <- list()
datalist_ocr_newNames <- list()
count <- 1
for (i in 1:length(datalist_ocr)){
  datalist_ocr_new[[count]] <- datalist_ocr[[i]][c(3:length(datalist_ocr[[i]]))]
  datalist_ocr_newNames[[count]] <- datalist_ocr[[i]][2] #list of column names that contained outliers
  count = count + 1
}

#converting the list of lists into a single one level list
datalist_ocr_new <- unlist(datalist_ocr_new)
datalist_ocr_newNames <- unlist(datalist_ocr_newNames)

#replacing the outliers step
#only loop through the list of outlier columns
for (row in datalist_ocr_newNames) { 
  #print(row)
  #print(ocr[[row]])
  for (i in 1:length(ocr[[row]])){
    #print(ocr[[row]][i])
    for (j in 1:length(datalist_ocr_new))
      if(is.nan(ocr[[row]][i]) == FALSE && (ocr[[row]][i] == datalist_ocr_new[j])) {
        ocr[[row]][i] <- "NaN" }
  }
}

#new dataframe without outliers
write.csv(ocr, paste(file_name,"_RMoutliers_7.10.19.csv"))


#calculate outlier values based on the IQR formula MANUALLY
#quantile(sub$X1.413664)[2]-1.5*IQR(sub$X1.413664) Q1-1.5*IQR
#quantile(sub$X1.413664)[4]+1.5*IQR(sub$X1.413664) Q3+1.5*IQR
