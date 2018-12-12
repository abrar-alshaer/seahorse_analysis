setwd("C:/Users/Abrar/Dropbox/UNC_OneDrive/Shaikh Lab/Seahorse/")
rm(list=ls())

#MAKE SURE TO PROVIDE THE FINAL TRANSPOSED DATA FILE - same input for python seahorse program (seahorse_analysis.py)
ocr <- read.csv("SeahorseOCR_Data_COMBO_DAY21.csv", header = TRUE) #load in OCR data 
ecar <- read.csv("SeahorseECAR_Data_COMBO_DAY21.csv", header = TRUE) #load in ECAR data

#####OCR ANALYSIS

groups = levels(ocr$Time..minutes.) #create a variable called groups for all levels (row groups) of the dataframe
columns = colnames(ocr)[2:length(colnames(ocr))] #store columns in variable

datalist_ocr = list()
c = 1

#loop through each group (rows) for OCR data
for(i in 1:length(groups)) {
  print(groups[i]) #print group
  sub_ocr <- subset(ocr, ocr$Time..minutes. %in% groups[i]) #subset dataframe by the biological group
  #loop through each column within each group
  for(j in 1:length(columns)){
    #print(paste("COLUMN (Time Point):",columns[j]))
    outliers_ocr <- boxplot.stats(sub_ocr[,columns[j]])$out #output of outlier values to a list, calculate outlier values based on the IQR formula (Q1-1.5*IQR, Q3+1.5*IQR)
    #if there are outliers (list isn't empty)
    if(length(outliers_ocr)>0)
      {
        datalist_ocr[[c]] <- c(groups[i],columns[j],boxplot.stats(sub_ocr[,columns[j]])$out)
        c = c+1
      }
  }
}

#If the matrix is uneven the command below won't work - then just run matrix_data = as.matrix(datalist_ocr) below
big_data_ocr = do.call(rbind, datalist_ocr) #combines all previous dataframes from for loop
list_data = data.frame(datalist = unlist(datalist_ocr))
matrix_data = as.matrix(datalist_ocr)
#Check to see if the list is empty (no outliers) - if not, add column headers & send dataframe to a file
if(length(big_data_ocr) < 1){
  print("There are no OCR outliers")
} else {
  colnames(big_data_ocr) <- c("Sample", "Timepoint", "OCR Value")
  write.csv(big_data_ocr, "seahorseOCR_outliers.csv")
}

write.csv(matrix_data, "seahorseOCR_outliers_Day21_COMBO_Data.csv")

#######replacing outliers with NaN

#extracting all OCR values from the outliers list and placing them in a new list (removing first 2 elements for group & timepoint)
datalist_ocr_new <- list()
count <- 1
for (i in 1:length(datalist_ocr)){
  datalist_ocr_new[[count]] <- datalist_ocr[[i]][c(3:length(datalist_ocr[[i]]))]
  count = count + 1
}

#converting the list of lists into a single one level list
datalist_ocr_new <- unlist(datalist_ocr_new)

#replacing the outliers step
for (row in 1:length(ocr)) {
  #print(row)
  #print(ocr[[row]])
  for (i in 1:length(ocr[[row]])){
    #print(ocr[[row]][i])
    for (j in 1:length(datalist_ocr_new))
      if(ocr[[row]][i] == datalist_ocr_new[j]) {
        ocr[[row]][i] <- "NaN" }
  }
}

#new dataframe without outliers
write.csv(ocr, "seahorseOCR_Day21_COMBO_Data_RMoutliers.csv")

#####ECAR ANALYSIS

groups = levels(ecar$Time..minutes.) #create a variable called groups for all levels (row groups) of the dataframe
columns = colnames(ecar)[2:length(colnames(ecar))] #store columns in variable

datalist_ecar = list()
c = 1

#loop through each group (rows) for OCR data
for(i in 1:length(groups)) {
  print(groups[i]) #print group
  sub_ecar <- subset(ecar, ecar$Time..minutes. %in% groups[i]) 
  #loop through each column within each group
  for(j in 1:length(columns)){
    outliers_ecar <- boxplot.stats(sub_ocr[,columns[j]])$out
    #if there are outliers (list isn't empty)
    if(length(outliers_ecar)>0)
    {
      datalist_ecar[[c]] <- c(groups[i],columns[j],boxplot.stats(sub_ecar[,columns[j]])$out)
      c = c+1
    }
  }
}

big_data_ecar = do.call(rbind, datalist_ecar) #combines all previous dataframes from for loop
#Check to see if the list is empty (no outliers) - if not, add column headers & send dataframe to a file
if(length(big_data_ecar) < 1){
  print("There are no ECAR outliers")
} else {
  colnames(big_data_ecar) <- c("Sample", "Timepoint", "ECAR Value")
  write.csv(big_data_ecar, "seahorseECAR_outliers.csv") 
}

#calculate outlier values based on the IQR formula MANUALLY
#quantile(sub$X1.413664)[2]-1.5*IQR(sub$X1.413664) Q1-1.5*IQR
#quantile(sub$X1.413664)[4]+1.5*IQR(sub$X1.413664) Q3+1.5*IQR