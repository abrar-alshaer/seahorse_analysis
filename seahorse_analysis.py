#Author: Abrar Al-Shaer, Shaikh Lab
#Date: 12/11/2018
#This program analyzes the output of seahose data from a mitochondrial stress test run.

#!/usr/bin/env python

#Initializing imported packages
import pandas as pd
from pandas import DataFrame
import numpy as np
import math
from scipy import stats
import os, sys

#Data must be exported from Seahorse Wave software into GraphPad Prism file. Be sure to deselect any wells that you want to exclude from the analysis.
#GraphPad file provides you with averages of each well for each timepoint
#Copy GraphPad OCR and ECAR data into excel file.
##PRE-PROCESSING: If you wish to exclude entire columns (wells) from the analysis be sure to remove them. DO NOT remove individual row values.
##If you wish to exclude any individual row value from the analysis, replace the value with "NaN". -> Python will skip NaN values and exclude them from analysis.
##To analyze whether a technical replicate is an outlier, run the seahorse_outlier_analysis.R program, or calculate outliers with the 1.5*IQR method.
##If each column does not contain same number of rows, the program will fail.
##Make sure to label each column with the corresponding biological replicate (ex: B cell Mouse 1)
##ALL BIOLOGICAL REPLICATES MUST HAVE THE SAME PREFIX (ex: B cell Mouse 1, B cell Mouse 2, etc...)
##Transpose the file and save it as a comma delimmited file (CSV)

#FILE INPUTS
#OCR & ECAR files MUST have identically ordered row and column identifiers (rows are samples, columns are timepoints) - or else program will output incorrect results
seaFileOCR = 'seahorseOCR_Day21_COMBO_Data_RMoutliers.csv' #make sure your file is transposed so that the timepoints are the column headers and groups are the rows
seaFileECAR = 'SeahorseECAR_Data_COMBO_DAY21.csv'

#specify biological REPLICATES (overall groups)
#The group names must be unique from each other, from example, you cannot have B cell HF and B cell HF_SPM because both contain the same prefix 'B cell HF'
#That is why you should change it to B cell HF and B cell SPM_HF to make sure all prefixes are unique. Or else the program will assign that sample to the wrong group.
#Day0 & Day21 bioRep
#'B cell Lean', 'B cell HF', 'B cell SPM_HF', 'CD4 Lean', 'CD4 HF', 'CD4 SPM_HF', 'CD8 Lean', 'CD8 HF', 'CD8 SPM_HF'
bioReps = ['B cell Lean', 'B cell HF', 'B cell SPM_HF', 'CD4 Lean', 'CD4 HF', 'CD4 SPM_HF', 'CD8 Lean', 'CD8 HF', 'CD8 SPM_HF']
#Day10 11/17/18 bioRep
#'B cell Lean', 'B cell HF', 'B cell SPM_HF', 'CD4 Lean', 'CD4 HF', 'CD4 SPM_HF', 'CD8 Lean', 'CD8 HF', 'CD8 SPM_HF', 'Lean Lung', 'HF Lung', 'HF_SPM Lung'
#bioReps = ['B cell Lean', 'B cell HF', 'B cell SPM_HF', 'CD4 Lean', 'CD4 HF', 'CD4 SPM_HF', 'CD8 Lean', 'CD8 HF', 'CD8 SPM_HF', 'Lean Lung', 'HF_SPM Lung'] #make sure the spelling matches your input data group names (from your CSV files)

from_csv = pd.read_csv(seaFileOCR, sep=None, header=0, engine='python') #reading in the OCR CSV file
#Checking to make sure that negative OCR values were filtered out
if (from_csv < 0).any().sum() > 0:
    print "\nYour dataset contains negative OCR values. Please filter them out by either removing the sample (row) completely or replacing that negative value with NaN\n"
    sys.exit(0)

print "\n-------------AVERAGED OCR TECHNICAL REPLICATES-------------\n"
avg_tech_reps_OCR = from_csv.groupby('Time (minutes)').mean().reset_index()
print avg_tech_reps_OCR #we now have biological replicates

sample = []
group = []
nonMitoResp = []
basalOCR = []
basalOCR_nonMitoResp = []
maxResp = []
protonLeak = []
ATP = []
spareRespCap = []
spareRespCapPercent = []
#acuteResponse = []
couplingEff = []
for index, row in avg_tech_reps_OCR.iterrows():
    sample.append(row[0])
    for i in bioReps:
        if i in row[0]:
            group.append(i) #assign the sample to a group (helpful for later calculating mean/standard dev for biological replicates)
    last = [row[len(row)-1], row[len(row)-2], row[len(row)-3]] #adding last 3 time measurements to a list
    median_nmr = np.nanmedian(last) #median non-mitochondrial respiration #ingore NaNs
    oligo_min = np.nanmin([row[4], row[5], row[6]]) #finding minimum measure after oligomycin injection #ingore NaNs
    nonMitoResp.append(median_nmr) #calculating non mitochondrial respiration from median of last 3 time measurements
    basalOCR.append(row[3]) #basal OCR
    basalOCR_nonMitoResp.append(row[3]-median_nmr) #basal OCR - non-mitochondrial respiration (from median calculated above)
    max_fccp = np.nanmax([row[7],row[8],row[9]]) #getting max fccp measurement #ignore NaNs
    max_respiration = max_fccp-median_nmr #Maximal respiration = max measurement after fccp injection - non-mitochondrial respiration
    maxResp.append(max_respiration) #appending max respiration to a list
    protonLeak.append(oligo_min-median_nmr) #proton leak = minimum rate after oligo injection - non-mitochondrial respiration
    atp_product = row[3]-oligo_min #ATP production = basal OCR - minimum oligo injection
    ATP.append(atp_product) #append to ATP list
    spareRespCap.append(max_respiration-row[3]) #spare respiratory capacity = max respiration - basal respiration (basal OCR)
    spareRespCapPercent.append((max_respiration/row[3])*100) #percentage for spare respiratory capacity taking basal OCR into account as a ratio
    #acuteResponse.append(row[3]-??)
    couplingEff.append((atp_product/row[3])*100) #Coupling efficiency = ATP production rate/basal OCR * 100 ##how much ATP is produced relative to oxygen consumed

from_csv2 = pd.read_csv(seaFileECAR, sep=None, header=0, engine='python') #reading in the ECAR CSV file
#Checking to make sure that negative ECAR values were filtered out
if (from_csv < 0).any().sum() > 0:
    print "\nYour dataset contains negative ECAR values. Please filter them out by either removing the sample (row) completely or replacing that negative value with NaN\n"
    sys.exit(0)

print "\n-------------AVERAGED ECAR TECHNICAL REPLICATES-------------\n"
avg_tech_reps_ECAR = from_csv2.groupby('Time (minutes)').mean().reset_index()
print avg_tech_reps_ECAR #we now have biological replicates

basalECAR = []
for index, row in avg_tech_reps_ECAR.iterrows():
    basalECAR.append(row[3]) #basal ECAR

OCR_ECAR_ratio = []
for i in range(0, len(basalOCR)):
    OCR_ECAR_ratio.append(basalOCR[i]/basalECAR[i]) #to determine if cell is oxidative or glycolytic calculate OCR/ECAR ratio

print "\n---------Seahorse Calculations Output File---------\n"
#creating calculations dataframe
print "Sample:", len(sample), "Group:", len(group), "Non mito", len(nonMitoResp), "ocr", len(basalOCR), "ecar", len(basalECAR), "ocr-nonMitoResp", len(basalOCR_nonMitoResp), "max resp", len(maxResp), "proton", len(protonLeak), "ATP", len(ATP), "spare resp", len(spareRespCap), "Spare %", len(spareRespCapPercent), "coupling", len(couplingEff), "ocr/ecar", len(OCR_ECAR_ratio)
print "SAMPLE\n", sample
print "GROUPS\n", group
df = pd.DataFrame({'1) Samples': sample, '2) Group': group, 'Non-Mitochondrial Respiraion': nonMitoResp, 'Basal OCR': basalOCR, 'Basal OCR - Non-mito respiration': basalOCR_nonMitoResp, 'Maximum Respiration': maxResp, 'Proton Leak': protonLeak, 'ATP Production': ATP, 'Spare Respiratory Capacity': spareRespCap, 'Spare Respiratory Capacity %': spareRespCapPercent, 'Coupling Efficiency': couplingEff, 'Basal ECAR': basalECAR, 'OCR/ECAR Ratio': OCR_ECAR_ratio})
print df
df.to_csv('seahorse_analysisResults_Day21_COMBO_Data_12.11.18.csv') #send dataframe to a file

print "\n\n---------Seahorse Calculations Merged Biological Replicates---------\n"

avg_bio_reps = df.groupby('2) Group').mean().reset_index() #averaging biological replicates
std_bio_reps = df.groupby('2) Group').std().reset_index() #standard deviation of biological replicates
sem_bio_reps = (df.drop(columns=['1) Samples'])).groupby('2) Group').sem().reset_index() #standard error of biological replicates

print "\n Mean of Biological Replicates\n",avg_bio_reps
print "\n Standard Deviation of Biological Replicates\n", std_bio_reps
print "\n SEMs of Biological Replicates\n", sem_bio_reps

avg_bio_reps.to_csv('seahorse_analysisResults_Means_Day21_COMBO_Data_12.11.18.csv') #send means dataframe to a file
std_bio_reps.to_csv('seahorse_analysisResults_StdDev_Day21_COMBO_Data_12.11.18.csv') #send standard deviations ataframe to a file
sem_bio_reps.to_csv('seahorse_analysisResults_SEMs_Day21_COMBO_Data_12.11.18.csv') #send SEMs dataframe to a file
