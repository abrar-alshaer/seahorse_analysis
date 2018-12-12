setwd("C:/Users/Abrar/Dropbox/UNC_OneDrive/Shaikh Lab/Seahorse/DHA PR8 Study Final Results/")
rm(list=ls())
library(ggplot2)
library(gplots)
library("dplyr")
library("ggpubr")
library(tidyverse)

sea_data <- read.csv("seahorse_analysisResults_Day21_COMBO_Data_12.11.18.csv", header = TRUE) #load in data with biological replicates

#Reordering dataframe to the groups we want to see on the x-axis of the graph
#Days 0 & 21 reference
reference = c('B cell Lean', 'B cell HF', 'B cell SPM_HF', 'CD4 Lean', 'CD4 HF', 'CD4 SPM_HF', 'CD8 Lean', 'CD8 HF', 'CD8 SPM_HF')
#Day 10 reference
#reference = c('B cell Lean', 'B cell HF', 'B cell SPM_HF', 'CD4 Lean', 'CD4 HF', 'CD4 SPM_HF', 'CD8 Lean', 'CD8 HF', 'CD8 SPM_HF', 'Lean Lung', 'HF Lung', 'HF_SPM Lung')
sea_data <- sea_data[order(sapply(sea_data$X2..Group, function(X2..Group) which(X2..Group == reference))), ] #ordering dataset based on reference variable
sea_data <- sea_data[,2:length(sea_data)]
sea_data$X2..Group <- factor(sea_data$X2..Group, levels = reference,ordered = TRUE)

#filter by cell type
sea_data <- dplyr::filter(sea_data, grepl("CD8", X2..Group)) #you can do !grepl for not "B cell" (or any other cell type)
#sea_data <- dplyr::filter(sea_data, !grepl("Lung", X2..Group)) 

####Stats & Graphing for Basal OCR####
shapiro.test(sea_data$Basal.OCR)
hist(sea_data$Basal.OCR)
qqnorm(sea_data$Basal.OCR)
qqline(sea_data$Basal.OCR)

kruskal.test(Basal.OCR ~ X2..Group,  data = sea_data)
compare_means(Basal.OCR ~ X2..Group,  data = sea_data, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 

#palette = c("lightskyblue", "lightblue2", "lightblue3", "lightblue4")
my_comparisons <- list(c("CD4 HF", "CD4 SPM_HF"))
p1 <- ggbarplot(sea_data, "X2..Group", y = "Basal.OCR", xlab = "Groups", ylab = "Basal OCR", label = FALSE, width = 0.5, add = c("mean_se"), palette = "jco", color = "X2..Group", title = "Basal OCR", legend = "None")+
  #ggrepel::geom_text_repel(aes(label = X1..Samples))+ #for labeling the graph with IDs/mice
  theme(axis.text=element_text(size=12))+ #size=8
  stat_compare_means(comparisons = my_comparisons, tip.length = c(0.01), method = "t.test", label = "p.format", hide.ns = TRUE, label.y = c(50)) #p.format to place actual p-values (uses wilcox.test)

#add = c("mean_se", "point") add = c("mean_se", "jitter")

####Stats & Graphing for Basal ECAR####
shapiro.test(sea_data$Basal.ECAR)
hist(sea_data$Basal.ECAR)
qqnorm(sea_data$Basal.ECAR)
qqline(sea_data$Basal.ECAR)

kruskal.test(Basal.ECAR ~ X2..Group,  data = sea_data)
compare_means(Basal.ECAR ~ X2..Group,  data = sea_data, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 

my_comparisons <- list(c("Lean Lung", "HF_SPM Lung"))
p2 <- ggbarplot(sea_data, "X2..Group", y = "Basal.ECAR", xlab = "Groups", ylab = "Basal ECAR", label = FALSE, width = 0.5, add = c("mean_se"), palette = "jco", color = "X2..Group", title = "Basal ECAR", legend = "None")+
  theme(axis.text=element_text(size=12))+
  #stat_compare_means(data = sea_data, method = "kruskal.test", label.y = 8.5, label.x = 1)+
  stat_compare_means(comparisons = my_comparisons, tip.length = c(0.05), method = "t.test", label = "p.format", hide.ns = TRUE, label.y = c(7.5)) #p.format to place actual p-values (uses wilcox.test)

####Stats & Graphing for ATP Production####
#must use sea_data file for calculating statistics so that it takes the biological replicates into account (sample size n)
shapiro.test(sea_data$ATP.Production)

kruskal.test(ATP.Production ~ X2..Group,  data = sea_data)
compare_means(ATP.Production ~ X2..Group,  data = sea_data, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 

my_comparisons <- list(c("CD8 Lean", "CD8 SPM_HF"))
p3 <- ggbarplot(sea_data, "X2..Group", y = "ATP.Production", xlab = "Groups", ylab = "ATP Production", label = FALSE, width = 0.5, add = c("mean_se"), palette = "jco", color = "X2..Group", title = "ATP Production", legend = "None")+
  theme(axis.text=element_text(size=12))+
  #stat_compare_means(data = sea_data, method = "kruskal.test", label.y = 52, label.x = 1)
  stat_compare_means(comparisons = my_comparisons, method = "t.test", tip.length = c(0.05), label = "p.format", hide.ns = TRUE, label.y = c(47)) #p.format to place actual p-values (uses wilcox.test)

####Stats & Graphing for Coupling Efficiency####
shapiro.test(sea_data$Coupling.Efficiency)

kruskal.test(Coupling.Efficiency ~ X2..Group,  data = sea_data)
compare_means(Coupling.Efficiency ~ X2..Group,  data = sea_data, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 

#my_comparisons <- list(c("HF B cell", "Lean CD4 cell"), c("HF B cell", "Lean CD8 cell"))
p4 <- ggbarplot(sea_data, "X2..Group", y = "Coupling.Efficiency", xlab = "Groups", ylab = "Coupling Efficiency", label = FALSE, width = 0.5, add = c("mean_se"), palette = "jco", color = "X2..Group", title = "Coupling Efficiency", legend = "None")+
  theme(axis.text=element_text(size=12))
  #stat_compare_means(data = sea_data, method = "kruskal.test", label.y = 75, label.x = 1)
#  stat_compare_means(comparisons = my_comparisons, tip.length = c(0.01), label = "p.format", hide.ns = TRUE, label.y = c(80, 85)) #p.format to place actual p-values (uses wilcox.test)

####Stats & Graphing for Max Respiration####
shapiro.test(sea_data$Maximum.Respiration)

kruskal.test(Maximum.Respiration ~ X2..Group,  data = sea_data)
compare_means(Maximum.Respiration ~ X2..Group,  data = sea_data, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 

my_comparisons <- list(c("Lean Lung", "HF Lung"), c("HF Lung", "HF_SPM Lung"))
p5 <- ggbarplot(sea_data, "X2..Group", y = "Maximum.Respiration", xlab = "Groups", ylab = "Maximum Respiration", label = FALSE, width = 0.5, add = c("mean_se"), palette = "jco", color = "X2..Group", title = "Maximum Respiration", legend = "None")+
  theme(axis.text=element_text(size=12))+
  #stat_compare_means(data = sea_data, method = "kruskal.test", label.y = 400, label.x = 1)
  stat_compare_means(comparisons = my_comparisons, method = "t.test", tip.length = c(0.03), label = "p.format", hide.ns = TRUE, label.y = c(155, 140)) #p.format to place actual p-values (uses wilcox.test)

####Stats & Graphing for Non-Mito Respiration####
shapiro.test(sea_data$Non.Mitochondrial.Respiraion)

kruskal.test(Non.Mitochondrial.Respiraion ~ X2..Group,  data = sea_data)
compare_means(Non.Mitochondrial.Respiraion ~ X2..Group,  data = sea_data, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 

#my_comparisons <- list(c("1", "2"), c("2", "3"), c("1", "3"))
p6 <- ggbarplot(sea_data, "X2..Group", y = "Non.Mitochondrial.Respiraion", xlab = "Groups", ylab = "Non-Mitochondrial Respiraion", label = FALSE, width = 0.5, add = c("mean_se"), palette = "jco", color = "X2..Group", title = "Non-Mitochondrial Respiraion", legend = "None")+
  theme(axis.text=element_text(size=12))+
  stat_compare_means(data = sea_data, method = "kruskal.test", label.y = 26, label.x = 1)
#  stat_compare_means(comparisons = my_comparisons, tip.length = c(0.001), label = "p.signif", hide.ns = TRUE, label.y = c(125, 127, 129, 131)) #p.format to place actual p-values (uses wilcox.test)

####Stats & Graphing for OCR/ECAR Ratio####
shapiro.test(sea_data$OCR.ECAR.Ratio)

kruskal.test(OCR.ECAR.Ratio ~ X2..Group,  data = sea_data)
compare_means(OCR.ECAR.Ratio ~ X2..Group,  data = sea_data, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 

#my_comparisons <- list(c("1", "2"), c("2", "3"), c("1", "3"))
p7 <- ggbarplot(sea_data, "X2..Group", y = "OCR.ECAR.Ratio", xlab = "Groups", ylab = "OCR/ECAR Ratio", label = FALSE, width = 0.5, add = c("mean_se"), palette = "jco", color = "X2..Group", title = "OCR/ECAR Ratio", legend = "None")+
  theme(axis.text=element_text(size=12))+
  stat_compare_means(data = sea_data, method = "kruskal.test", label.y = 15, label.x = 1)
#  stat_compare_means(comparisons = my_comparisons, tip.length = c(0.001), label = "p.signif", hide.ns = TRUE, label.y = c(125, 127, 129, 131)) #p.format to place actual p-values (uses wilcox.test)

####Stats & Graphing for Proton Leak####
shapiro.test(sea_data$Proton.Leak)

kruskal.test(Proton.Leak ~ X2..Group,  data = sea_data)
compare_means(Proton.Leak ~ X2..Group,  data = sea_data, method = "wilcox.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 

my_comparisons <- list(c("CD4 Lean", "CD4 SPM_HF"))
p8 <- ggbarplot(sea_data, "X2..Group", y = "Proton.Leak", xlab = "Groups", ylab = "Proton Leak", label = FALSE, width = 0.5, add = c("mean_se"), palette = "jco", color = "X2..Group", title = "Proton Leak", legend = "None")+
  theme(axis.text=element_text(size=12))+
  #stat_compare_means(data = sea_data, method = "kruskal.test", label.y = 15, label.x = 1)
  stat_compare_means(comparisons = my_comparisons, tip.length = c(0.01), label = "p.format", hide.ns = TRUE, label.y = c(17)) #p.format to place actual p-values (uses wilcox.test)

####Stats & Graphing for Spare Respiratory Capacity####
shapiro.test(sea_data$Spare.Respiratory.Capacity) 
hist(sea_data$Spare.Respiratory.Capacity)
qqnorm(sea_data$Spare.Respiratory.Capacity)
qqline(sea_data$Spare.Respiratory.Capacity)

kruskal.test(Spare.Respiratory.Capacity ~ X2..Group,  data = sea_data)
compare_means(Spare.Respiratory.Capacity ~ X2..Group,  data = sea_data, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 

my_comparisons <- list(c("Lean Lung", "HF Lung"), c("HF Lung", "HF_SPM Lung"))
p9 <- ggbarplot(sea_data, "X2..Group", y = "Spare.Respiratory.Capacity", xlab = "Groups", ylab = "Spare Respiratory Capacity", label = FALSE, width = 0.5, add = c("mean_se"), palette = "jco", color = "X2..Group", title = "Spare Respiratory Capacity", legend = "None")+
  theme(axis.text=element_text(size=12))+
  #stat_compare_means(data = sea_data, method = "kruskal.test", label.y = 350, label.x = 1)
  stat_compare_means(comparisons = my_comparisons, tip.length = c(0.01), method = "t.test", label = "p.format", hide.ns = TRUE, label.y = c(100,110)) #p.format to place actual p-values (uses wilcox.test)

####Stats & Graphing for Spare Respiratory Capacity %####
shapiro.test(sea_data$Spare.Respiratory.Capacity..) 

kruskal.test(Spare.Respiratory.Capacity.. ~ X2..Group,  data = sea_data)
compare_means(Spare.Respiratory.Capacity.. ~ X2..Group,  data = sea_data, method = "t.test") #if you wanted to make all pairwise comparisons to one reference group, then use parameter, ref.group. 

my_comparisons <- list(c("CD4 Lean", "CD4 SPM_HF"))
p10 <- ggbarplot(sea_data, "X2..Group", y = "Spare.Respiratory.Capacity..", xlab = "Groups", ylab = "Spare Respiratory Capacity %", label = FALSE, width = 0.5, add = c("mean_se"), palette = "jco", color = "X2..Group", title = "Spare Respiratory Capacity Percentage", subtitle = "Maximum Respiration / Basal OCR", legend = "None")+
  theme(axis.text=element_text(size=12))+font("subtitle", size = 12)+ #change subtitle font size
  #stat_compare_means(data = sea_data, method = "kruskal.test", label.y = 480, label.x = 1)
  stat_compare_means(comparisons = my_comparisons, tip.length = c(0.02),  method = "t.test", label = "p.format", hide.ns = TRUE, label.y = c(430)) #p.format to place actual p-values (uses wilcox.test)


####Combining All Plots####
ggarrange(p1,p2,p3,p4,p5,p6, ncol = 2, nrow = 3, legend = "none")
ggarrange(p7,p8,p9,p10, ncol = 2, nrow = 3, legend = "none")

#outputting all results to a PDF
pdf("Seahorse_Day21_COMBOplots_CD8_spleen_cells_12.11.18.pdf")
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
print(p8)
print(p9)
print(p10)
dev.off()

####Adding error bars manually from SEM file####
ggbarplot(means, "X2..Group", y = "Basal.ECAR", xlab = "Groups", ylab = "Basal ECAR", label = FALSE, width = 0.5, fill = "X2..Group", color = "X2..Group", title = "Basal ECAR", legend = "None")+
  stat_compare_means(data = sea_data, method = "kruskal.test", label.y = 8, label.x = 1)+      # Add global p-value
  geom_errorbar(aes(ymin=means$Basal.ECAR-sem$Basal.ECAR, ymax=means$Basal.ECAR+sem$Basal.ECAR), width=.12,
                position=position_dodge(.9))
