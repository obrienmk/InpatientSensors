clc
clear all
close all


AoM = import_AoM('./Feature_Matrix_AoM_Admission_6MWT.csv')
GF = import_GF('./General_Feature_Matrix_Admission_6MWT.csv')

GF.ID = []

FM = [AoM GF]

writetable(FM,'Feature_Matrix_Admission_6MWT.csv','Delimiter',',','QuoteStrings',true)
