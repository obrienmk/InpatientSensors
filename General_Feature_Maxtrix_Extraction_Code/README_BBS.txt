Processing Stream for Extracting General and Balance Features from BBS Inpatient Sensor Data
Adam P. Horin
June 25, 2020

------------------------------------
Matlab Processing
------------------------------------
# Matlab script description: 
The main_General_Feature_Matrix_Admission_BBS.m extracts general and balance features from the BBS inpatient sensor data.
The output is a feature matrix containing groups (CVA and HC), BBS activities (1-14), and cutoff times (15-120, and 200 seconds).

The main_General_Feature_Matrix_AllSessions_BBS.m extracts general and balance features from the BBS inpatient sensor data.
The output is a feature matrix containing groups (CVA and HC), session (1-4), BBS activities (1-14), and cutoff times (15-120, and 200 seconds).
Note: The sessions are defined SN1, SN2, SN3, and SN4; SN1 is admission, SN4 is discharge. If a participant has only 2 or 3 sessions, their last session was defined as SN4.

# Input files and locations:
BBS CVA and HC clean extracted data: "Y:\Inpatient Sensors -Stroke\MC10 Study\Data analysis\2_Clean_Data_Extracted\BBS"
BBS CVA Scores: "\\fs2.smpp.local\RTO\Inpatient Sensors -Stroke\MC10 Study\Outcome Measures\individual_tests_CVA\scores_cva_BBS.xlsx"
BBS CVA Scores: "\\fs2.smpp.local\RTO\Inpatient Sensors -Stroke\MC10 Study\Outcome Measures\individual_tests_HC\scores_hc_BBS.xlsx"

# Main script:
main_General_Feature_Matrix_Admission_BBS.m

# Function scripts:
FFeatures.m
FFT_data.m
acctransformation.m
ellipsoid.m
sampen.m
staticBalance.m

--------------------------
Feature Matrix
--------------------------
# Feature matrix output description:
The feature matrix output from the main matlab script above contains general features and balance features for all of the BBS
activities at 15 second increments ranging from 15-120 seconds.

# Output file and location:

------------------------------------------------
Adding Demographic Information
------------------------------------------------
# Steps to add demographic information to the feature matrix:
1. To add demographic information to the feature matrix, add appropriate columns and
add column headings for the information to be included.
2. Use the following formula (example below) in the appropriate cells of the feature matrix spreadsheet
to match the index of the subject ID of the demographic spreadsheet.

=INDEX('[CVA_Clinical Outcome Measure Scores Master Shee - Yaar HARARI.xlsm]Final'!$E$2:$E$56,MATCH($A2,'[CVA_Clinical Outcome Measure Scores Master Shee - Yaar HARARI.xlsm]Final'!$A$2:$A$56,0))

# Input files and locations:
CVA Demographic Information: "\\fs2.smpp.local\RTO\Inpatient Sensors -Stroke\MC10 Study\Outcome Measures\CVA_Clinical Outcome Measure Scores Master Shee - Yaar HARARI.xlsm"
HC Demographic Information: "\\fs2.smpp.local\RTO\Inpatient Sensors -Stroke\MC10 Study\Outcome Measures\HC_Clinical Outcome Measure Scores Master Sheet-updated 1-11-18.xlsx"

# Final data file and location:

------------------------------------------------
Adding Clinical Outcome Scores
------------------------------------------------
# Steps to add clinical outcome scores:
1. These may be matched based on multiple criteria in order to add values from different session for each participant. Add appropriate columns to the feature matrix.
2. Use the following formula (example below) to match based on multiple criteria. To execute this formula you must press ctrl + shift + enter after you have finished entering the formula.

=INDEX(BBS_ClinicalScores.csv!$G$2:$G$127,MATCH(1,(A2=BBS_ClinicalScores.csv!$A$2:$A$3480)*(D2=BBS_ClinicalScores.csv!$D$2:$D$3480)*(E2=BBS_ClinicalScores.csv!$E$2:$E$3480),0))

