
Processing for Clinical Outcome Scores from the Inpatient Sensor Data

Adam P. Horin
June 30, 2020

-----------------------------------------------
General Description of Matlab Processing
-----------------------------------------------
All main scripts in this folder process the clinical scores for their corresponding test.
Input files:
  CVA "\\fs2.smpp.local\RTO\Inpatient Sensors -Stroke\MC10 Study\Outcome Measures\individual_tests_CVA"
  HC "\\fs2.smpp.local\RTO\Inpatient Sensors -Stroke\MC10 Study\Outcome Measures\individual_tests_HC"

Every clinical outcome has a maximum of 4 sessions (SN). SN1 refers to the admission session and SN4 refers to the discharge session.
If a participant had only 2 or 3 sessions, their last session was changed to SN4. Every script identifies the date for each session and checks
to make the sure the order of the sessions is correct. A matrix that is generated during processing that identifies any dates that were out of 
order (List_WrongDate) so you can manually confirm the corrected order. If there were multiple trials for a single session, only the first trial
was processed.


-----------------------------------------------
BBS
-----------------------------------------------
This script processes the BBS subscores and total scores for each BBS trial for every session.
Input files:
  CVA "\\fs2.smpp.local\RTO\Inpatient Sensors -Stroke\MC10 Study\Outcome Measures\individual_tests_CVA\scores_cva_BBS"
  HC "\\fs2.smpp.local\RTO\Inpatient Sensors -Stroke\MC10 Study\Outcome Measures\individual_tests_HC\scores_hc_BBS"

Output matrix:
subject: subject ID
group: CVA or CONTROLS
activity: BBS
session: 1-4
trial_no: 1-14
BBS_total: total score of all 14 tasks
BBS_Subscore: score for each task


-----------------------------------------------
6MWT
-----------------------------------------------




-----------------------------------------------
10MWT
-----------------------------------------------





-----------------------------------------------
TUG
-----------------------------------------------






