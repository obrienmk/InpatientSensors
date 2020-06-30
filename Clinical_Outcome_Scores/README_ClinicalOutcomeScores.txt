
Processing Stream for Clinical Outcome Scores from the Inpatient Sensor Data

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
This script processes the 6MWT scores for every minute for every session.
Input files:
  CVA "\\fs2.smpp.local\RTO\Inpatient Sensors -Stroke\MC10 Study\Outcome Measures\individual_tests_CVA\scores_cva_6MWT"
  HC "\\fs2.smpp.local\RTO\Inpatient Sensors -Stroke\MC10 Study\Outcome Measures\individual_tests_HC\scores_hc_6MWT"

Output matrix:
  subject: subject ID
  group: CVA or CONTROLS
  activity: 6MWT
  session: 1-4
  minute: 1-6
  step count: number of steps at each minute (cumulative) 
  distance_m: distance in meters at each minute (cumulative)

-----------------------------------------------
10MWT
-----------------------------------------------
This script processes the 10MWT scores for each test (SSV,FV) for every session.
Input files:
  CVA "\\fs2.smpp.local\RTO\Inpatient Sensors -Stroke\MC10 Study\Outcome Measures\individual_tests_CVA\scores_cva_10MWT"
  HC "\\fs2.smpp.local\RTO\Inpatient Sensors -Stroke\MC10 Study\Outcome Measures\individual_tests_HC\scores_hc_10MWT"

Output matrix: 
Note: This output matrix is unstacked by the test column (SSV, FV) before exported.
  subject: subject ID
  group: CVA or CONTROLS
  activity: 10MWT
  session: 1-4
  test: SSV or FV
  step count: number of steps
  distance_m: distance in meters
  time_s: total time
  velocity: gait speed in m/s calculated from distance_m/time_s

-----------------------------------------------
TUG
-----------------------------------------------
This script processes the TUG scores for every session.
Input files:
  CVA "\\fs2.smpp.local\RTO\Inpatient Sensors -Stroke\MC10 Study\Outcome Measures\individual_tests_CVA\scores_cva_TUG"
  HC "\\fs2.smpp.local\RTO\Inpatient Sensors -Stroke\MC10 Study\Outcome Measures\individual_tests_HC\scores_hc_TUG"

Output matrix:
  subject: subject ID
  group: CVA or CONTROLS
  activity: TUG
  session: 1-4
  step count: number of steps
  time_s: total time

