function [ tasks ] = get_tasks( Activity )
%**************************************************************************
% get_tasks.m
%
% Written by Megan O'Brien
% version 6/1/2018
%
% Get list of relevant tasks (post-label survey) for a given labelled
% activity in the BiostampRC app
%**************************************************************************

% ***** MAS *****
if strcmp(Activity,'Clinical - MAS')
    MAS_activities = {'HIP FLEXION','HIP EXTENSION','KNEE FLEXION','KNEE EXTENSION','DORSIFLEXION','PLANTARFLEXION'};
    
    tasks = MAS_activities;
end

% ***** MMT *****
if strcmp(Activity,'Clinical - MMT')
    MMT_activities = {'HIP FLEXION','HIP EXTENSION','KNEE FLEXION','KNEE EXTENSION','DORSIFLEXION','PLANTARFLEXION',...
        'ELBOW FLEXION','SHOULDER FLEXION','SUPINATION'};
    
    tasks = MMT_activities;
end

% ***** BBS *****
if strcmp(Activity,'Clinical - BBS')
    BBS_activities = {'1. SIT TO STAND','2. STAND UNSUPPORTED','3. SIT W/ BACK','4. STAND TO SIT','5. TRANSFERS',...
        '6. STAND W/ EYES CLOSED','7. STAND W/ FEET TOGETHER','8. REACH FORWARD','9. PICK UP OBJECT','10. TURN',...
        '11. TURN 360 DEGS','12. FOOT ON STEP/STOOL','13. STAND W/ ONE FOOT IN FRONT','14. STAND ON ONE LEG'};
    
    tasks = BBS_activities;
end
    
% ***** 10MWT *****
if strcmp(Activity,'Clinical - 10MWT SSV') || strcmp(Activity,'Clinical - 10MWT FV')
    MWT10_activities = {'Trial'};
    
    tasks = MWT10_activities;
end

% ***** TUG *****
if strcmp(Activity,'Clinical - TUG')
    TUG_activities = {'Trial'};
    
    tasks = TUG_activities;
end

% ***** 6MWT *****
if strcmp(Activity,'Clinical - 6MWT')
    MWT6_activities = {'Trial'};
    
    tasks = MWT6_activities;
end


% ***** Resting ECG *****
if strcmp(Activity,'Resting ECG')
    ECGrest_activities = {'PRE EXERCISE','POST EXERCISE'};
    
    tasks = ECGrest_activities;
end


% ***** Activity Recognition *****
if strcmp(Activity,'Activity Recognition')
    AR_activities = {'LYING','SITTING','STANDING','WALKING','STAIRS UP','STAIRS DOWN','SIT TO STAND','STAND TO SIT'};
    
    tasks = AR_activities;
end


% ***** Physical Therapy *****
if strcmp(Activity,'Physical Therapy')
    PT_activities = {'WALKING ON A TREADMILL','WALKING UP AND DOWN STAIRS','CYCLING'};
    
    tasks = PT_activities;
end

end

