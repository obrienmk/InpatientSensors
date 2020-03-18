function  transformAcc= acctransformation(test,sbname,session_no)
% -------------------------------------------------------------------------
% METHOD'S SUMMARY:
% An accelerometer positioned over the lumbar spine of a moving subject,
% may be at a tilt because of inaccurate positioning of the device in the
% frontal plane, and because of lumbar curvature and forward leaning
% in the sagittal plane. The algorithms removes and the static gravity
% component estimating the acceleration vectors in a HORIZONTAL-VERTICAL
% ortoganl coordinate system. ASSUMING THAT ONCE THE SENSOR IS PLACED it
% will not move so the angles will BE CONSTANT.

% CITATIONS:
% [1]. Moe-Nielssen, R., Moe-Nilssen, R., & Moe-Nielssen, R. (1998).
% A new method for evaluating motor control in gait under real-life
% environmental conditions. Part 2: Gait analysis. Clinical Biomechanics
% (Bristol, Avon), 13(4–5), 328–335. https://doi.org/10.1016/S0268-0033(98)00089-8

% Notes:
% In the paper the AP direction is POSITIVE going forward
% so in this study the AP is negative going forward so it will be
% multiplied by (-)

% Inputs:
% aap = measured anteroposterior acceleration
% aml = measured medio lateral acceleration
% av  = measured vertical acceleration
% sbname = name of the subject
% session_no = session #

% Outputs:
% Estimations or transform accelerations
% aAP = Estimated anteroposterior acceleration
% aML = Estimated medio lateral acceleration
% aV  = Estimated vertical acceleration
% aV_prime  = Estimated TEMPORAL vertical acceleration
% theta_ap = Angle between aap and aAP
% theta_ml = Angle between aml and aML
% -------------------------------------------------------------------------

syms   aAP aML aV aV_prime theta_ap theta_ml real


if strcmp(sbname, 'HC06') || (strcmp(sbname, 'CVA20') && session_no == 2)  % CVA20 on second day of testing
    
    % ***** NEED TO FIX EQNS FOR HC06 & CVA20 b/c sensor was vertical
    aap=-test.Acc.z';
    aml=test.Acc.y';
    av=-test.Acc.x';
    
    angY=test.Gyr.z';
    angX=test.Gyr.y';
    angZ=test.Gyr.x';
    
else
    
    aap=test.Acc.z';
    aml=test.Acc.x';
    av=test.Acc.y';
    
    angY=test.Gyr.y';
    angX=test.Gyr.x';
    angZ=test.Gyr.z';
    
end

%Dynamic Accelerations (No Gravity Contribution)
aAP= -aap*cos(theta_ap)+av*sin(theta_ap);
aV_prime = aap*sin(theta_ap)+ av*cos(theta_ap);
aML = -aml*cos(theta_ml) + aV_prime*sin(theta_ml);
aV = aml*sin(theta_ml)+aV_prime*cos(theta_ml)-1; % -1 means substracting 1g

% According to the study the mean of a given acceleration(a_j) is equal to the
% sin(theta_j) so by arcsin(theta_j) the angle can be estimated!

% Mean Acc values
aap_mean= mean(aap);
aml_mean= mean(aml);

theta_ap= asin(aap_mean);
theta_ml= asin(aml_mean);

aAP= eval(aAP);
aV_prime = eval(aV_prime);
aML = eval(aML);
aV = eval(aV);

%% STORE values
transformAcc.AP=aAP*9.80665;
transformAcc.ML=aML*9.80665;
transformAcc.V=aV*9.80665;
transformAcc.angX=angX;
transformAcc.angY=angY;
transformAcc.angZ=angZ;
end