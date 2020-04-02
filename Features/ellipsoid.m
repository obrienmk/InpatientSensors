function [center,angles,R,area,axis]=ellipsoid(data)

% Citation:
% [1]. Schubert, P., & Kirchner, M. (2014). Ellipse area calculations and their
% applicability in posturography. Gait and Posture, 39(1), 518–522.
% https://doi.org/10.1016/j.gaitpost.2013.09.001
% -------------------------------------------------------------------------
% Inputs
% data= [Mediolateral Acceleration (x), Anteroposterior Acceleration (y)]

% Outputs
% center = Center of the Ellipsoid
% angles = Ellipsoid Angles
% R = Rotation matrix || eigen vectors of the covariance matrix of "data"
% area = Ellipsoid Area
% axis = Ellipsod axis || Eigen Values of the covariance matrix of "data"
%--------------------------------------------------------------------------

[n, p] = size(data);                          % 2-D array dimensions
covar = cov(data);                            % Covariance matrix of data
[~, EigVal, EigVec] = svd(covar);             % singular value decomposition S =eigen values V=eigen vectors
f95 = finv(.95,p,n-p)*(n-1)*p*(n+1)/n/(n-p);  % F 95 percent point function
axis = sqrt(diag(EigVal)*f95);                % Axis
area = pi^(p/2)/gamma(p/2+1)*prod(axis);      % Area

R=EigVec;                                     % Rotation Matrix
center=mean(data,1);
angles=[rad2deg(atan2(R(2,1),R(1,1))),90-(rad2deg(atan2(R(2,1),-R(1,1))))];

end