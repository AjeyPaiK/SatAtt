%Orbit ECEF ECI Coordinate conversion
%
%
% function [ECI, V_ECI] = ecef2eci(ECEF, GST, V_ECEF)
% 
% Purpose:  This function rotates Earth Centered Earth Fixed (ECEF) coordinates
%           to Earth Centered Inertial (ECI) coordinates via the Greenwich Sideral Time
%           hour angle (GST).
% 
% Inputs:  ECEF   - A 3 x n matrix of position vectors in the ECEF frame in km
%          GST    - A vector of length n providing the Greenwich hour angle for each of
%                   the above ECI position vectors in radians
%          V_ECEF - A 3 x n matrix of the velocity vectors in km/s [OPTIONAL]
% 
% Outputs:  ECI   - A 3 x n matrix of position vectors in the ECI coordinate system in km
%           V_ECI - A 3 x n matrix of velocity vectors in the ECI coordinate system
%                   in km/s [OPTIONAL]
%
% NOTE:  This function requires the use of the subfunction 'R3.m' which creates a
%        rotation matrix about the 3-axis (Z-axis).

function [ECI, V_ECI] = ecef2eci(ECEF, GST, V_ECEF)

%Checking number of inputs for errors
if nargin < 2 || nargin > 3
    error('Incorrect number of inputs.  See help eci2ecef.')
end

[b,n] = size(ECEF);
%Checking to see if length of ECEF matrix is the same as the length of the GST vector
if n ~= length(GST)
    error('Size of ECEF vector not equal to size of GST vector.  Check inputs.')
end

%Checking to see if the ECEF vector has 3 elements
if b ~= 3
    error('ECEF vector must have 3 elements to the vector (X,Y,Z) coordinates')
end

ECI = zeros(3,n);

if nargin == 3
    V_ECI = zeros(3,n);
end

for j = 1:n  %Iterating thru the number of positions provided by user
    % Rotating the ECEF vector into the ECEF frame via the GST angle about the Z-axis
    ECI(:,j) = R3(-GST(j))*ECEF(:,j);
    if nargin == 3
        V_ECI(:,j) = R3(-GST(j))*V_ECEF(:,j);
    end
end