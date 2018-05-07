%Orbit ECEF ECI Coordinate conversion
% Richard Rieber
% October 1, 2006
% rrieber@gmail.com
%
% Revision 10/1/09: Preallocated memory for ECEF and V_ECEF vectors before
%                   the for-loop
%
% function [ECEF,V_ECEF] = eci2ecef(ECI, GST, V_ECI)
%
% Purpose:  This function rotates Earth Centered Inertial (ECI) coordinates to Earth
%           Centered Earth Fixed (ECEF) coordinates via the Greenwich Sideral Time
%           hour angle (GST).
% 
% Inputs:  ECI     - A 3 x n matrix of position vectors in km
%          GST     - A vector of length n providing the Greenwich hour angle for each of
%                    the above ECI position vectors in radians
%          V_ECI   - A 3 x n matrix of the velocity vectors in km/s [OPTIONAL]
% 
% Outputs:  ECEF   - A 3 x n matrix of position vectors in the ECEF coordinate system in km
%           V_ECEF - A 3 x n matrix of velocity vectors in the ECEF coordinate system in km/s
%
% NOTE:  This function requires the use of the subfunction 'R3.m' which creates a
%        rotation matrix about the 3-axis (Z-axis).


function [ECEF,V_ECEF] = eci2ecef(ECI, GST, V_ECI)

%Checking number of inputs for errors
if nargin < 2 || nargin > 3
    error('Incorrect number of inputs.  See help eci2ecef.')
end

[b,n] = size(ECI);
%Checking to see if length of ECI matrix is the same as the length of the GST vector
if n ~= length(GST)
    error('Size of ECI vector not equal to size of GST vector.  Check inputs.')
end

%Checking to see if the ECI vector has 3 elements
if b ~= 3
    error('ECI vector must have 3 elements to the vector (X,Y,Z) coordinates')
end

ECEF = zeros(3,n);

if nargin == 3
    V_ECEF = zeros(3,n);
end

for j = 1:n  %Iterating thru the number of positions provided by user
    % Rotating the ECI vector into the ECEF frame via the GST angle about the Z-axis
    ECEF(:,j) = R3(GST(j))*ECI(:,j);
    if nargin == 3
        V_ECEF(:,j) = R3(GST(j))*V_ECI(:,j);
    end
end