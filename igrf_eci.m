function B = igrf_eci(lat,long,h,datetime)

%%%%%%%%%%%%%% Calculates Magnetic field in Spherical Coordinates %%%%%%%%%%%

%
% Inputs
% r Geocentric radius
% theta Latitude measured in degrees positive from equator
% phi Longitude measured in degrees positive east from Greenwich
% datestring in 'dd-mm-yy HH:MM:SS' format Eg. '01-03-2017 15:45:17'
%
% Outputs - magnetic field strength in local tangential coordinates
% Br B in radial direction
% Bt B in theta direction
% Bp B in phi direction
% Checks to see if located at either pole to avoid singularities

if (lat>-0.00000001 && lat<0.00000001)
    lat=0.00000001;
elseif(lat<180.00000001 && lat>179.99999999)
    lat=179.99999999;
end

costheta = cos((90 - lat)*pi/180);
sintheta = sin((90 - lat)*pi/180);

 ra = 6378.137; f = 1/298.257223563; b = ra*(1 - f);
    rho = hypot(ra*sintheta, b*costheta);
    r = sqrt( h.^2 + 2*h.*rho + ...
        (ra^4*sintheta.^2 + b^4*costheta.^2) ./ rho.^2 );
    cd = (h + rho) ./ r;
    sd = (ra^2 - b^2) ./ rho .* costheta.*sintheta./r;
    oldcos = costheta;
    costheta = costheta.*cd - sintheta.*sd;
    sintheta = sintheta.*cd + oldcos.*sd;
    
 phi = long*pi/180;


% [X,Y,Z] = geodetic2ecef(referenceEllipsoid('earth'),lat,long,h);
% [phi,theta,r] = cart2sph(X,Y,Z);
% r=r/1000;


givendate=datetime;
setdate=[2015 1 1 0 0 0];
seconds=etime(givendate,setdate);
days=seconds/86400;


% The angles must be converted from degrees into radians
% theta=(pi/2)-theta;

a=6371.2; % Reference radius used in IGRF
% This section of the code simply reads in the g and h Schmidt
% quasi-normalized coefficients

[gn, gm, gvali, gsvi] = textread('igrfSg_sq_normalised2015.txt','%f %f %f %f');
[hn, hm, hvali, hsvi] = textread('igrfSh_sq_normalised2015.txt','%f %f %f %f');

% igrfSg = importfile('igrfSg_sq_normalised2015.txt', 1, 104);
% gn=igrfSg(:,1);
% gm=igrfSg(:,2);
% gvali=igrfSg(:,3);
% gsvi=igrfSg(:,4);
% 
% igrfSh = importfile('igrfSh_sq_normalised2015.txt', 1, 104);
% hn=igrfSh(:,1);
% hm=igrfSh(:,2);
% hvali=igrfSh(:,3);
% hsvi=igrfSh(:,4);

N=max(gn);
g=zeros(N,N+1);
hd=zeros(N,N+1);
for x=1:length(gn)
    g(gn(x),gm(x)+1) = gvali(x) + gsvi(x)*days/365;
    hd(hn(x),hm(x)+1) = hvali(x) + hsvi(x)*days/365;
end
% Initialize each of the variables
% Br B in the radial driection
% Bt B in the theta direction
% Bp B in the phi direction
% P The associated Legendre polynomial evaluated at cos(theta)
% The nomenclature for the recursive values generally follows
% the form P10 = P(n-1,m-0)
% dP The partial derivative of P with respect to theta
Br=0; Bt=0; Bp=0;
P11=1; P10=P11;
dP11=0; dP10=dP11;
for m=0:N
    for n=1:N
        if m<=n
            % Calculate Legendre polynomials and derivatives recursively
            if n==m
                P2 = sintheta*P11;
                dP2 = sintheta*dP11 + costheta*P11;
                P11=P2; P10=P11; P20=0;
                dP11=dP2; dP10=dP11; dP20=0;
            elseif n==1
                  P2 = costheta*P10;
                  dP2 = costheta*dP10 - sintheta*P10;
                  P20=P10; P10=P2;
                  dP20=dP10; dP10=dP2;
            else
                K = ((n-1)^2-m^2)/((2*n-1)*(2*n-3));
                P2 = costheta*P10 - K*P20;
                dP2 = costheta*dP10 - sintheta*P10 - K*dP20;
                P20=P10; P10=P2;
                dP20=dP10; dP10=dP2;
            end
            
            % Calculate Br, Bt, and Bp
            Br = Br + (a/r)^(n+2)*(n+1)*...
                ((g(n,m+1)*cos(m*phi) + hd(n,m+1)*sin(m*phi))*P2);
            Bt = Bt + (a/r)^(n+2)*...
                ((g(n,m+1)*cos(m*phi) + hd(n,m+1)*sin(m*phi))*dP2);
            Bp = Bp + (a/r)^(n+2)*...
                (m*(-g(n,m+1)*sin(m*phi) + hd(n,m+1)*cos(m*phi))* P2);
        end
    end
end
Bt =-Bt;
Bp= -(Bp/(sintheta));

% %NED in geocentric coordinates
 Bx = -Bt;
 By = Bp;
 Bz = -Br;
% %NED in geodectic coordinates
 Bx_old = Bx;
 Bx = Bx.*cd + Bz.*sd;
 By=By;
 Bz = Bz.*cd - Bx.*sd;
 
 B_ned = [Bx By Bz];
 
 %%%%% Co-ordinates Conversion %%%%%%%%%%%
 B_ecef = ned_ecef(B_ned, lat,long);
 jd = jday(givendate(1) ,givendate(2), givendate(3), givendate(4), givendate(5), givendate(6));
 gst = JD2GAST(jd);
 gst_radians = gst*pi/180;
 B = ecef2eci(B_ecef,gst_radians);
end