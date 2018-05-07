load gpsdata.mat;
mu=398600.436233;
TESTPREFILTER

%---------------- Initializing Orbital Elemental Vectors --------------%

oev=[20095,0.00401,0.3517*pi,0,0,0];
%  oev(1) = semimajor axis (kilometers)
%  oev(2) = orbital eccentricity (non-dimensional)
%           (0 <= eccentricity < 1)
%  oev(3) = orbital inclination (radians)
%           (0 <= inclination <= pi)
%  oev(4) = argument of perigee (radians)
%           (0 <= argument of perigee <= 2 pi)
%  oev(5) = right ascension of ascending node (radians)
%           (0 <= raan <= 2 pi)
%  oev(6) = true anomaly (radians)
%           (0 <= true anomaly <= 2 pi)

%----------------- End Initialization -----------------------------------%

%------------------- Simulation of magnetic field -----------------------%

B=wgn(1,25,0.2);               % Assuming these as Magnetometer values
E=B+wgn(1,25,0.2);             % Measured signal with  measurement noise
E=E';
[~,v_eci] = orb2eci(mu,oev);   % Generating Velocity Vector and Position Vector
%------------------------------------------------------------------------%

%------------------ Initializing state vectors --------------------------%
x_pf_old = [0.534;0.5678;0.3456;0.543;0.52;0.35;0;0;0];
P_pf_old = diag([10^4 10^4 10^4 10^6 10^6 10^6 10^6 10^6 10^6]);
%------------------------------------------------------------------------%

for i=1:20
    lat=gpsdata(i,1);
    long=gpsdata(i,2);
    h_km=gpsdata(i,3);
    h_m=1000*h_km;
    if i<=20
        datestring=gpsdata(i,4:end);
    end
    %--------------Generating the IGRF values in ECI frame--------------%
    
    igrf(:,i)=igrf_eci(lat,long,h_km,datestring);
    coord(i,:) = [lat, long, h_m, datestring];
    [grad_B,doB_dot] = bGrad(coord);
    B_eci_dot = grad_B*v_eci + doB_dot;
        Bi  = [E(i:i+2); B_eci_dot]; % 6x1
        [x_pf(:,i),~] = prefilter(E(i:i+2), x_pf_old, P_pf_old); %Invoking Prefilter
        
    %-----------Initializing state vectors for Attitude Filter---------%
    q_old=[0.28222;0.56443;0.18814;0.75258];
    w_old=[2.2;5.2;3.3];
    x_att_old = [q_old;w_old];  % State Vector for Attitude Filter
    P_att_old = diag([1*10^-3 1*10^-3 1*10^-3 1*10^-3 1*10^-4 1*10^-4 1*10^-4]);
    [x_new(:,i),~] = attitude_filter(x_pf(:,i), x_att_old, P_att_old,Bi); %Invoking Attitude filter
    q=x_new(1:4,i);
    q=q';
    [Yaw,Pitch,Roll]=quat2angle(q); % Conversion of Quaternions to Euler Angles
    att(i,:)=[Yaw,Pitch,Roll];
end

%----------------- Plotting Parameters ------------------------------------%
AngMom=x_new(5:7,:);
Yaw=att(:,1);
Pitch=att(:,2);
Roll=att(:,3);
figure(1);
subplot(2,1,1);
title('3-Axis parameters');
plot(Yaw,'g');
hold on;
plot(Pitch,'b');
hold on;
plot(Roll,'r');
hold off;
grid on;
legend('Yaw','Pitch','Roll','Location','east');
xlabel('Time(s)');
ylabel('Angle (degrees)');
subplot(2,1,2);
title('Angular Momentum of the Satellite');
plot(AngMom(1,:),'g');
hold on;
plot(AngMom(2,:),'b');
hold on;
plot(AngMom(3,:),'r');
hold off;
grid on;
xlabel('Time(s)');
ylabel('Angular momentum(kgm/s).');
legend('Momentum along X','Momentum along Y','Momentum along Z');