function att = att_det(B_body,gps,ang_vel)
mu=398600.436233;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Prefilter
x_pf_old = [0.534;0.5678;0.3456;0.543;0.52;0.35;0;0;0];
P_pf_old = diag([10^4 10^4 10^4 10^6 10^6 10^6 10^6 10^6 10^6]);
[x_pf ,~] = prefilter(B_body, x_pf_old, P_pf_old);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Position and Velocity of the Satellite in ECI frame
oev=[700000,0,1.715658655,0,0,1.047197551]; %Generated from orbit design
[~,v_eci] = orb2eci(mu,oev);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%B_eci using IGclcRF
% lat, long,h, datestring obtained from GPS

lat = gps(1); % degrees
long = gps(2); % degrees
h_km = gps(3); % km
h_m = h_km*1000;
datestring = [2017,2,22,5,38,0]; %'dd-mm-yyyy HH:MM:SS'
B = igrf_eci(lat,long,h_km,datestring); %in geodetic coordinate system(nT)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%B_eci_dot
coord = [lat, long, h_m, datestring];
[grad_B,doB_dot] = bGrad(coord);
B_eci_dot = grad_B*v_eci + doB_dot;
Bi  = [B; B_eci_dot]; % 6x1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%attitude filter
q_old=[0.28222;0.56443;0.18814;0.75258];
w_old=[2.2;5.2;3.3];
x_att_old = [q_old;w_old];
P_att_old = diag([1*10^-3 1*10^-3 1*10^-3 1*10^-3 1*10^-4 1*10^-4 1*10^-4]);
[x_new,~] = attitude_filter(x_pf(1:7), x_att_old, P_att_old,Bi);
q=x_new(1:4,1);
q=q';
[Yaw,Pitch,Roll]=quat2angle(q);
att=[Yaw,Pitch,Roll]
end