function [grad_B,doB_dot] = bGrad(coord)
time = [2017 2 22 5 38 0];
datetime = coord(4:9);
lla = coord(1:3);
r_eci = lla_eci(lla,time);

x = r_eci(1);
y = r_eci(2);
z = r_eci(3);

delx = 100;
dely = 100;
delz = 100; %in metres
delt = 30; %minutes

dox_plus = r_eci' + [delx,0,0];
doy_plus = r_eci' + [0,dely,0];
doz_plus = r_eci' + [0,0,delz];
%dot_plus = time + delt;

dot_plus = [2017 2 22 5 5 0];

dox_minus = r_eci' - [delx,0,0];
doy_minus = r_eci' - [0,dely,0];
doz_minus = r_eci' - [0,0,delz];
dot_minus = datetime - [0,0,0,0,delt,0];

a1 = eci_lla(dox_plus,datetime);
b1 = eci_lla(doy_plus,datetime);
c1 = eci_lla(doz_plus,datetime);
d1 = eci_lla(r_eci',dot_plus);

a2 = eci_lla(dox_minus,datetime);
b2 = eci_lla(doy_minus,datetime);
c2 = eci_lla(doz_minus,datetime);
d2 = eci_lla(r_eci',dot_minus);

%p = wrldmagm(680000,12.9716,77.5946,decyear(2016,09,02),'2015')
q = igrf_eci(a1(1),a1(2),a1(3),datetime);
r = igrf_eci(a1(1),a1(2),a1(3),datetime);
doBx =  (q - r)/(2*delx);

q = igrf_eci(b1(1),b1(2),b1(3),datetime);
r = igrf_eci(b2(1),b2(2),b2(3),datetime);
doBy =  (q - r)/(2*dely);

q = igrf_eci(c1(1),c1(2),c1(3),datetime);
r = igrf_eci(c2(1),c2(2),c2(3),datetime);
doBz =  (q - r)/(2*delz);

q = igrf_eci(d1(1),d1(2),d1(3),datetime);
r = igrf_eci(d2(1),d2(2),d2(3),datetime);
doBt =  (q - r)/(2*delt);


grad_B = [doBx , doBy , doBz];
doB_dot = doBt;

end