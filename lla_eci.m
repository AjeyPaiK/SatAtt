function eci = lla_eci(lla,utc_time)
    jd = jday(utc_time(1),utc_time(2),utc_time(3),utc_time(4),utc_time(5),utc_time(6));
    gst = JD2GAST(jd)/360*2*pi;
    ecef = geod2ecef(lla(1),lla(2), lla(3));
    eci  = ecef2eci(ecef',gst);
end
