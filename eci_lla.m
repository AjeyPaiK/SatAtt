function lla = eci_lla(eci,utc_time)
    jd = jday(utc_time(1),utc_time(2),utc_time(3),utc_time(4),utc_time(5),utc_time(6));
    gst = JD2GAST(jd)/360*2*pi;
    ecef  = eci2ecef(eci',gst);
    lla = ecef2geod(ecef);
end
