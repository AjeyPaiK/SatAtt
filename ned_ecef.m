function ecef_vec = ned_ecef(ned_vec,lat,lon)
    rot=[-cosd(lon)*sind(lat) -sind(lon) -cosd(lon)*cosd(lat);
        -sind(lon)*sind(lat) cosd(lon) -sind(lon)*cosd(lat);
        cosd(lat) 0 -sind(lat)];

    ecef_vec = rot * ned_vec';
end