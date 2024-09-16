function distance = calculate_distance(lat1, lon1, lat2, lon2)
    R = 6371000;  % 地球半徑（米）
    phi1 = deg2rad(lat1);
    phi2 = deg2rad(lat2);
    delta_phi = deg2rad(lat2 - lat1);
    delta_lambda = deg2rad(lon2 - lon1);

    a = sin(delta_phi / 2) ^ 2 + cos(phi1) * cos(phi2) * sin(delta_lambda / 2) ^ 2;
    c = 2 * atan2(sqrt(a), sqrt(1 - a));

    distance = R * c;
end