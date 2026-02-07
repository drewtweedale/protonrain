function B = sph2cart_field(Br, Btheta, Bphi, r)
    x = r(1); y = r(2); z = r(3);
    r_mag = norm(r);
    theta = acos(z/r_mag);
    phi = atan2(y,x);
    
    Bx = Br*sin(theta)*cos(phi) + Btheta*cos(theta)*cos(phi) - Bphi*sin(phi);
    By = Br*sin(theta)*sin(phi) + Btheta*cos(theta)*sin(phi) + Bphi*cos(phi);
    Bz = Br*cos(theta) - Btheta*sin(theta);
    B = [Bx, By, Bz];
end