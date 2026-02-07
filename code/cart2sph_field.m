function B = cart2sph_field(Bx, By, Bz, theta, phi)
    Br = Bx*sin(theta)*cos(phi) + By*sin(theta)*sin(phi) + Bz*cos(theta);
    Btheta = Bx*cos(theta)*cos(phi) + By*cos(theta)*sin(phi) - Bz*sin(theta);
    Bphi = By*cos(phi) - Bx*sin(phi);
    B = [Br, Btheta, Bphi];
end