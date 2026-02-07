function [x, y, z] = sph2cartcoord(phi, theta, r_mag)
    x = r_mag.*sin(theta).*cos(phi);
    y = r_mag.*sin(theta).*sin(phi);
    z = r_mag.*cos(theta);
end