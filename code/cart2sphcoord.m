function [r, phi, theta] = cart2sphcoord(x,y,z)
    r = norm([x,y,z]);
    phi = atan2(y, x);
    theta = acos(z/r);
end