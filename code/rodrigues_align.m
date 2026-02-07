function R = rodrigues_align(u, v)
    % R such that R*u = v
    u = u(:) / norm(u);
    v = v(:) / norm(v);
    c = dot(u, v);

    axis = cross(u, v); 
    axis = axis / norm(axis);
    ang  = acos(c);
    K = [   0       -axis(3)  axis(2);
          axis(3)      0     -axis(1);
         -axis(2)   axis(1)     0    ];
    R = eye(3) + sin(ang)*K + (1 - cos(ang))*(K*K);
end
