function [Br, Btheta, Bphi] = new_comb(r, theta, phi, RN, NMag)
    r = r(:);
    rmag = norm(r);
    % Keep the NLS angles to project the final B back
    theta_nls = theta;
    phi_nls   = phi;

    % =========================
    %  Schmidt (Gauss) Coeffs
    % =========================
    sf = 1e-4;
    
    % dipole (n=1)
    g10 =  0.09732 * sf;
    g11 =  0.03220 * sf;
    h11 = -0.09889 * sf;

    % quadrupole (n=2)
    g20 =  0.07448 * sf;
    g21 =  0.00664 * sf;
    h21 =  0.11230 * sf;
    g22 =  0.04499 * sf;
    h22 = -0.00070 * sf;

    % octupole (n=3)
    g30 = -0.06592 * sf;
    g31 =  0.04098 * sf;
    h31 = -0.03669 * sf;
    g32 = -0.03581 * sf;
    h32 =  0.01791 * sf;
    g33 =  0.00484 * sf;
    h33 = -0.00770 * sf;

    % =========================
    %  Rodrigues rotation
    % =========================
    use_rotation = NMag;
    if use_rotation
        tilt = deg2rad(46.9);         % colatitude from +Z
        lonW = deg2rad(252);          % west longitude
        lonE = -lonW;                 % convert to east-positive azimuth

        % Dipole axis in NLS frame
        m_axis = [sin(tilt)*cos(lonE);
                  sin(tilt)*sin(lonE);
                  cos(tilt)];

        % R such that R*m_axis = +z_hat (aligns dipole with +z)
        R = rodrigues_align(m_axis, [0;0;1]);

        % Rotate position into dipole frame
        r_work = R * r;
        rmag_work = norm(r_work);    

        % Spherical angles in dipole frame (radians)
        theta_work = acos(r_work(3)/rmag_work);         % colatitude
        phi_work   = atan2(r_work(2), r_work(1));       % azimuth (east+)
    else
        R = eye(3);
        r_work = r;
        rmag_work = rmag;
        theta_work = theta_nls;
        phi_work   = phi_nls;
    end

    % ======================================
    %  Associated Legendre 
    % ======================================
    P  = zeros(4,4);
    dP = zeros(4,4);
    
    % dipole (n=1)
    P(2,1) = cos(theta_work);
    P(2,2) = sin(theta_work);
    
    % quadrupole (n=2)
    P(3,1) = (3/2) * (cos(theta_work)^2 - 1/3);
    P(3,2) = sqrt(3) * cos(theta_work) * sin(theta_work);
    P(3,3) = (sqrt(3)/2) * sin(theta_work)^2;
    
    % octupole (n=3)
    P(4,1) = (5/2) * cos(theta_work) * (cos(theta_work)^2 - 9/15);
    P(4,2) = (5*sqrt(3))/(2*sqrt(2)) * sin(theta_work) * (cos(theta_work)^2 - 3/15);
    P(4,3) = (sqrt(3)*sqrt(5))/2 * cos(theta_work) * sin(theta_work)^2;
    P(4,4) = sqrt(5)/(2*sqrt(2)) * sin(theta_work)^3;
    
    % derivatives with respect to theta
    if sin(theta_work) > 1e-10
        % dipole (n=1)
        dP(2,1) = 1.0;
        dP(2,2) = -1/tan(theta_work);
        
        % quadrupole (n=2)
        dP(3,1) = 3*cos(theta_work);
        dP(3,2) = -sqrt(3)/sin(theta_work) * (cos(theta_work)^2 - sin(theta_work)^2);
        dP(3,3) = -sqrt(3) * cos(theta_work);
        
        % octupole (n=3)
        dP(4,1) = (3/2) * (5*cos(theta_work)^2 - 1);
        dP(4,2) = -(5*sqrt(3))/(2*sqrt(2)) * (1/tan(theta_work)) * (-2*sin(theta_work)^2 + cos(theta_work)^2 - 1/5);
        dP(4,3) = sqrt(15)/2 * (sin(theta_work)^2 - 2*cos(theta_work)^2);
        dP(4,4) = -(3/2) * sqrt(5/2) * sin(theta_work) * cos(theta_work);
    else
        % at poles, derivatives = zero
        dP(2,1) = 0;
        dP(2,2) = 0;
        dP(3,1) = 0;
        dP(3,2) = 0;
        dP(3,3) = 0;
        dP(4,1) = 0;
        dP(4,2) = 0;
        dP(4,3) = 0;
        dP(4,4) = 0;
    end

    % ======================================
    %  Magnetic field calculation
    % ======================================
    Br_w = 0.0;
    Bth_w = 0.0;
    Bph_w = 0.0;
    
    at_pole = abs(sin(theta_work)) < 1e-10;
    
    % DIPOLE TERMS (n=1)
    n = 1;
    
    % Br
    Br_factor_n1 = (n+1) * (RN/rmag_work)^(n+2);
    % m=0
    Br_w = Br_w + Br_factor_n1 * g10 * P(2,1);
    % m=1
    Br_w = Br_w + Br_factor_n1 * (g11*cos(1*phi_work) + h11*sin(1*phi_work)) * P(2,2);
    
    % Btheta
    Btheta_factor_n1 = (RN/rmag_work)^(n+2);
    % m=0
    Bth_w = Bth_w + Btheta_factor_n1 * sin(theta_work) * g10 * dP(2,1);
    % m=1
    Bth_w = Bth_w + Btheta_factor_n1 * sin(theta_work) * (g11*cos(1*phi_work) + h11*sin(1*phi_work)) * dP(2,2);
    
    % Bphi
    if ~at_pole
        Bphi_factor_n1 = (RN/rmag_work)^(n+2) / sin(theta_work);
        % m=1
        m = 1;
        Bph_w = Bph_w + Bphi_factor_n1 * m * (g11*sin(m*phi_work) - h11*cos(m*phi_work)) * P(2,2);
    end
    
    % QUADRUPOLE TERMS (n=2)
    n = 2;
    
    % Br
    Br_factor_n2 = (n+1) * (RN/rmag_work)^(n+2);
    % m=0
    Br_w = Br_w + Br_factor_n2 * g20 * P(3,1);
    % m=1
    Br_w = Br_w + Br_factor_n2 * (g21*cos(1*phi_work) + h21*sin(1*phi_work)) * P(3,2);
    % m=2
    Br_w = Br_w + Br_factor_n2 * (g22*cos(2*phi_work) + h22*sin(2*phi_work)) * P(3,3);
    
    % Btheta
    Btheta_factor_n2 = (RN/rmag_work)^(n+2);
    % m=0
    Bth_w = Bth_w + Btheta_factor_n2 * sin(theta_work) * g20 * dP(3,1);
    % m=1
    Bth_w = Bth_w + Btheta_factor_n2 * sin(theta_work) * (g21*cos(1*phi_work) + h21*sin(1*phi_work)) * dP(3,2);
    % m=2
    Bth_w = Bth_w + Btheta_factor_n2 * sin(theta_work) * (g22*cos(2*phi_work) + h22*sin(2*phi_work)) * dP(3,3);
    
    % Bphi
    if ~at_pole
        Bphi_factor_n2 = (RN/rmag_work)^(n+2) / sin(theta_work);
        % m=1
        m = 1;
        Bph_w = Bph_w + Bphi_factor_n2 * m * (g21*sin(m*phi_work) - h21*cos(m*phi_work)) * P(3,2);
        % m=2
        m = 2;
        Bph_w = Bph_w + Bphi_factor_n2 * m * (g22*sin(m*phi_work) - h22*cos(m*phi_work)) * P(3,3);
    end
    
    % OCTUPOLE TERMS (n=3)
    n = 3;
    
    % Br
    Br_factor_n3 = (n+1) * (RN/rmag_work)^(n+2);
    % m=0
    Br_w = Br_w + Br_factor_n3 * g30 * P(4,1);
    % m=1
    Br_w = Br_w + Br_factor_n3 * (g31*cos(1*phi_work) + h31*sin(1*phi_work)) * P(4,2);
    % m=2
    Br_w = Br_w + Br_factor_n3 * (g32*cos(2*phi_work) + h32*sin(2*phi_work)) * P(4,3);
    % m=3
    Br_w = Br_w + Br_factor_n3 * (g33*cos(3*phi_work) + h33*sin(3*phi_work)) * P(4,4);
    
    % Btheta
    Btheta_factor_n3 = (RN/rmag_work)^(n+2);
    % m=0
    Bth_w = Bth_w + Btheta_factor_n3 * sin(theta_work) * g30 * dP(4,1);
    % m=1
    Bth_w = Bth_w + Btheta_factor_n3 * sin(theta_work) * (g31*cos(1*phi_work) + h31*sin(1*phi_work)) * dP(4,2);
    % m=2
    Bth_w = Bth_w + Btheta_factor_n3 * sin(theta_work) * (g32*cos(2*phi_work) + h32*sin(2*phi_work)) * dP(4,3);
    % m=3
    Bth_w = Bth_w + Btheta_factor_n3 * sin(theta_work) * (g33*cos(3*phi_work) + h33*sin(3*phi_work)) * dP(4,4);
    
    % Bphi
    if ~at_pole
        Bphi_factor_n3 = (RN/rmag_work)^(n+2) / sin(theta_work);
        % m=1
        m = 1;
        Bph_w = Bph_w + Bphi_factor_n3 * m * (g31*sin(m*phi_work) - h31*cos(m*phi_work)) * P(4,2);
        % m=2
        m = 2;
        Bph_w = Bph_w + Bphi_factor_n3 * m * (g32*sin(m*phi_work) - h32*cos(m*phi_work)) * P(4,3);
        % m=3
        m = 3;
        Bph_w = Bph_w + Bphi_factor_n3 * m * (g33*sin(m*phi_work) - h33*cos(m*phi_work)) * P(4,4);
    end

    % ======================================
    %  Rotate B back to NLS and project
    % ======================================
    if use_rotation
        % Convert working-frame spherical field to Cartesian (dipole frame),
        % rotate back to NLS, then project onto caller's (theta_nls, phi_nls).
        B_cart_work = sph2cart_field(Br_w, Bth_w, Bph_w, r_work.');  % [Bx,By,Bz]
        B_cart_nls  = R.' * B_cart_work(:);                          % rotate back
        B_sph_nls   = cart2sph_field(B_cart_nls(1), B_cart_nls(2), B_cart_nls(3), ...
                                     theta_nls, phi_nls);
        Br     = B_sph_nls(1);
        Btheta = B_sph_nls(2);
        Bphi   = B_sph_nls(3);
    else
        Br     = Br_w;
        Btheta = Bth_w;
        Bphi   = Bph_w;
    end
end