function boris_algo(B0, RN, q, m, pitch_angle_deg, KE_total, NMag)
    % Initial conditions
    current_R = 3 * RN;
    current_lon = 0;
    x_pos = current_R * cosd(current_lon);
    y_pos = current_R * sind(current_lon);
    z_pos = 0;  % Start in equatorial plane
            
    r0 = [x_pos, y_pos, z_pos];
    
    r_mag = norm(r0);
    theta = acos(r0(3)/r_mag);
    phi = atan2(r0(2), r0(1));
    [Br, Btheta, Bphi] = new_comb(r0, theta, phi, RN, 1);
    B_cart = sph2cart_field(Br, Btheta, Bphi, r0);
    B_dir = B_cart / norm(B_cart);  % Unit vector in direction of B
    if abs(B_dir(3)) < 0.9  % If B is not mostly in z-direction
        perp_dir = cross(B_dir, [0, 0, 1]);
    else
        perp_dir = cross(B_dir, [1, 0, 0]);
    end
    perp_dir = perp_dir / norm(perp_dir);

    % Initial velocity calculations
    pitch_angle = deg2rad(pitch_angle_deg);
    v_mag = sqrt(2 * KE_total / m);
    v_parallel = -v_mag * cos(pitch_angle) * B_dir;
    v_perp = v_mag * sin(pitch_angle) * perp_dir;
    v0 = v_parallel + v_perp;

    % Gyro period / Drift Period calculations
    T_g = 2*pi*m/(q*B0); 
    dt = T_g/20;
    t_end = 1000000*T_g;
    t = 0:dt:t_end; 
    n_steps = length(t);
    check_after = 100000; % Start checking for full orbit after this many steps
    check_after_bounce = 5000;
    check_after_boundary = 1000;
    
    % Arrays
    x = zeros(n_steps, 3); 
    v = zeros(n_steps, 3);
    E = [0,0,0];
    x(1,:) = r0; 
    v(1,:) = v0;
    hits = 0;

    % --- NEW: incremental/unwrapped longitude tracking ---
    initial_longitude = atan2d(r0(2), r0(1));
    prev_long = initial_longitude;                 % previous longitude (deg)
    cumulative_drift_signed = 0;                   % signed accumulated drift (deg)
    cumulative_drift = NaN;                        % final reported value (deg)
    full_orbit = false;
    lost = false;
    % ---------------------------------------------------

    % Add bounce tracking variables
    bounce_half_completed = false;
    bounce_return_point = [];
    bounce_left_return_point = false;
    bounce_check_active = false;
    bounce_count = 0;
    bounce_detected = false;
    
    % Variables for tracking boundary crossings
    entry_point = [];
    exit_point = [];
    in_boundary = false;
    
    % Integration
    for i = 1:n_steps-1
        
        x_mid = x(i, :) + 0.5 * dt * v(i, :);
        
        % Relativistic Correction
        v_curr = v(i,:);
        v_curr_sq = v_curr(1)^2 + v_curr(2)^2 + v_curr(3)^2;
        gamma = 1 / sqrt(1 - v_curr_sq / (3e8)^2);

        % Define theta and phi for magnetic field
        phi = atan2(x_mid(2), x_mid(1));
        r_mag = norm(x_mid);
        theta = acos(x_mid(3)/r_mag);
        
        % Field calculation
        [Br, Btheta, Bphi] = new_comb(x_mid, theta, phi, RN, NMag);
        B = sph2cart_field(Br, Btheta, Bphi, x_mid);

        % Boris push
        t_b = (q./(m.*gamma)) .* 0.5 .* dt .* B;
        s = 2 .* t_b ./ (1 + norm(t_b)^2);
        v_minus = v(i,:) + (q ./ m) .* E .* 0.5 .* dt;
        v_prime = v_minus + cross(v_minus,t_b);
        v_plus = v_minus + cross(v_prime,s);
        v(i+1,:) = v_plus + (q ./ m) .* E .* 0.5 .* dt;
        x(i+1,:) = x(i,:) + v(i+1,:) .* dt;
        
        % --- NEW: update incremental longitude and accumulate drift ---
        curr_long = atan2d(x(i+1,2), x(i+1,1));
        % signed smallest-angle step in (-180,180]
        step_delta = atan2d( sind(curr_long - prev_long), cosd(curr_long - prev_long) );
        cumulative_drift_signed = cumulative_drift_signed + step_delta;
        prev_long = curr_long;
        % --------------------------------------------------------------
        
        % To use the latest bounce return point as reference:
        if ~isempty(bounce_return_point)
            current_distance = norm(x(i+1,:) - bounce_return_point);
        else
            current_distance = norm(x(i+1,:) - r0);
        end

        current_iter = 0;
        
        % Detect entry into boundary
        if ~in_boundary && current_distance < RN * 0.1 && i > check_after_bounce
            in_boundary = true;
            entry_point = x(i+1,:);
            current_iter = i;
        end
        
        % Detect exit from boundary
        if in_boundary && current_distance > RN * 0.1 && (i - current_iter) > check_after_boundary
            in_boundary = false;
            exit_point = x(i+1,:);
            
            % Calculate average bounce point
            if ~isempty(entry_point) && ~isempty(exit_point)
                bounce_return_point = (entry_point + exit_point) / 2;
                bounce_count = bounce_count + 1;
                
                if bounce_count == 2
                    bounce_detected = true;
                end
                
                % NOTE: Do NOT recompute drift by comparing initial->avg_longitude.
                % We now compute drift incrementally above (cumulative_drift_signed).
                entry_point = [];
                exit_point = [];
            end
        end

        % Hit detection - check if particle hits Neptune's surface
        if norm(x(i+1,:)) <= RN
            hits = hits + 1;
            if ~bounce_detected
                fprintf('Particle hit before bouncing - no drift counted\n');
            end
            lost = true;
            break;
        end
        
        if norm(x(i+1,:)) >= 2 * current_R
            if ~bounce_detected
                cumulative_drift = 0;
                fprintf('Particle was lost via gradient intensity.\n');
            end
            lost = true;
            break;
        end


        % Check for full orbit completion (use accumulated value)
        if i > check_after && bounce_detected
            if abs(cumulative_drift_signed) > 360
                full_orbit = true;
                cumulative_drift = 360;  % sentinel for full orbit (as you used earlier)
                fprintf('Full orbit completed! Cumulative drift: %.2f°\n', cumulative_drift);
                break;
            end
        end
    end
    
    % If not full orbit, set final cumulative drift from accumulated sum
    if ~full_orbit
        if isnan(cumulative_drift)  % wasn't set to 0 on early loss
            cumulative_drift = abs(cumulative_drift_signed);
        end
    end

    % Plot
    figure(4); 
    hold on;
    [x_neptune, y_neptune, z_neptune] = sphere;
    surf(x_neptune, y_neptune, z_neptune, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'blue');
    plot3(x(:,1)/RN, x(:,2)/RN, x(:,3)/RN, 'r-', 'LineWidth', 1.5);
    plot3([0 0], [0 0], [-3 3], 'Color', [0 0.5 0], 'LineWidth', 1);
    xlabel('x (R_N)'); 
    ylabel('y (R_N)'); 
    zlabel('z (R_N)');
    title('Particle Trajectory'); 
    axis equal; 
    view(3); 
    grid on;
    fprintf('Particle collisions: %d\n', hits);
    fprintf('Final cumulative drift reported: %.2f°\n', cumulative_drift);
end
