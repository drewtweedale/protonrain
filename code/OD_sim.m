function OD_sim()
    % Parameters
    
    q_proton = 1.6e-19;       % (C)
    m_proton = 1.67e-27;      % (kg)
    RN = 24765e3;             % Neptune radius (m)
    energy_mev = 3;           % Fixed energy (MeV)
    KE_total = energy_mev * 1e6 * q_proton; % Total kinetic energy (J)
    
    % Radial distances and longitudes to test
    radial_distances = [1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3] * RN; 
    longitudes = 0:15:345;    % Degrees
    
    % Initialize results storage
    results = struct('radial_distance', num2cell(radial_distances/RN), ...
                    'longitudes', longitudes, ...
                    'loss_cones', zeros(length(longitudes), 1), ...
                    'drift_angles', zeros(length(longitudes), 1), ...
                    'loss_fac', zeros(length(longitudes), 1));
    
    % Main simulation loop
    fprintf('Starting simulation for %d MeV protons\n', energy_mev);
    fprintf('Testing %d positions (%.0f RN at %d longitudes)\n', ...
            length(radial_distances) * length(longitudes), ...
            radial_distances(1)/RN, length(longitudes));
    
    for r_idx = 1:length(radial_distances)
        current_R = radial_distances(r_idx);
        fprintf('\nTesting R = %.1f R_N:\n', current_R/RN);
        
        for lon_idx = 1:length(longitudes)
            current_lon = longitudes(lon_idx);
            fprintf('  Longitude %.1f°... ', current_lon);
            
            % Calculate initial positoin (magnetic equator)
            x_pos = current_R * cosd(current_lon);
            y_pos = current_R * sind(current_lon);
            z_pos = 0;  
            
            r0 = [x_pos, y_pos, z_pos];      

            % Loss cone finding
            found_bounce = false;
            test_pitch_angle = 0;
            max_test_angle = 90;
            loss_cone_deg = NaN;
            cumulative_drift = 0;

            while test_pitch_angle <= max_test_angle && ~found_bounce
                fprintf('    Testing pitch angle %.1f°... ', test_pitch_angle);
    
                % Simulation returns the total drift angle and a boolean to
                % show if the particle bounced.
                [cumulative_drift, bounced] = run_boris_simulation(r0, test_pitch_angle, KE_total, q_proton, m_proton, RN);
    
                if bounced
                    fprintf('Bounced (drift: %.2f°)\n', cumulative_drift);
                    found_bounce = true;
                    loss_cone_deg = test_pitch_angle - 1; % Loss cone is previous angle
                    if loss_cone_deg == -1
                        loss_cone_deg = 0;
                    end
                else
                    test_pitch_angle = test_pitch_angle + 1;
                end
            end

            if ~found_bounce
                fprintf('No bounce found up to 90°, using 90° as loss cone\n');
                loss_cone_deg = 90;
                % For this case, we still need to run the simulation to get drift angle
                cumulative_drift = run_boris_simulation(r0, 90, KE_total, q_proton, m_proton, RN);
            end

            fprintf('Pitch angle: %.2f° (loss cone: %.2f°)\n', loss_cone_deg + 1, loss_cone_deg);
            fprintf('    Cumulative drift: %.2f°\n', cumulative_drift);

            % Store both loss cone and drift angle
            results(r_idx).loss_cones(lon_idx) = loss_cone_deg;
            results(r_idx).drift_angles(lon_idx) = cumulative_drift;

            % Loss factor search
            fprintf('Performing Loss Factor Calculation\n');

            loss_fac = NaN;  
            step = 10;  
            current_pitch = loss_cone_deg + step;

            while found_bounce
                if current_pitch >= 90
                    loss_fac = 999;
                    fprintf('Particle never completes full orbit.');
                    break
                end
                fprintf('Testing pitch angle %.1f°\n', current_pitch);
                [cumulative_drift1, ~] = run_boris_simulation(r0, current_pitch, KE_total, q_proton, m_proton, RN);

                if cumulative_drift1 > 355
                    % Full orbit -> refine search downwards
                    fprintf('Full orbit at %.1f° -> refining downwards\n', current_pitch);

                    % Step down in 5° increments
                    test_pitch = current_pitch - 5;
                    while test_pitch >= loss_cone_deg
                        fprintf('  Testing %.1f°\n', test_pitch);
                        [cumulative_drift2, ~] = run_boris_simulation(r0, test_pitch, KE_total, q_proton, m_proton, RN);
                        if cumulative_drift2 <= 355
                            % Boundary found between test_pitch and test_pitch+5
                            % Refine with 1° increments
                            for fine_pitch = test_pitch+1 : current_pitch
                                fprintf('    Refining at %.1f°\n', fine_pitch);
                                [cumulative_drift3, ~] = run_boris_simulation(r0, fine_pitch, KE_total, q_proton, m_proton, RN);
                                if cumulative_drift3 > 355
                                    loss_fac = fine_pitch - loss_cone_deg;
                                    break;
                                end
                            end
                            break;
                        end
                        test_pitch = test_pitch - 5;
                    end

                    break;  % done refining
                else
                    % No full orbit -> move outward by +10°
                    current_pitch = current_pitch + step;
                end
            end

            % Store result
            if isnan(loss_fac)
                loss_fac = 0;  % default if no orbit found
            end
            results(r_idx).loss_fac(lon_idx) = loss_fac;
            fprintf('  ==> Loss Factor = %.1f°\n', loss_fac);
        end
    end
  
    % Create polar heat map
    create_polar_heatmap(radial_distances, longitudes, results, RN);
    create_losscone_heatmap(radial_distances, longitudes, results, RN);
    create_loss_fac_heatmap(radial_distances, longitudes, results, RN);
    % Save results to CSV
    save_results_to_csv(results, radial_distances, longitudes, RN);
end

function [cumulative_drift, bounced] = run_boris_simulation(r0, pitch_angle_deg, KE_total, q, m, RN)
    % Modified boris_algo function to return cumulative drift angle
    
    % Initial velocity calculations
    B0 = 1.42e-5;

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

    T_g = 2*pi*m/(q*B0); 
    dt = T_g/20;
    t_end = 5000000*T_g;
    t = 0:dt:t_end; 
    n_steps = length(t);
    check_after = 100000; % Start checking for full orbit after this many steps
    check_after_bounce = 5000;
    check_after_boundary = 1000;
    
    % Arrays for tracking
    x = zeros(n_steps, 3);
    v = zeros(n_steps, 3);
    E = [0, 0, 0];
    x(1,:) = r0;
    v(1,:) = v0;
    
    % ---incremental/unwrapped longitude tracking ---
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
        phi_mid = atan2(x_mid(2), x_mid(1));
        r_mag = norm(x_mid);
        theta_mid = acos(x_mid(3)/r_mag);
        
        % Field calculation
        [Br, Btheta, Bphi] = new_comb(x_mid, theta_mid, phi_mid, RN, 1);
        B = sph2cart_field(Br, Btheta, Bphi, x_mid);

        % Boris push
        t_b = (q./(m.*gamma)) .* 0.5 .* dt .* B;
        s = 2 .* t_b ./ (1 + norm(t_b)^2);
        v_minus = v(i,:) + (q ./ m) * E * 0.5 * dt;
        v_prime = v_minus + cross(v_minus,t_b);
        v_plus = v_minus + cross(v_prime,s);
        v(i+1,:) = v_plus + (q ./ m) * E * 0.5 * dt;
        x(i+1,:) = x(i,:) + v(i+1,:) * dt;
        
        % --- update incremental longitude and accumulate drift ---
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
            % Calculate average bounce point
            if ~isempty(entry_point) && ~isempty(exit_point)
                bounce_return_point = (entry_point + exit_point) / 2;
                bounce_count = bounce_count + 1;
                
                if bounce_count == 2
                    bounce_detected = true;
                end
                entry_point = [];
                exit_point = [];
            end
        end

        % Hit detection - check if particle hits Neptune's surface
        if norm(x(i+1,:)) <= RN
            if ~bounce_detected
                cumulative_drift = 0;
                fprintf('Particle hit before bouncing - no drift counted\n');
            end
            lost = true;
            break;
        end
        
        % Check for full orbit completion
        if i > check_after && bounce_detected
            if cumulative_drift > 355
                fprintf('Full orbit completed! Cumulative drift: %.2f°\n', cumulative_drift);
                break;
            end
        end

        if norm(x(i+1,:)) >= 2 * r_mag
            if ~bounce_detected
                cumulative_drift = 0;
            end
            fprintf('Particle was lost via gradient intensity.\n');
            cumulative_drift = 999;
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
    
    if ~lost && cumulative_drift ~= 360
        fprintf('Simulation ended without completion\n');
        cumulative_drift = NaN;
    end

    
    bounced = bounce_detected;
end

function create_polar_heatmap(radial_distances, longitudes, results, RN)
    % Create polar heat map visualization
    
    figure('Position', [100, 100, 800, 800], 'Name', 'Cumulative Drift Distance Heat Map');
    
    ax = polaraxes;  % Polar axes
    hold on;
    
    % Convert radial distances to RN units
    r_values = radial_distances / RN;
    
    % Reverse the hot colormap
    cmap = hot(256);
    cmap = flipud(cmap);
    colormap(ax, cmap);
    clim([0, 270]);
    
    % Plot each point
    for r_idx = 1:length(r_values)
        for lon_idx = 1:length(longitudes)
            drift = results(r_idx).drift_angles(lon_idx);
            original_lon = results(r_idx).longitudes(lon_idx);  % Get the original longitude
 
            if ~isnan(drift)
                r = r_values(r_idx);
                theta = deg2rad(original_lon);  % convert degrees to radians for polarplot
                
                if drift == 360
                    % Full orbit - white circle with black outline
                    polarscatter(theta, r, 300, 'w', 'filled', ...
                        'MarkerEdgeColor', 'k', 'LineWidth', 1);
                else
                    % Partial drift - colored circle with black outline
                    polarscatter(theta, r, 300, drift, 'filled', ...
                        'MarkerEdgeColor', 'k', 'LineWidth', 1);
                end
            end
        end
    end
    
    % Add Neptune at center (plot a big dot at r=0)
    polarscatter(0, 0, 27000, [0, 0, 0.5], 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    
    % Customize the plot
    colorbar;
    title(ax, sprintf('Cumulative Drift Distance Heat Map (%d MeV Protons)', 3));
    rlim([0, max(r_values) + 0.5]);

    ax.RTick = 1:floor(max(r_values));
    ax.RTickLabel{1} = '';
    
    hold off;
end

function create_losscone_heatmap(radial_distances, longitudes, results, RN)
    % Create polar heat map visualization
    
    figure('Position', [100, 100, 800, 800], 'Name', 'Loss Cone Angle Heatmap');
    
    ax = polaraxes;  % Polar axes
    hold on;
    
    % Convert radial distances to RN units
    r_values = radial_distances / RN;
    
    % Reverse the hot colormap
    cmap = hot(256);
    cmap = flipud(cmap);
    colormap(ax, cmap);
    clim([0, 90]);
    
    % Plot each point
    for r_idx = 1:length(r_values)
        for lon_idx = 1:length(longitudes)
            losscone = results(r_idx).loss_cones(lon_idx);
            original_lon = results(r_idx).longitudes(lon_idx);  % Get the original longitude
 
            if ~isnan(losscone)
                r = r_values(r_idx);
                theta = deg2rad(original_lon);  % convert degrees to radians for polarplot
                % Partial drift - colored circle with black outline
                polarscatter(theta, r, 300, losscone, 'filled', ...
                    'MarkerEdgeColor', 'k', 'LineWidth', 1);
            end
        end
    end
    
    % Add Neptune at center (plot a big dot at r=0)
    polarscatter(0, 0, 27000, [0, 0, 0.5], 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    
    % Customize the plot
    colorbar;
    title(ax, sprintf('Loss Cone Heat Map (%d MeV Protons)', 3));
    rlim([0, max(r_values) + 0.5]);

    ax.RTick = 1:floor(max(r_values));
    ax.RTickLabel{1} = '';
    
    hold off;
end

function save_results_to_csv(results, radial_distances, longitudes, RN)
    % save_results_to_csv - Saves simulation results to a CSV file.
    %
    %   Columns:
    %   1. Starting longitude (deg)
    %   2. Starting radial distance (R_N)
    %   3. Loss cone angle (deg)
    %   4. Cumulative drift angle (deg)

    % Preallocate arrays
    n_r = length(radial_distances);
    n_lon = length(longitudes);
    n_total = n_r * n_lon;

    data = zeros(n_total, 5);  % (longitude, r, loss cone, drift)

    idx = 1;
    for r_idx = 1:n_r
        for lon_idx = 1:n_lon
            data(idx,1) = results(r_idx).longitudes(lon_idx);              % Longitude
            data(idx,2) = results(r_idx).radial_distance;                  % R (in R_N)
            data(idx,3) = results(r_idx).loss_cones(lon_idx);              % Loss cone
            data(idx,4) = results(r_idx).drift_angles(lon_idx);            % Drift angle
            data(idx,5) = results(r_idx).loss_fac(lon_idx);            % Drift angle
            idx = idx + 1;
        end
    end

    % Convert to table for nice column headers
    T = array2table(data, ...
        'VariableNames', {'Longitude_deg','RadialDistance_RN','LossCone_deg','DriftAngle_deg', 'LossFac_deg'});

    % Create filename with timestamp
    timestamp = datestr(now,'yyyymmdd_HHMMSS');
    filename = sprintf('OD_sim_results_%s.csv', timestamp);

    % Write to CSV in current directory
    writetable(T, filename);

    fprintf('Results saved to %s\n', filename);
end

function create_loss_fac_heatmap(radial_distances, longitudes, results, RN)
    % Create polar heat map visualization
    
    figure('Position', [100, 100, 800, 800], 'Name', 'Loss Factor Heatmap');
    
    ax = polaraxes;  % Polar axes
    hold on;
    
    % Convert radial distances to RN units
    r_values = radial_distances / RN;
    
    % Reverse the hot colormap
    cmap = hot(256);
    cmap = flipud(cmap);
    colormap(ax, cmap);
    clim([0, 90]);
    
    % Plot each point
    for r_idx = 1:length(r_values)
        for lon_idx = 1:length(longitudes)
            loss_fac = results(r_idx).loss_fac(lon_idx);
            original_lon = results(r_idx).longitudes(lon_idx);  % Get the original longitude
 
            if ~isnan(loss_fac)
                r = r_values(r_idx);
                theta = deg2rad(original_lon);  % convert degrees to radians for polarplot
                % Partial drift - colored circle with black outline
                if loss_fac == 999
                    polarscatter(theta, r, 300, 'w', 'filled', ...
                    'MarkerEdgeColor', 'k', 'LineWidth', 1);
                    % Overlay an X at the same location
                    polarscatter(theta, r, 300, 'kx', 'LineWidth', 2); 
                else
                    polarscatter(theta, r, 300, loss_fac, 'filled', ...
                    'MarkerEdgeColor', 'k', 'LineWidth', 1);
                end
            end
        end
    end
    
    % Add Neptune at center
    polarscatter(0, 0, 27000, [0, 0, 0.5], 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    
    % Customize the plot
    colorbar;
    title(ax, sprintf('Loss Factor Map (%d MeV Protons)', 3));
    rlim([0, max(r_values) + 0.5]);

    ax.RTick = 1:floor(max(r_values));
    ax.RTickLabel{1} = '';
    
    hold off;
end
