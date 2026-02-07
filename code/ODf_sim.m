function ODf_sim() 
    csv_file = "OD_sim_results_20251001_120034.csv";
    data = load_csv(csv_file);
    
    % Parameters
    q_proton = 1.6e-19;       % (C)
    m_proton = 1.67e-27;      % (kg)
    RN = 24765e3;             % Neptune radius (m)
    energy_mev = 3;           % Fixed energy (MeV)
    KE_total = energy_mev * 1e6 * q_proton; % Total kinetic energy (J)
    
    % Initialize arrays to store precipitation locations
    precip_longitudes = [];
    precip_latitudes = [];
    precip_colatitudes = [];
    start_positions = [];
    
    fprintf('Tracking particle precipitation locations...\n');
    
    % Process each starting position from the CSV
    for i = 1:length(data.longitudes)
        current_lon = data.longitudes(i);
        current_R = data.radial_distances(i) * RN; % Convert RN to meters
        loss_cone = data.loss_cones(i);
        
        fprintf('Processing position %d/%d: Lon=%.1f°, R=%.2f RN, LossCone=%.1f°\n', ...
                i, length(data.longitudes), current_lon, current_R/RN, loss_cone);
        
        % Calculate initial position (magnetic equator)
        x_pos = current_R * cosd(current_lon);
        y_pos = current_R * sind(current_lon);
        z_pos = 0;  
        
        r0 = [x_pos, y_pos, z_pos];
        
        % Initialize particle 1 degree outside loss cone
        pitch_angle = loss_cone + 1;
        
        % Run simulation to track precipitation location
        [precip_location, lost] = track_precipitation_location(r0, pitch_angle, KE_total, q_proton, m_proton, RN);
        
        if lost && ~isempty(precip_location)
            % Convert precipitation location to spherical coordinates
            x_precip = precip_location(1);
            y_precip = precip_location(2);
            z_precip = precip_location(3);
            
            % Calculate spherical coordinates
            r_precip = norm(precip_location);
            precip_colatitude = acosd(z_precip / r_precip); % 0° at north pole, 180° at south pole
            precip_latitude = 90 - precip_colatitude; % -90° to +90°
            precip_longitude = atan2d(y_precip, x_precip);
            if precip_longitude < 0
                precip_longitude = precip_longitude + 360; % Convert to 0-360 range
            end
            
            % Store results
            precip_longitudes = [precip_longitudes; precip_longitude];
            precip_latitudes = [precip_latitudes; precip_latitude];
            precip_colatitudes = [precip_colatitudes; precip_colatitude];
            start_positions = [start_positions; current_lon, current_R/RN, loss_cone];
            
            fprintf('  Particle precipitated at: Lon=%.1f°, Lat=%.1f°\n', ...
                    precip_longitude, precip_latitude);
        else
            fprintf('  Particle did not precipitate\n');
        end
    end
    
    % Create magnetic field heatmap with precipitation points
    create_precipitation_heatmap(RN, precip_longitudes, precip_latitudes, start_positions);
    
    % Save precipitation data to CSV
    save_precipitation_data(precip_longitudes, precip_latitudes, precip_colatitudes, start_positions);
end

function [precip_location, lost] = track_precipitation_location(r0, pitch_angle_deg, KE_total, q, m, RN)
    % Track particle until it precipitates and return the location
    
    % Initial velocity calculations
    B0 = 1.42e-5;

    r_mag = norm(r0);
    theta = acos(r0(3)/r_mag);
    phi = atan2(r0(2), r0(1));
    [Br, Btheta, Bphi] = new_comb(r0, theta, phi, RN, 1);
    B_cart = sph2cart_field(Br, Btheta, Bphi, r0);
    B_dir = B_cart / norm(B_cart);
    if abs(B_dir(3)) < 0.9
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
    
    % Arrays for tracking
    x = zeros(n_steps, 3);
    v = zeros(n_steps, 3);
    E = [0, 0, 0];
    x(1,:) = r0;
    v(1,:) = v0;
    
    lost = false;
    precip_location = [];
    
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
        
        % Check for precipitation (hit Neptune's surface)
        if norm(x(i+1,:)) <= RN
            lost = true;
            precip_location = x(i+1,:);
            break;
        end
        
        % Early termination if particle escapes
        if norm(x(i+1,:)) >= 5 * RN
            break;
        end
    end
end

function create_precipitation_heatmap(RN, precip_longitudes, precip_latitudes, start_positions)
    % Create magnetic field heatmap with precipitation points overlaid
    
    % Define grid of longitudes and latitudes
    longitudes = 0:0.5:355;
    latitudes = -90:0.5:90;
    
    % Initialize magnetic field strength matrix
    B_strength = zeros(length(latitudes), length(longitudes));
    
    fprintf('Creating precipitation heatmap...\n');
    
    % Calculate magnetic field strength at each point
    for lat_idx = 1:length(latitudes)
        for lon_idx = 1:length(longitudes)
            latitude = latitudes(lat_idx);
            colatitude = 90 - latitude;
            
            theta = deg2rad(colatitude);
            phi = deg2rad(longitudes(lon_idx));
            
            x = RN * sin(theta) * cos(phi);
            y = RN * sin(theta) * sin(phi);
            z = RN * cos(theta);
            r = [x; y; z];
            
            [Br, Btheta, Bphi] = new_comb(r, theta, phi, RN, 1);
            B_mag = sqrt(Br^2 + Btheta^2 + Bphi^2) * 1e6;
            
            B_strength(lat_idx, lon_idx) = B_mag;
        end
    end
    
    % Create the heatmap figure
    figure('Position', [100, 100, 1200, 700], 'Name', 'Neptune Precipitation Map');
    
    % Create heatmap
    imagesc(longitudes, latitudes, B_strength);
    colormap(flipud(hot));
    colorbar;
    c = colorbar;
    c.Label.String = 'Magnetic Field Strength (μT)';
    
    hold on;
    
    % Plot precipitation points
    if ~isempty(precip_longitudes)
        scatter(precip_longitudes, precip_latitudes, 50, 'cyan', 'filled', ...
                'MarkerEdgeColor', 'blue', 'LineWidth', 1.5);
        
        % Add labels for some points to show starting positions
        num_labels = min(10, length(precip_longitudes));
        for i = 1:num_labels
            text(precip_longitudes(i), precip_latitudes(i) + 2, ...
                 sprintf('Start: %.0f°', start_positions(i,1)), ...
                 'Color', 'white', 'FontSize', 8, 'HorizontalAlignment', 'center');
        end
    end
    
    hold off;
    
    % Set axis labels and title
    xlabel('Longitude (degrees)');
    ylabel('Latitude (degrees)');
    title('Neptune Magnetic Field Strength with Particle Precipitation Locations');
    
    % Set axis properties
    axis xy;
    set(gca, 'XTick', 0:30:360);
    set(gca, 'YTick', -90:30:90);
    grid on;
    
    % Add legend
    if ~isempty(precip_longitudes)
        legend('Precipitation Points', 'Location', 'best');
    end
    
    fprintf('Precipitation heatmap completed. %d precipitation events plotted.\n', ...
            length(precip_longitudes));
end

function save_precipitation_data(precip_longitudes, precip_latitudes, precip_colatitudes, start_positions)
    % Save precipitation data to CSV file
    
    % Create data table
    n_points = length(precip_longitudes);
    data = zeros(n_points, 6);
    
    for i = 1:n_points
        data(i,1) = start_positions(i,1);  % Start longitude
        data(i,2) = start_positions(i,2);  % Start radial distance (R_N)
        data(i,3) = start_positions(i,3);  % Loss cone
        data(i,4) = precip_longitudes(i);  % Precipitation longitude
        data(i,5) = precip_latitudes(i);   % Precipitation latitude
        data(i,6) = precip_colatitudes(i); % Precipitation colatitude
    end
    
    % Convert to table
    T = array2table(data, ...
        'VariableNames', {'StartLongitude_deg', 'StartRadialDistance_RN', 'LossCone_deg', ...
                         'PrecipLongitude_deg', 'PrecipLatitude_deg', 'PrecipColatitude_deg'});
    
    % Create filename with timestamp
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    filename = sprintf('precipitation_data_%s.csv', timestamp);
    
    % Write to CSV
    writetable(T, filename);
    
    fprintf('Precipitation data saved to %s\n', filename);
end