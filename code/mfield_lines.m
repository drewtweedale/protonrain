function mfield_lines(B0, RN, NMag)
    % Create grid
    [x_grid, y_grid, z_grid] = meshgrid(-20:1:20, -20:1:20, -20:1:20);
    x_grid = x_grid * RN; y_grid = y_grid * RN; z_grid = z_grid * RN;
    % Compute field
    Bx_grid = zeros(size(x_grid));
    By_grid = zeros(size(y_grid));
    Bz_grid = zeros(size(z_grid));
    for i = 1:size(x_grid, 1)
        for j = 1:size(x_grid, 2)
            for k = 1:size(x_grid, 3)
                r = [x_grid(i,j,k), y_grid(i,j,k), z_grid(i,j,k)];
                phi = atan2(r(2), r(1));
                r_mag = norm(r);
                theta = acos(r(3)/r_mag);
                [Br, Btheta, Bphi] = new_comb(r, theta, phi, RN, NMag);
                B = sph2cart_field(Br, Btheta, Bphi, r);  % Get single output
                Bx_grid(i,j,k) = B(1);  % Assign components separately
                By_grid(i,j,k) = B(2);
                Bz_grid(i,j,k) = B(3);
            end
        end
    end

    % Plot
    figure(4); hold on;
    [x_neptune, y_neptune, z_neptune] = sphere;
    surf(x_neptune, y_neptune, z_neptune, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'blue');
    plot3([0 0], [0 0], [-3 3], 'Color', [0 0.5 0], 'LineWidth', 1);
    
    % Define starting points for magnetic field lines
    num_streamlines = 100;
    theta = linspace(0, 2*pi, num_streamlines);
    radius = 3 * RN;
    startx = radius * cos(theta);
    starty = radius * sin(theta);
    startz = .005 * RN * zeros(size(startx)); % Slightly above equator

    % Plot magnetic field lines
    streamline_color = 'r';

    for i = 1:num_streamlines
        h = streamline(x_grid / RN, y_grid / RN, z_grid / RN, Bx_grid, By_grid, Bz_grid, startx(i) / RN, starty(i) / RN, startz(i) / RN);
        g = streamline(x_grid / RN, y_grid / RN, z_grid / RN, -1 * Bx_grid, -1 * By_grid, -1 * Bz_grid, startx(i) / RN, starty(i) / RN, startz(i) / RN);
        hnegz = streamline(x_grid / RN, y_grid / RN, z_grid / RN, Bx_grid, By_grid, Bz_grid, startx(i) / RN, starty(i) / RN, -1 * startz(i) / RN);
        gnegz = streamline(x_grid / RN, y_grid / RN, z_grid / RN, -1 * Bx_grid, -1 * By_grid, -1 * Bz_grid, startx(i) / RN, starty(i) / RN, -1 * startz(i) / RN);
    
        % Set properties for the streamlines
        set(h, 'Color', streamline_color, 'LineWidth', 1);
        set(g, 'Color', streamline_color, 'LineWidth', 1);
        set(hnegz, 'Color', streamline_color, 'LineWidth', 1);
        set(gnegz, 'Color', streamline_color, 'LineWidth', 1);
    end
    
    xlabel('x (R_N)'); ylabel('y (R_N)'); zlabel('z (R_N)');
    title('Magnetic Field Lines'); 
    axis equal;
    view(3)
    grid on;
end