format long

% range of parameters
A_values = linspace(0, 35, 36);
Phi_m_values = linspace(1e-10, 1e-6, 50); 
d_values = linspace(0, 15, 50);
% Tumor volume without treatment 
tumor_volume_no_treatment = 1.85e9;

% Time span 
t_span = [0, 22];
y0 = [1e6, 1e6,  1e6]; 

% Initialize the efficacy matrix
efficacy_matrix = zeros(length(Phi_m_values), length(A_values), length(d_values));

% Loop over all combinations 
for i = 1:length(Phi_m_values)
    for j = 1:length(A_values)
        for k = 1:length(d_values)
            A = A_values(j);
            Phi_m = Phi_m_values(i);
            d = d_values(k);
            
            % Solve 
            [~, y] = ode23s(@(t, y) g(t, y, A, Phi_m, d), t_span, y0);
            
            %  tumor volume with treatment at the last time point (day 22)
            tumor_volume_with_treatment = y(end, 1) + y(end, 2) + y(end, 3);
            
            % Calculate efficacy
            efficacy = ((tumor_volume_no_treatment - tumor_volume_with_treatment) / tumor_volume_no_treatment) * 100;
            
            % Store efficacy in the matrix
            efficacy_matrix(i, j, k) = efficacy;
        end
    end
end

% original 3D plot
[A_mesh, Phi_m_mesh, d_mesh] = meshgrid(A_values, Phi_m_values, d_values);
figure;
hold on;

% Create isosurfaces
p1 = patch(isosurface(A_mesh, Phi_m_mesh, d_mesh, efficacy_matrix, 25));
set(p1, 'FaceColor', [0.5 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.7); % Dark Red

p2 = patch(isosurface(A_mesh, Phi_m_mesh, d_mesh, efficacy_matrix, 50));
set(p2, 'FaceColor', [0.5 0 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.6); % Dark Purple

p3 = patch(isosurface(A_mesh, Phi_m_mesh, d_mesh, efficacy_matrix, 75));
set(p3, 'FaceColor', [0.9 0.9 0], 'EdgeColor', 'none', 'FaceAlpha', 0.5); % Dark Yellow

p4 = patch(isosurface(A_mesh, Phi_m_mesh, d_mesh, efficacy_matrix, 95));
set(p4, 'FaceColor', [0 0.5 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4); % Dark Green

% y-slice values for multiple layers
yslice_values = [7.46e-7, 9.07e-7, 4.34e-8, 1.43e-8, 9.63e-9];

% Define light colors for each y-slice
slice_colors = [
    0.9 0.2 0.2;     % Light Red
    0.2 0.2 0.9;     % Light Blue
    0.2 0.9 0.2;     % Light Green
    0.9 0.6 0.2;   % Light Orange
    0.6 0.2 0.9;   % Light Purple
];

% Interpolate and plot each y-slice as a surface with light colors
for idx = 1:length(yslice_values)
    yslice = yslice_values(idx);
    
    % Find the index closest to the y-slice value in Phi_m_values
    [~, slice_idx] = min(abs(Phi_m_values - yslice));
    
    % Extract the corresponding slice from the efficacy matrix
    slice_data = squeeze(efficacy_matrix(slice_idx, :, :));
    
    % Create a 2D meshgrid for the slice to plot as a surface
    [A_grid, D_grid] = meshgrid(A_values, d_values);
    
    % Plot the surface for the y-slice with the assigned light color
    surf(A_grid, yslice * ones(size(A_grid)), D_grid, slice_data', ...
         'EdgeColor', 'none', 'FaceAlpha', 0.3, 'FaceColor', slice_colors(idx, :));
end

%
view(3);


xlabel('A', 'FontSize', 18);
zlabel('Dose (Gy)', 'FontSize', 18);
ylabel('\phi_m', 'FontSize', 18);


title('3D Efficacy Cuboid with Highlighted Points', 'FontSize', 14);

% Add grid
grid on;
camlight headlight;
lighting phong; % Smoother shading model
alpha(0.8); 
% Add legend
legend('Efficacy Level 25', 'Efficacy Level 50', 'Efficacy Level 75', 'Efficacy Level 95', 'Location', 'best');

hold off;

% Create a new figure for 2D projections
figure;
for idx = 1:length(yslice_values)
    yslice = yslice_values(idx);
    
    % Find the index closest to the y-slice value in Phi_m_values
    [~, slice_idx] = min(abs(Phi_m_values - yslice));
    
    % Extract the corresponding slice from the efficacy matrix
    slice_data = squeeze(efficacy_matrix(slice_idx, :, :));
    
    % Create a new subplot
    subplot(1, length(yslice_values), idx);
    
    % Set background color to match the Phi_m slice color
    set(gca, 'Color', slice_colors(idx, :));
    
    % Hold the plot to add multiple contour levels
    hold on;
    
    % Define contour levels and corresponding dark colors
    efficacy_levels = [25, 50, 75, 95];
    contour_colors = [
        0.7 0 0;     % Dark Red for 25%
        0 0 0.7;     % Dark Blue for 50%
        0 0.7 0;     % Dark Green for 75%
        0.7 0.5 0;   % Dark Orange for 95%
    ];
    
    % Efficacy 
    for c_idx = 1:length(efficacy_levels)
        contour(A_values, d_values, slice_data', [efficacy_levels(c_idx) efficacy_levels(c_idx)], ...
                'LineWidth', 1.5, 'LineColor', contour_colors(c_idx, :));
    end
    
 
    title(['\phi_m = ' num2str(yslice)], 'FontSize', 12);
    xlabel('A', 'FontSize', 10);
    ylabel('Dose (Gy)', 'FontSize', 10);
    grid on;
    
    hold off;
end

function dydt =  g(t, y, A, Phi_m, d)
    dydt = zeros(3,1);
        t0 = d/1.2;
        if (d==0)
            u=1;
        else
            u = 2 * (2.56 *t0 + exp(-2.56 *t0) - 1) / ((  2.56^2)*t0^2);
        end
       S=exp(-5.15e-3 *d -6.58e-3*d*d*u);
    
    dydt(1) = 0.453 *y(1)*(1 - y(1)/1.77e9) -4.5e-08 *y(1)*y(2)-Phi_m*y(1)*y(3);
    dydt(2) = 1.65e4-7.41e-2*y(2)- 1.98e-9*y(1)*y(2);
    dydt(3) =4.39e-7   -3.949e-1 *y(3)-8.15338e-7*y(1)*y(3);
    
    if (t >= 12)
        dydt(1) = dydt(1)*S-A*(1-S)*dydt(1); 
    else
        dydt(1) = dydt(1); 
    end
    end