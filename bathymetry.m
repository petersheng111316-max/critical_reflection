clear; close all; clc;

%% Load files and bathymetry data
bathy_data = load('Bathymetry_Profiles.mat');
x = bathy_data.dist_from_shore(:);
x_km = x / 1000;
paths = {
    'LJ_PL_depth', 'La Jolla to Point Loma', 'b';
    'Car_LJ_depth', 'Cardiff to La Jolla', 'r';
    'Car_PL_depth', 'Cardiff to Point Loma', 'g'
};
n2_data = load('Seasonal_N_2.mat');
N2_smooth = n2_data.N_2_smooth;

% Average across seasons to get a fixed N (used for plot)
N2_mean_profile = mean(N2_smooth, 2);
N2_mean_depth = mean(N2_mean_profile);
N_rad = sqrt(N2_mean_depth);
N_cph = N_rad * 3600 / (2 * pi); % Unit conversion

%% Input all parameters
w = 7.2921e-5; % Earth's rotation
latitude = 33.5;
f_rad = 2 * w * sin(deg2rad(latitude));
f_cph = f_rad * 3600 / (2 * pi); % convert into cph
M2_cph = 1.93 / 24; % Principal lunar semidiurnal constituent, M2
fprintf('The M2 frequency is %.4f cph\n', M2_cph);

% Give other important tidal frequencies for analysis
tidal_freqs = struct();
tidal_freqs.M2 = 1.93 / 24;        % Principal lunar semidiurnal (1.93 cpd)
tidal_freqs.S2 = 2.00 / 24;        % Principal solar semidiurnal
tidal_freqs.K1 = 1.00 / 24;        % Lunar-solar diurnal
tidal_freqs.O1 = 0.93 / 24;        % Lunar diurnal
tidal_freqs.N2 = 1.90 / 24;        % Larger lunar elliptic semidiurnal
tidal_freqs.near_inertial = f_cph; % Local inertial frequency

% Preallocation of struct
results = struct('name', cell(size(paths, 1), 1), ...
                 'raw_data', cell(size(paths, 1), 1), ...
                 'smooth_data', cell(size(paths, 1), 1), ...
                 's', cell(size(paths, 1), 1), ...
                 'sigma_c_cph', cell(size(paths, 1), 1), ...
                 'color', cell(size(paths, 1), 1));

for t = 1:size(paths, 1)
    depth_data = paths{t, 1};
    name = paths{t, 2};
    color = paths{t, 3};
    raw_data = bathy_data.(depth_data)(:);
    
    % Want depth negative and 0 at the surface
    if mean(raw_data) > 0
        raw_data = -raw_data;
    end
    
    smooth_window = max(round(length(raw_data)/20), 10);
    smooth_data = smoothdata(raw_data, 'movmean', smooth_window);
    
    % Define bottom slope s = -dD/dx from the paper
    s = -gradient(smooth_data, x);
    s = max(s, 0.0005);
    
    % Define critical frequency sigma_c
    sigma_c_rad = sqrt((N_rad^2 * s.^2 + f_rad^2) ./ (1 + s.^2));
    sigma_c_cph = sigma_c_rad * 3600 / (2 * pi);
    
    % Store results
    results(t).name = name;
    results(t).raw_data = raw_data;
    results(t).smooth_data = smooth_data;
    results(t).s = s;
    results(t).sigma_c_cph = sigma_c_cph;
    results(t).color = color;
end

%% Generate plots

% Bathymetry curves (depth versus distance from shore): raw + smooth
figure('Position', [50, 50, 1200, 900]);

subplot(2,1,1);
hold on;
for t = 1:length(results)
    plot(x_km, results(t).raw_data, results(t).color, 'LineWidth', 1.5);
end
xlabel('Distance from shore (km)', 'FontSize', 12);
ylabel('Depth (m)', 'FontSize', 12);
title('(a) Raw Bathymetry Data', 'FontSize', 12);
legend({results.name}, 'Location', 'best');
grid on;

subplot(2,1,2);
hold on;
for t = 1:length(results)
    plot(x_km, results(t).smooth_data, results(t).color, 'LineWidth', 2);
end
xlabel('Distance from shore (km)', 'FontSize', 12);
ylabel('Depth (m)', 'FontSize', 12);
title('(b) Smoothed Bathymetry Data', 'FontSize', 12);
legend({results.name}, 'Location', 'best');
grid on;

sgtitle(sprintf('Distance from Shore vs Depth for All Transects'), ...
    'FontSize', 14);

% Parameter sensitivity plot
t = 2; % Choose Cardiff to La Jolla

fprintf('\nTransect: %s\n', results(t).name);
fprintf('------------------------------------------------------------\n');
fprintf('%-20s %-12s %-12s %-20s\n', 'Constituent', 'Freq (cph)', ...
    'Freq (cpd)', 'Critical Distance (km)');
fprintf('%-20s %-12s %-12s %-20s\n', '----------', '-----------', ...
    '-----------', '---------------------');

% For each tidal frequency, find where sigma_c crosses that frequency
freq_names = fieldnames(tidal_freqs);
critical_locations = struct();

for i = 1:length(freq_names)
    fname = freq_names{i};
    fval = tidal_freqs.(fname);
    fval_cpd = fval * 24;  % convert to cycles per day for display
    
    % Find crossings where sigma_c equals this frequency
    sigma_c = results(t).sigma_c_cph;
    crossings = find(diff(sign(sigma_c - fval)) ~= 0);
    
    if isempty(crossings)
        % No crossings
        fprintf('No critical reflection at all frequencies for %s internal tides!\n', ...
            fname);
    else
        % Store crossing locations
        crossing_location = x_km(crossings);
        critical_locations.(fname) = crossing_location;
        for j = 1:length(crossing_location)
            fprintf('%-20s %-12.4f %-12.2f %-20.1f\n', ...
                fname, fval, fval_cpd, crossing_location(j));
        end
    end
end

% Comparison plot
% Choose Cardiff to La Jolla transect (second path)
t = 2;
figure('Position', [50, 50, 1200, 900]);

subplot(2,1,1);
plot(x_km, results(t).smooth_data, 'b-', 'LineWidth', 2);
ylabel('Depth (m)', 'FontSize', 12);
xlabel('Distance from shore (km)', 'FontSize', 12);
title(sprintf('Bathymetry: %s', results(t).name), 'FontSize', 12);
grid on;

subplot(2,1,2);
hold on;
plot(x_km, results(t).sigma_c_cph, 'r-', 'LineWidth', 2);

% Define colors for each tidal constituent
freq_names = fieldnames(tidal_freqs);
freq_colors = {'k', 'r', 'g', 'm', 'c', 'y'};

% Indicate all tidal frequencies (allowing rotation)
for i = 1:length(freq_names)
    fname = freq_names{i};
    fval = tidal_freqs.(fname);
    color = freq_colors{mod(i-1, length(freq_colors))+1}; 
    plot(x_km, fval * ones(size(x_km)), '--', ...
        'Color', color, 'LineWidth', 1);
end

ylabel('\sigma_c (cph)', 'FontSize', 12);
xlabel('Distance from shore (km)', 'FontSize', 12);
title('Critical Frequency vs Distance with Tidal Frequencies', 'FontSize', 12);
grid on;

% Add legend for bottom plot
legend_labels = {'\sigma_c', ...
    'M2 (1.93 cpd)', ...
    'S2 (2.00 cpd)', ...
    'K1 (1.00 cpd)', ...
    'O1 (0.93 cpd)', ...
    'N2 (1.90 cpd)', ...
    'near-inertial (1.10 cpd)'};
legend(legend_labels, 'Location', 'best', 'FontSize', 8);