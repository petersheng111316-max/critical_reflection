clear; close all; clc;

%% Input all constants and Parameters (from Eriksen 1982, Muir Seamount)

% Earth's rotation rate (rad/s)
w = 7.2921e-5;

% Parameters from Eriksen 1982 + buoyancy(N)/Coriolis(F) frequencies
latitude = 33.5;  % degrees North
f = 2 * w * sin(deg2rad(latitude));
f_cph = f * 3600 / (2 * pi);  % Convert to cycles per hour (cph)
N_cph = 0.605;  % cph
N = N_cph * 2 * pi / 3600;  % Convert to rad/s

% Bottom slope s=tan(alpha)
s = 0.3;
alpha = atan(s);  % slope angle in radians

% Compute critical frequency from formula
sigma_c_rad = sqrt((N^2 * s^2 + f^2) / (1 + s^2));
sigma_c_cph = sigma_c_rad * 3600 / (2 * pi); % Convert to cph

fprintf('============================================================\n');
fprintf('Parameters (Muir Seamount, from Eriksen 1982)\n');
fprintf('============================================================\n');
fprintf('Latitude: %.1f°N\n', latitude);
fprintf('Coriolis frequency f: %.4f cph\n', f_cph);
fprintf('Buoyancy frequency N: %.3f cph\n', N_cph);
fprintf('Critical frequency sigma_c: %.3f cph\n', sigma_c_cph);
fprintf('============================================================\n');

% Frequency range (cph)
sigma_cph = logspace(log10(1.01 * f_cph), log10(0.99 * N_cph), 200);
sigma_rad = sigma_cph * 2 * pi / 3600;

% Compute wavenumber ratio and energy density for 2D case
m_ratio_2d = NaN(size(sigma_rad));
theta_i_vals = NaN(size(sigma_rad));

for i = 1:length(sigma_rad)
    sig = sigma_rad(i);
    [theta, isvalid] = compute_theta(sig, f, N);
    if isvalid
        theta_i_vals(i) = theta;
        m_ratio_2d(i) = wavenumber(theta, alpha);
    end
end

% Convert to dB according to Eriksen
E_ratio = m_ratio_2d.^2;
E_total = 1 + E_ratio;
E_total_dB = 10 * log10(E_total);

%% Generate all plots for visualization

% Figure 1: Total Energy vs Frequency (2D)
figure('Position', [100, 100, 700, 500]);
semilogx(sigma_cph, E_total_dB, 'b-', 'LineWidth', 2);
hold on;
plot([sigma_c_cph, sigma_c_cph], ylim, 'r--', 'LineWidth', 1.5);
plot([f_cph, f_cph], ylim, 'k:', 'LineWidth', 1);
plot([N_cph, N_cph], ylim, 'k:', 'LineWidth', 1);
xlabel('Frequency (cph)', 'FontSize', 12);
ylabel('10 log_{10}(E_{total}) (dB)', 'FontSize', 12);
title('Total Energy vs Frequency for Single Wave Reflection', ...
    'FontSize', 12);
legend('Trend', '\sigma_c', 'f', 'N', 'Location', 'best');
grid on;
xlim([0.9 * f_cph, 1.1 * N_cph]); % no trend below f or above N

% Figure 2: Wavenumber Amplification vs Frequency (2D)
figure('Position', [100, 100, 700, 500]);
semilogx(sigma_cph, m_ratio_2d, 'g-', 'LineWidth', 2);
hold on;
plot([sigma_c_cph, sigma_c_cph], ylim, 'r--', 'LineWidth', 1.5);
plot([f_cph, f_cph], ylim, 'k:', 'LineWidth', 1);
plot([N_cph, N_cph], ylim, 'k:', 'LineWidth', 1);
xlabel('Frequency (cph)', 'FontSize', 12);
ylabel('m_r / m_i', 'FontSize', 12);
title('Wavenumber Amplification vs Frequency for Single Wave Reflection', ...
    'FontSize', 12);
legend('Trend', '\sigma_c', 'Location', 'best');
grid on;
xlim([0.9 * f_cph, 1.1 * N_cph]); % no trend below f or above N

% Figure 3: Critical Frequency vs latitude for varying slopes

% Latitude range
latitudes = linspace(0, 90, 200);

% Slope range
slopes = linspace(0.01, 0.3, 30);
cmap = jet(length(slopes));

figure('Position', [100, 100, 700, 500]);
hold on;

% Plot sigma_c vs latitude for each slope
for i = 1:length(slopes)
    s_val = slopes(i);
    sigma_c_lat = zeros(size(latitudes));
    
    for j = 1:length(latitudes)
        lat = latitudes(j);
        f_lat = 2 * w * sin(deg2rad(lat));
        sigma_c = sqrt((N^2 * s_val^2 + f_lat^2) / (1 + s_val^2));
        sigma_c_lat(j) = sigma_c * 3600 / (2 * pi);
    end
    
    plot(latitudes, sigma_c_lat, 'Color', cmap(i,:), 'LineWidth', 1.5);
end

% Add in M2 reference
M2_cph = 1.93 / 24; % M2 tidal frequency (cph)
M2_plot = plot(latitudes, M2_cph * ones(size(latitudes)), 'k--', ...
    'LineWidth', 2);

xlabel('Latitude (degrees)', 'FontSize', 12);
ylabel('\sigma_c (cph)', 'FontSize', 12);
title('Critical Frequency vs Latitude With Varying Bottom Slopes', ...
    'FontSize', 12);

% Add colorbar to show slope values
colormap(cmap);
c = colorbar;
c.Label.String = 'Bottom slope s';
c.Ticks = linspace(0, 1, 5);
c.TickLabels = cellstr(num2str(linspace(0.01, 0.3, 5)', '%.2f'));

legend(M2_plot, 'M2 (1.93 cpd)', 'Location', 'best', 'FontSize', 9);
grid on;

%% Local helper functions

function [theta, isvalid] = compute_theta(sigma, f, N)
    % Compute incident theta_i from wave frequency from the 
    % dispersion relation
    if sigma <= f || sigma >= N
        isvalid = false;
        theta = NaN;
        return;
    end
    tan_2_theta = (N^2 - sigma^2) / (sigma^2 - f^2);
    if tan_2_theta < 0
        isvalid = false;
        theta = NaN;
        return;
    end
    isvalid = true;
    theta = atan(sqrt(tan_2_theta));
end

function m_ratio = wavenumber(theta_i, alpha)
    % Compute wavenumber amplification for 2D case
    m_ratio = -cos(alpha - theta_i) / cos(alpha + theta_i);
end
