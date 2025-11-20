%% plot numerical conical frustum %%

clear all; clc;
set(0,'DefaultAxesFontSize',14)

load('inv_results_analytical_conical_frustum.mat');
results1 = combined_results.section1;
results2 = combined_results.section2;
gap_time = combined_results.gap_time;

%load crater_profiles.mat

% colorbar limits of misfit colorplots
c_min = 0;
c_max = 0.05;

%% figure 1 geometry misfit %%

% Find best geometry depth indices
[~, i_gd_best1] = min(results1.best_misfit_per_geom(:));
[~, ~, gd_idx1] = ind2sub(size(results1.best_misfit_per_geom), i_gd_best1);
[~, i_gd_best2] = min(results2.best_misfit_per_geom(:));
[~, ~, gd_idx2] = ind2sub(size(results2.best_misfit_per_geom), i_gd_best2);

figHand1 = figure(1); clf;
set(figHand1,'Position', [100,100,1200,1100]); 

% Period 1 misfit heatmap (at best geometry depth)
subplot(2, 2, 1);
results1.best_misfit_per_geom(isinf(results1.best_misfit_per_geom)) = NaN;
imagesc(results1.r_inlet, results1.r_outlet, results1.best_misfit_per_geom(:, :, gd_idx1),...
        'AlphaData', ~isnan(results1.best_misfit_per_geom(:, :, gd_idx1)));
colorbar;
clim([c_min c_max]);
xlabel('Inlet Radius (m)');
ylabel('Outlet Radius (m)');
title(sprintf('(a) Period 1: Misfit\n Crater Depth = %.0f m', results1.best_geom_depth));
set(gca, 'YDir', 'normal');
hold on;
plot(results1.best_r_inlet, results1.best_r_outlet, 'r*', 'MarkerSize', 15, 'LineWidth', 2);

% find range of best fitting misfit
misfit1 = results1.best_misfit_per_geom(:, :, gd_idx1);   % 2-D slice
r_in = results1.r_inlet;
r_out = results1.r_outlet;
misfit1(isinf(misfit1)) = NaN;
mmin = min(misfit1(:), [], 'omitnan');
tol = 0.10;                      % percent
threshold = mmin * (1 + tol);    % misfit must be < this
accept = misfit1 <= threshold;
[r_in_grid, r_out_grid] = meshgrid(r_in, r_out);
r_in_accepted  = r_in_grid(accept);
r_out_accepted = r_out_grid(accept);
r_in_range  = [min(r_in_accepted),  max(r_in_accepted)];
r_out_range = [min(r_out_accepted), max(r_out_accepted)];
fprintf('Period 1 acceptable inlet radii:  %.2f – %.2f m\n', r_in_range);
fprintf('Period 1 acceptable outlet radii: %.2f – %.2f m\n', r_out_range);

% Period 2 misfit heatmap (at best geometry depth)
subplot(2, 2, 2);
results2.best_misfit_per_geom(isinf(results2.best_misfit_per_geom)) = NaN;
imagesc(results2.r_inlet, results2.r_outlet, results2.best_misfit_per_geom(:, :, gd_idx2),...
        'AlphaData', ~isnan(results2.best_misfit_per_geom(:, :, gd_idx2)));
colorbar;
clim([c_min c_max]);
xlabel('Inlet Radius (m)');
ylabel('Outlet Radius (m)');
title(sprintf('(b) Period 2: Misfit\nCrater Depth = %.0f m', results2.best_geom_depth));
set(gca, 'YDir', 'normal');
hold on;
plot(results2.best_r_inlet, results2.best_r_outlet, 'r*', 'MarkerSize', 15, 'LineWidth', 2);

% find range of best fitting misfit
misfit2 = results2.best_misfit_per_geom(:, :, gd_idx2);   % 2-D slice
r_in = results2.r_inlet;
r_out = results2.r_outlet;
misfit2(isinf(misfit2)) = NaN;
mmin = min(misfit2(:), [], 'omitnan');
threshold = mmin * (1 + tol);    % misfit must be < this
accept = misfit2 <= threshold;
[r_in_grid, r_out_grid] = meshgrid(r_in, r_out);
r_in_accepted  = r_in_grid(accept);
r_out_accepted = r_out_grid(accept);
r_in_range  = [min(r_in_accepted),  max(r_in_accepted)];
r_out_range = [min(r_out_accepted), max(r_out_accepted)];
fprintf('Period 2 acceptable inlet radii:  %.2f – %.2f m\n', r_in_range);
fprintf('Period 2 acceptable outlet radii: %.2f – %.2f m\n', r_out_range);

% period 1
subplot(2,2,[3 4]);
z_vals = (0:results1.best_geom_depth)';
r_vals1 = results1.best_r_outlet + (results1.best_r_inlet - results1.best_r_outlet) * (z_vals / results1.best_geom_depth);
plot(r_vals1, z_vals , 'b-', 'LineWidth', 2); hold on;
plot(-r_vals1, z_vals, 'b-', 'LineWidth', 2, 'HandleVisibility','off');
depths1 = results1.best_magma_depths_timeseries;
%plot([min(r_vals1) max(r_vals1)], [min(depths1) min(depths1)], 'r--', 'LineWidth', 1.5);
%plot([min(r_vals1) max(r_vals1)], [max(depths1) max(depths1)], 'b--', 'LineWidth', 1.5);
xlabel('Radius (m)');
ylabel('Depth (m)');
set(gca,'YDir','Reverse')
title('(c) Crater Geometry')
axis equal
grid on

% magma depth for period 1
r_shade = [-70 -40];
fill([r_shade(1) r_shade(2) r_shade(2) r_shade(1)], ...
     [min(depths1) min(depths1) max(depths1) max(depths1)], ...
     [0.7 0.85 1], ...          
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.5,...
     'HandleVisibility','off');        % transparency
text(mean(r_shade), mean([min(depths1) max(depths1)]), ...
     sprintf('Period 1\nMagma Range'), ...   % <-- two lines
     'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'middle', ...
     'Rotation', 90, ...
     'FontWeight', 'bold', ...
     'Color', 'b', 'FontSize', 11);

% period 2
subplot(2,2,[3 4]);
z_vals = (0:results2.best_geom_depth)';
r_vals2 = results2.best_r_outlet + (results2.best_r_inlet - results2.best_r_outlet) * (z_vals / results2.best_geom_depth);
plot(r_vals2, z_vals , 'r-', 'LineWidth', 2); hold on;
plot(-r_vals2, z_vals, 'r-', 'LineWidth', 2, 'HandleVisibility','off');
depths2 = results2.best_magma_depths_timeseries;
xlabel('Radius (m)');
ylabel('Depth (m)');
legend('Period 1','Period 2','Location','SouthEast')

% magma depth for period 2
r_shade = [40 70];
fill([r_shade(1) r_shade(2) r_shade(2) r_shade(1)], ...
     [min(depths2) min(depths2) max(depths2) max(depths2)], ...
     [1 0.8 0.8], ...         
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.5,...
     'HandleVisibility','off');        % transparency
text(mean(r_shade), mean([min(depths2) max(depths2)]),...
      sprintf('Period 2\nMagma Range'),...
     'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'middle', ...
     'Rotation', 90, ...       % vertical text
     'FontWeight', 'bold', ...
     'Color', 'r','FontSize',11);

% print out crater dimensions
text(-150, 25, 'Period 1:', 'FontSize', 16, 'FontWeight', 'bold');
text(-150, 35, sprintf('Outlet radius: %.0f m', results1.best_r_outlet), 'FontSize', 14);
text(-150, 45, sprintf('Inlet radius: %.0f m', results1.best_r_inlet), 'FontSize', 14);
text(-150, 55, sprintf('Crater depth: %.0f m', results1.best_geom_depth), 'FontSize', 14);
text(-150, 65, sprintf('Magma range: %.0f - %.0f m', min(depths1), max(depths1)), 'FontSize', 14);

text(80, 25, 'Period 2:', 'FontSize', 16, 'FontWeight', 'bold');
text(80, 35, sprintf('Outlet radius: %.0f m', results2.best_r_outlet), 'FontSize', 14);
text(80, 45, sprintf('Inlet radius: %.0f m', results2.best_r_inlet), 'FontSize', 14);
text(80, 55, sprintf('Crater depth: %.0f m', results2.best_geom_depth), 'FontSize', 14);
text(80, 65, sprintf('Magma range: %.0f - %.0f m', min(depths2), max(depths2)), 'FontSize', 14);



%% figure 2 time series peak frequency and depth and misfit

figure('Position', [100,100,1200,1000]); clf;
colors = lines(3);

depths1 = results1.best_magma_depths_timeseries;
depths2 = results2.best_magma_depths_timeseries;
mode_weights = [1.0, 1.0, 1.0];
t1 = results1.peak_data.time;
t2 = results2.peak_data.time;


%%% ---------------------- RESULTS 1 ---------------------- %%%
for mode = 1:3
    % Observed & modeled frequencies (Section 1)
    f_obs_1 = results1.peak_data.(['f' num2str(mode)]);
    f_model_1 = results1.modeled_freqs(:, mode);

    % Mask valid (non-NaN) observed points
    valid_idx_1 = ~isnan(f_obs_1);
    t_valid_1 = t1(valid_idx_1);
    f_obs_valid_1 = f_obs_1(valid_idx_1);
    f_model_valid_1 = f_model_1(valid_idx_1);

    % --- Compute time-varying misfit ---
    w1 = mode_weights(mode);
    rel_error_1 = ((f_obs_valid_1 - f_model_valid_1) ./ f_obs_valid_1).^2;
    misfit_t_1 = w1 .* rel_error_1; % misfit per time step
    
    % --- Plot observed and modeled frequencies ---
    subplot(2,1,1); hold on
    scatter(t_valid_1, f_obs_valid_1, 200, colors(mode,:), ...
        'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none', ...
        'DisplayName', sprintf('f_%d Observed', mode));
    scatter(t_valid_1, f_model_valid_1, 120, colors(mode,:), ...
        'filled', 'MarkerFaceAlpha', 1, ...
        'DisplayName', sprintf('f_%d Modeled', mode),...
        'MarkerEdgeColor','k','Marker','s');
    
    ylabel('Frequency (Hz)');
    title('Observed vs Modeled Frequencies');
    legend('Location','NorthWest');
    grid on; box on;

end

%%% ---------------------- RESULTS 2 ---------------------- %%%
for mode = 1:3
    % Observed & modeled frequencies (Section 2)
    f_obs_2 = results2.peak_data.(['f' num2str(mode)]);
    f_model_2 = results2.modeled_freqs(:, mode);

    % Mask valid (non-NaN) observed points
    valid_idx_2 = ~isnan(f_obs_2);
    t_valid_2 = t2(valid_idx_2);
    f_obs_valid_2 = f_obs_2(valid_idx_2);
    f_model_valid_2 = f_model_2(valid_idx_2);

    % --- Compute time-varying misfit ---
    w2 = mode_weights(mode);
    rel_error_2 = ((f_obs_valid_2 - f_model_valid_2) ./ f_obs_valid_2).^2;
    misfit_t_2 = w2 .* rel_error_2; % misfit per time step

    % --- Plot observed and modeled frequencies ---
    subplot(2,1,1); hold on
    scatter(t_valid_2, f_obs_valid_2, 200, colors(mode,:), ...
        'filled', 'MarkerFaceAlpha', 0.25, ...
        'HandleVisibility', 'off', ...
        'MarkerEdgeColor','none');
    scatter(t_valid_2, f_model_valid_2, 100, colors(mode,:), ...
        'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 0.75, ...
        'HandleVisibility', 'off',...
        'MarkerEdgeColor','k','Marker','s');
end
title('(a) Observed and Modelled Frequencies')


%%% DEPTH %%%

subplot(2, 1, 2); 
scatter(results1.peak_data.time, depths1, 120, colors(1,:), 's', ...
    'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 1,...
    'MarkerEdgeColor','k');
hold on;
scatter(results2.peak_data.time, depths2, 120, colors(1,:), 's', 'filled', ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 1,...
    'MarkerEdgeColor','k');
ylabel('Magma Level Depth (m)');
title('(b) Inferred Magma Depth Over Time');
grid on; box on;
set(gca, 'YDir', 'reverse');
ylim([80 160])


t1_start = results1.peak_data.time(1);
t1_end = results1.peak_data.time(end);
t2_start = results2.peak_data.time(1);
t2_end = results2.peak_data.time(end);

fill([t1_start t1_end t1_end t1_start], ...
    [80 80 90 90],...
    [0.7 0.85 1], ...          % light blue color
     'FaceAlpha', 0.5, ...       % transparency
     'EdgeColor', 'none');       % no border

fill([t2_start t2_end t2_end t2_start], ...
    [80 80 90 90],...
    [1 0.8 0.8], ...          % light blue color
     'FaceAlpha', 0.5, ...       % transparency
     'EdgeColor', 'none');       % no border

text(mean([t1_start, t1_end]), 85, ...
     sprintf('Period 1'), ...   
     'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'middle', ...
     'FontWeight', 'bold', ...
     'Color', 'b', 'FontSize', 14);

text(mean([t2_start, t2_end]), 85, ...
     sprintf('Period 2'), ...   
     'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'middle', ...
     'FontWeight', 'bold', ...
     'Color', 'r', 'FontSize', 14);