% Grid Search Inversion with Conical Frustum Model
% Uses analytical transcendental equation with lookup table

clear all; clc;

%% ========================================
%% USER-DEFINED PARAMETERS
%% ========================================

% Data processing parameters
gap_time = datetime(2022, 5, 20);  % Date to split data into two sections
downsample_hours = 24;             % Time interval for downsampling (hours)

% Frustum model parameters
c_sound = 387;  % Speed of sound [m/s]

% Geometry search parameters (linear taper: alpha = 1)
r_outlet_range = 8:10:80;         % Outlet radius (top) [m]
r_inlet_range = 2:10:42;          % Inlet radius (bottom) [m]
geom_depth_range = 60:10:160;     % Maximum depth of geometry (where inlet is defined) [m]

% Magma level depth search parameters
magma_depth_range = 60:10:160;    % Depth of magma surface to search [m]

%% ========================================
%% LOAD DATA
%% ========================================
peak_data = readtable('averaged_peak_frequencies_manual_edit.csv');

%% SPLIT DATA INTO SECTIONS
section1_idx = peak_data.time < gap_time;
section2_idx = peak_data.time >= gap_time;

section1_data = peak_data(section1_idx, :);
section2_data = peak_data(section2_idx, :);

fprintf('\n========================================\n');
fprintf('DATA SECTIONS\n');
fprintf('========================================\n');
fprintf('Section 1 (pre-gap): %d observations\n', height(section1_data));
fprintf('Section 2 (post-gap): %d observations\n', height(section2_data));
fprintf('========================================\n\n');

%% SUBSAMPLE DATA
section1_sub = subsample_hr(section1_data, downsample_hours);
section2_sub = subsample_hr(section2_data, downsample_hours);

fprintf('After %d-hour subsampling:\n', downsample_hours);
fprintf('Section 1: %d time windows\n', height(section1_sub));
fprintf('Section 2: %d time windows\n', height(section2_sub));
fprintf('Total time windows: %d\n\n', height(section1_sub) + height(section2_sub));

%% PRE-COMPUTE FORWARD MODEL LOOKUP TABLE
fprintf('\n========================================\n');
fprintf('BUILDING FORWARD MODEL LOOKUP TABLE\n');
fprintf('========================================\n');
lookup_table = build_lookup_table_frustum(r_outlet_range, r_inlet_range, ...
    geom_depth_range, magma_depth_range, c_sound);

%% INVERT SECTION 1
fprintf('\n========================================\n');
fprintf('INVERTING SECTION 1 (PRE-GAP)\n');
fprintf('========================================\n');
[geom1, depths1, misfit1, results1] = invert_geometry_depth_fast(section1_sub, lookup_table);

%% INVERT SECTION 2
fprintf('\n========================================\n');
fprintf('INVERTING SECTION 2 (POST-GAP)\n');
fprintf('========================================\n');
[geom2, depths2, misfit2, results2] = invert_geometry_depth_fast(section2_sub, lookup_table);

%% SAVE COMBINED RESULTS
combined_results = struct();
combined_results.section1 = results1;
combined_results.section2 = results2;
combined_results.gap_time = gap_time;
combined_results.lookup_table = lookup_table;

save('inv_results_analytical_conical_frustum.mat', 'combined_results');
fprintf('\nCombined results saved to inv_results_analytical_conical_frustum.mat\n\n');

%% VISUALIZE BOTH SECTIONS
visualize_split_results(results1, results2, gap_time);

%% ========================================
%% FUNCTIONS
%% ========================================

%% SUBSAMPLE FUNCTION
function subsampled_data = subsample_hr(data, downsample_hours)
    if isempty(data)
        subsampled_data = data;
        return;
    end
    
    % Sort by time
    data = sortrows(data, 'time');
    
    % Initialize with first observation
    subsampled_data = data(1, :);
    last_time = data.time(1);
    
    % Loop through and select observations at least N hours apart
    for i = 2:height(data)
        time_diff = hours(data.time(i) - last_time);
        if time_diff >= downsample_hours
            subsampled_data = [subsampled_data; data(i, :)];
            last_time = data.time(i);
        end
    end
end

%% BUILD LOOKUP TABLE - CONICAL FRUSTUM MODEL
function lookup = build_lookup_table_frustum(r_outlet_range, r_inlet_range, ...
    geom_depth_range, magma_depth_range, c_sound)
    
    n_outlet = length(r_outlet_range);
    n_inlet = length(r_inlet_range);
    n_geom_depth = length(geom_depth_range);
    n_magma_depth = length(magma_depth_range);
    
    fprintf('Computing conical frustum model for all parameter combinations...\n');
    fprintf('Outlet radii: %d values\n', n_outlet);
    fprintf('Inlet radii: %d values\n', n_inlet);
    fprintf('Geometry depths: %d values\n', n_geom_depth);
    fprintf('Magma level depths: %d values\n', n_magma_depth);
    fprintf('Total evaluations: %d\n', n_outlet * n_inlet * n_geom_depth * n_magma_depth);
    fprintf('----------------------------------------\n\n');
    
    % Initialize lookup table
    % Store as 5D array: (outlet_idx, inlet_idx, geom_depth_idx, magma_depth_idx, mode_number)
    lookup.frequencies = nan(n_outlet, n_inlet, n_geom_depth, n_magma_depth, 3);
    lookup.acoustic_inlet_radii = nan(n_outlet, n_inlet, n_geom_depth, n_magma_depth);
    lookup.r_outlet = r_outlet_range;
    lookup.r_inlet = r_inlet_range;
    lookup.geom_depth = geom_depth_range;
    lookup.magma_depth = magma_depth_range;
    lookup.valid_geom = false(n_outlet, n_inlet, n_geom_depth);
    lookup.c_sound = c_sound;
    
    tic_total = tic;
    
    % Progress tracking
    total_combinations = n_outlet * n_inlet * n_geom_depth * n_magma_depth;
    completed = 0;
    last_percent = 0;
    
    fprintf('Progress: [');
    
    % Loop over all geometries
    for i_out = 1:n_outlet
        for i_in = 1:n_inlet
            for i_gd = 1:n_geom_depth
                r_outlet = r_outlet_range(i_out);
                r_inlet = r_inlet_range(i_in);
                geom_depth = geom_depth_range(i_gd);
                
                % Check if geometry is valid (outlet must be larger than inlet)
                if r_outlet <= r_inlet
                    completed = completed + n_magma_depth;
                    continue;
                end
                
                lookup.valid_geom(i_out, i_in, i_gd) = true;
                
                % Loop over magma depths
                for i_md = 1:n_magma_depth
                    magma_depth = magma_depth_range(i_md);
                    
                    % Skip if magma depth is deeper than geometry depth
                    if magma_depth > geom_depth
                        completed = completed + 1;
                        continue;
                    end
                    
                    % Calculate acoustic resonator inlet radius at this magma depth
                    % The crater geometry goes from r_outlet (at depth 0) to r_inlet (at geom_depth)
                    % At magma_depth, the radius is interpolated:
                    acoustic_inlet = r_outlet - (r_outlet - r_inlet) * (magma_depth / geom_depth);
                    acoustic_outlet = r_outlet;  % Always at the top
                    acoustic_length = magma_depth;  % From top to magma surface
                    
                    % Store acoustic inlet radius
                    lookup.acoustic_inlet_radii(i_out, i_in, i_gd, i_md) = acoustic_inlet;
                    
                    % Compute frequencies using analytical frustum formula
                    peak_freqs = frustum_frequencies(acoustic_inlet, acoustic_outlet, ...
                        acoustic_length, c_sound);
                    
                    % Store in lookup table
                    lookup.frequencies(i_out, i_in, i_gd, i_md, :) = peak_freqs;
                    
                    % Update progress
                    completed = completed + 1;
                    current_percent = floor(100 * completed / total_combinations);
                    if current_percent > last_percent
                        fprintf('%d%% ', current_percent);
                        last_percent = current_percent;
                    end
                end
            end
        end
    end
    fprintf('] Complete!\n');
    
    total_time = toc(tic_total);
    fprintf('\nLookup table complete!\n');
    fprintf('Total time: %s\n', formatTime(total_time));
    fprintf('Total evaluations: %d\n', total_combinations);
    fprintf('Average time per evaluation: %.6f s\n', total_time / total_combinations);
    fprintf('========================================\n');
end

%% CONICAL FRUSTUM MODEL - COMPUTE FREQUENCIES
function peak_freqs = frustum_frequencies(a, b, L, c)
    % Solve transcendental equation for conical frustum resonance
    % tan(pi*f/f0) = -B/(1-B*pi*f/f0)
    % where f0 = c/(2*L) and B = a/b
    %
    % Inputs:
    %   a - inlet radius at closed end (magma surface) [m]
    %   b - outlet radius at open end (top) [m]
    %   L - length of resonator (magma depth) [m]
    %   c - speed of sound [m/s]
    % Output:
    %   peak_freqs - frequencies for modes 1, 2, 3 [Hz]
    
    f0 = c / (2 * L);
    B = a / b;
    
    peak_freqs = zeros(1, 3);
    
    % For each mode, find the root of the transcendental equation
    for n = 1:3
        % Initial guess based on mode number
        f_guess = n * f0;
        
        % Define the transcendental equation
        transcendental_eq = @(f) tan(pi*f/f0) + B/(1 - B*pi*f/f0);
        
        % Search range around the guess
        f_min = max(0.5 * f_guess, 0.01);
        f_max = 1.5 * f_guess;
        
        try
            % Use fzero to find the root
            peak_freqs(n) = fzero(transcendental_eq, [f_min, f_max]);
        catch
            % If fzero fails, try alternative starting points
            try
                peak_freqs(n) = fzero(transcendental_eq, f_guess);
            catch
                % If all else fails, return NaN
                peak_freqs(n) = NaN;
            end
        end
    end
end

%% FAST INVERSION USING LOOKUP TABLE
function [best_geom, best_depths, best_misfit, results] = invert_geometry_depth_fast(...
    peak_data, lookup)
    
    n_outlet = length(lookup.r_outlet);
    n_inlet = length(lookup.r_inlet);
    n_geom_depth = length(lookup.geom_depth);
    n_magma_depth = length(lookup.magma_depth);
    n_times = height(peak_data);
    
    fprintf('Using pre-computed lookup table for inversion...\n');
    fprintf('Time steps: %d\n', n_times);
    fprintf('Valid geometries: %d\n', sum(lookup.valid_geom(:)));
    fprintf('----------------------------------------\n\n');
    
    tic_inversion = tic;
    
    % Initialize results storage
    results = struct();
    results.r_outlet = lookup.r_outlet;
    results.r_inlet = lookup.r_inlet;
    results.geom_depth = lookup.geom_depth;
    results.magma_depth = lookup.magma_depth;
    results.best_magma_depths = nan(n_outlet, n_inlet, n_geom_depth, n_times);
    results.best_misfit_per_geom = inf(n_outlet, n_inlet, n_geom_depth);
    
    % Progress tracking
    total_geometries = n_outlet * n_inlet * n_geom_depth;
    completed_geom = 0;
    last_percent_inv = 0;
    
    fprintf('Inversion progress: [');
    
    % Loop over geometries
    for i_out = 1:n_outlet
        for i_in = 1:n_inlet
            for i_gd = 1:n_geom_depth
                
                % Skip invalid geometries
                if ~lookup.valid_geom(i_out, i_in, i_gd)
                    completed_geom = completed_geom + 1;
                    continue;
                end
                
                % Track total misfit for this geometry across all time steps
                total_misfit_geom = 0;
                valid_times = 0;
                
                % Loop over time steps
                for t = 1:n_times
                    
                    % Get observed frequencies
                    f_obs = [peak_data.f1(t), peak_data.f2(t), peak_data.f3(t)];
                    valid_idx = ~isnan(f_obs);
                    f_obs = f_obs(valid_idx);
                    
                    if isempty(f_obs)
                        continue;
                    end
                    
                    best_misfit_t = inf;
                    best_magma_depth_t = nan;
                    
                    % Loop over magma depths - just lookup, no computation!
                    for i_md = 1:n_magma_depth
                        
                        % Get modeled frequencies from lookup table
                        f_model = squeeze(lookup.frequencies(i_out, i_in, i_gd, i_md, :))';
                        
                        % Skip if any frequencies are NaN (invalid or failed to solve)
                        if any(isnan(f_model))
                            continue;
                        end
                        
                        f_model = f_model(valid_idx);
                        
                        % Compute misfit (relative error)
                        mode_weights = [1.0, 1.0, 1.0];
                        weights = mode_weights(valid_idx);
                        rel_error = ((f_obs - f_model) ./ f_obs).^2;
                        misfit = sum(weights .* rel_error) / length(f_obs);
                        
                        % Track best magma depth for this time step
                        if misfit < best_misfit_t
                            best_misfit_t = misfit;
                            best_magma_depth_t = lookup.magma_depth(i_md);
                        end
                    end
                    
                    % Store best magma depth for this time step
                    if ~isnan(best_magma_depth_t)
                        results.best_magma_depths(i_out, i_in, i_gd, t) = best_magma_depth_t;
                        total_misfit_geom = total_misfit_geom + best_misfit_t;
                        valid_times = valid_times + 1;
                    end
                end
                
                % Average misfit for this geometry
                if valid_times > 0
                    results.best_misfit_per_geom(i_out, i_in, i_gd) = total_misfit_geom / valid_times;
                end
                
                % Update progress
                completed_geom = completed_geom + 1;
                current_percent_inv = floor(100 * completed_geom / total_geometries);
                if current_percent_inv > last_percent_inv
                    fprintf('%d%% ', current_percent_inv);
                    last_percent_inv = current_percent_inv;
                end
            end
        end
    end
    fprintf('] Complete!\n');
    
    inversion_time = toc(tic_inversion);
    
    fprintf('Inversion complete!\n');
    fprintf('Inversion time: %s\n', formatTime(inversion_time));
    
    % Find best geometry
    [best_misfit, min_idx] = min(results.best_misfit_per_geom(:));
    [i_out_best, i_in_best, i_gd_best] = ind2sub(size(results.best_misfit_per_geom), min_idx);
    
    best_r_outlet = lookup.r_outlet(i_out_best);
    best_r_inlet = lookup.r_inlet(i_in_best);
    best_geom_depth = lookup.geom_depth(i_gd_best);
    best_magma_depths = squeeze(results.best_magma_depths(i_out_best, i_in_best, i_gd_best, :));
    
    best_geom = struct('r_outlet', best_r_outlet, 'r_inlet', best_r_inlet, 'geom_depth', best_geom_depth);
    best_depths = best_magma_depths;
    
    fprintf('\nBEST-FIT RESULTS:\n');
    fprintf('  Outlet radius: %.0f m\n', best_r_outlet);
    fprintf('  Inlet radius: %.0f m\n', best_r_inlet);
    fprintf('  Geometry depth: %.0f m\n', best_geom_depth);
    fprintf('  Average misfit: %.6e\n', best_misfit);
    fprintf('  Magma depth range: %.1f - %.1f m\n', ...
        min(best_magma_depths(~isnan(best_magma_depths))), ...
        max(best_magma_depths(~isnan(best_magma_depths))));
    
    % Extract modeled frequencies and acoustic inlet radii for best geometry at each time step
    fprintf('Extracting modeled frequencies for best-fit geometry...\n');
    
    modeled_freqs = nan(n_times, 3);
    acoustic_inlet_radii = nan(n_times, 1);
    
    for t = 1:n_times
        if ~isnan(best_magma_depths(t))
            % Find magma depth index
            i_md = find(lookup.magma_depth == best_magma_depths(t), 1);
            if ~isempty(i_md)
                % Get frequencies from lookup table
                modeled_freqs(t, :) = squeeze(lookup.frequencies(i_out_best, i_in_best, i_gd_best, i_md, :));
                % Get acoustic inlet radius from lookup table
                acoustic_inlet_radii(t) = lookup.acoustic_inlet_radii(i_out_best, i_in_best, i_gd_best, i_md);
            end
        end
    end
    
    % Store best-fit parameters
    results.best_r_outlet = best_r_outlet;
    results.best_r_inlet = best_r_inlet;
    results.best_geom_depth = best_geom_depth;
    results.best_magma_depths_timeseries = best_magma_depths;
    results.acoustic_inlet_radii = acoustic_inlet_radii;
    results.best_misfit = best_misfit;
    results.peak_data = peak_data;
    results.modeled_freqs = modeled_freqs;
    results.c_sound = lookup.c_sound;
end

%% VISUALIZATION FOR SPLIT SECTIONS
function visualize_split_results(results1, results2, gap_time)
    
    %% FIGURE 1: MISFIT AND GEOMETRY
    figure('Position', [100, 100, 1400, 900]);
    
    % Find best geometry depth indices
    [~, i_gd_best1] = min(results1.best_misfit_per_geom(:));
    [~, ~, gd_idx1] = ind2sub(size(results1.best_misfit_per_geom), i_gd_best1);
    [~, i_gd_best2] = min(results2.best_misfit_per_geom(:));
    [~, ~, gd_idx2] = ind2sub(size(results2.best_misfit_per_geom), i_gd_best2);
    
    % Section 1 misfit heatmap (at best geometry depth)
    subplot(3, 3, 1);
    imagesc(results1.r_inlet, results1.r_outlet, results1.best_misfit_per_geom(:, :, gd_idx1));
    colorbar;
    xlabel('Inlet Radius (m)');
    ylabel('Outlet Radius (m)');
    title(sprintf('Section 1: Misfit (Geom Depth = %.0fm)', results1.best_geom_depth));
    set(gca, 'YDir', 'normal');
    hold on;
    plot(results1.best_r_inlet, results1.best_r_outlet, 'r*', 'MarkerSize', 15, 'LineWidth', 2);
    
    % Section 1 geometry
    subplot(3, 3, 2);
    z_vals = (0:results1.best_geom_depth)';
    r_vals1 = results1.best_r_outlet - (results1.best_r_outlet - results1.best_r_inlet) * (z_vals / results1.best_geom_depth);
    plot(r_vals1, z_vals, 'b-', 'LineWidth', 2); hold on;
    plot(-r_vals1, z_vals, 'b-', 'LineWidth', 2);
    depths1 = results1.best_magma_depths_timeseries;
    valid_depths1 = depths1(~isnan(depths1));
    plot([min(r_vals1) max(r_vals1)], [min(valid_depths1) min(valid_depths1)], 'r--', 'LineWidth', 1.5);
    plot([min(r_vals1) max(r_vals1)], [max(valid_depths1) max(valid_depths1)], 'b--', 'LineWidth', 1.5);
    set(gca, 'YDir', 'reverse');
    xlabel('Radius (m)');
    ylabel('Depth (m)');
    title(sprintf('Section 1 Geometry\nR_{out}=%.0fm, R_{in}=%.0fm, D_{geom}=%.0fm', ...
        results1.best_r_outlet, results1.best_r_inlet, results1.best_geom_depth));
    grid on;
    axis equal;
    ylim([0 results1.best_geom_depth]);
    legend('', '', sprintf('Min magma: %.0fm', min(valid_depths1)), ...
        sprintf('Max magma: %.0fm', max(valid_depths1)), 'Location', 'southeast');
    
    % Section 1 statistics
    subplot(3, 3, 3);
    axis off;
    text(0.1, 0.9, 'Section 1 Statistics:', 'FontSize', 12, 'FontWeight', 'bold');
    text(0.1, 0.75, sprintf('Outlet radius: %.0f m', results1.best_r_outlet), 'FontSize', 10);
    text(0.1, 0.65, sprintf('Inlet radius: %.0f m', results1.best_r_inlet), 'FontSize', 10);
    text(0.1, 0.55, sprintf('Geometry depth: %.0f m', results1.best_geom_depth), 'FontSize', 10);
    text(0.1, 0.45, sprintf('Average misfit: %.2e', results1.best_misfit), 'FontSize', 10);
    text(0.1, 0.35, sprintf('Magma range: %.0f - %.0f m', min(valid_depths1), max(valid_depths1)), 'FontSize', 10);
    
    % Section 2 misfit heatmap (at best geometry depth)
    subplot(3, 3, 4);
    imagesc(results2.r_inlet, results2.r_outlet, results2.best_misfit_per_geom(:, :, gd_idx2));
    colorbar;
    xlabel('Inlet Radius (m)');
    ylabel('Outlet Radius (m)');
    title(sprintf('Section 2: Misfit (Geom Depth = %.0fm)', results2.best_geom_depth));
    set(gca, 'YDir', 'normal');
    hold on;
    plot(results2.best_r_inlet, results2.best_r_outlet, 'r*', 'MarkerSize', 15, 'LineWidth', 2);
    
    % Section 2 geometry
    subplot(3, 3, 5);
    z_vals2 = (0:results2.best_geom_depth)';
    r_vals2 = results2.best_r_outlet - (results2.best_r_outlet - results2.best_r_inlet) * (z_vals2 / results2.best_geom_depth);
    plot(r_vals2, z_vals2, 'r-', 'LineWidth', 2); hold on;
    plot(-r_vals2, z_vals2, 'r-', 'LineWidth', 2);
    depths2 = results2.best_magma_depths_timeseries;
    valid_depths2 = depths2(~isnan(depths2));
    plot([min(r_vals2) max(r_vals2)], [min(valid_depths2) min(valid_depths2)], 'r--', 'LineWidth', 1.5);
    plot([min(r_vals2) max(r_vals2)], [max(valid_depths2) max(valid_depths2)], 'b--', 'LineWidth', 1.5);
    set(gca, 'YDir', 'reverse');
    xlabel('Radius (m)');
    ylabel('Depth (m)');
    title(sprintf('Section 2 Geometry\nR_{out}=%.0fm, R_{in}=%.0fm, D_{geom}=%.0fm', ...
        results2.best_r_outlet, results2.best_r_inlet, results2.best_geom_depth));
    grid on;
    axis equal;
    ylim([0 results2.best_geom_depth]);
    legend('', '', sprintf('Min magma: %.0fm', min(valid_depths2)), ...
        sprintf('Max magma: %.0fm', max(valid_depths2)), 'Location', 'southeast');
    
    % Section 2 statistics
    subplot(3, 3, 6);
    axis off;
    text(0.1, 0.9, 'Section 2 Statistics:', 'FontSize', 12, 'FontWeight', 'bold');
    text(0.1, 0.75, sprintf('Outlet radius: %.0f m', results2.best_r_outlet), 'FontSize', 10);
    text(0.1, 0.65, sprintf('Inlet radius: %.0f m', results2.best_r_inlet), 'FontSize', 10);
    text(0.1, 0.55, sprintf('Geometry depth: %.0f m', results2.best_geom_depth), 'FontSize', 10);
    text(0.1, 0.45, sprintf('Average misfit: %.2e', results2.best_misfit), 'FontSize', 10);
    text(0.1, 0.35, sprintf('Magma range: %.0f - %.0f m', min(valid_depths2), max(valid_depths2)), 'FontSize', 10);
    
    % Comparison subplot
    subplot(3, 3, 7:9);
    bar_data = [results1.best_r_outlet, results2.best_r_outlet; ...
                results1.best_r_inlet, results2.best_r_inlet; ...
                results1.best_geom_depth, results2.best_geom_depth];
    bar(bar_data);
    set(gca, 'XTickLabel', {'Outlet Radius', 'Inlet Radius', 'Geometry Depth'});
    ylabel('Value (m)');
    legend('Section 1', 'Section 2', 'Location', 'best');
    title('Geometry Comparison Between Sections');
    grid on;
    
    sgtitle('Grid Search Results: Misfit and Geometry (Analytical Frustum)', 'FontSize', 14, 'FontWeight', 'bold');
    
    %% FIGURE 2: FREQUENCIES AND DEPTH TIME SERIES
    figure('Position', [150, 150, 1400, 700]);
    
    colors = lines(3);
    
    % Top subplot: All three modes together
    subplot(2, 1, 1);
    hold on;
    
    for mode = 1:3
        % Section 1 - observed
        f_obs1 = results1.peak_data.(['f' num2str(mode)]);
        plot(results1.peak_data.time, f_obs1, 'o', 'Color', colors(mode,:), ...
            'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', sprintf('f%d Obs (Sec 1)', mode));
        
        % Section 1 - modeled
        f_model1 = results1.modeled_freqs(:, mode);
        plot(results1.peak_data.time, f_model1, '-', 'Color', colors(mode,:), ...
            'LineWidth', 2, 'DisplayName', sprintf('f%d Model (Sec 1)', mode));
        
        % Section 2 - observed
        f_obs2 = results2.peak_data.(['f' num2str(mode)]);
        plot(results2.peak_data.time, f_obs2, 's', 'Color', colors(mode,:), ...
            'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', sprintf('f%d Obs (Sec 2)', mode));
        
        % Section 2 - modeled
        f_model2 = results2.modeled_freqs(:, mode);
        plot(results2.peak_data.time, f_model2, '--', 'Color', colors(mode,:), ...
            'LineWidth', 2, 'DisplayName', sprintf('f%d Model (Sec 2)', mode));
    end
    
    % Add gap line
    xline(gap_time, 'k--', 'LineWidth', 2, 'DisplayName', 'Gap');
    
    xlabel('Time');
    ylabel('Frequency (Hz)');
    title('Peak Frequencies: Observed vs Modeled (All Modes)');
    grid on;
    
    % Bottom subplot: Depth time series
    subplot(2, 1, 2);
    plot(results1.peak_data.time, depths1, 'b-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Section 1');
    hold on;
    plot(results2.peak_data.time, depths2, 'r-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Section 2');
    xline(gap_time, 'k--', 'LineWidth', 2, 'DisplayName', 'Gap');
    xlabel('Time');
    ylabel('Magma Level Depth (m)');
    title('Inferred Magma Depth Over Time');
    legend('Location', 'best');
    grid on;
    set(gca, 'YDir', 'reverse');
    
    sgtitle('Peak Frequencies and Inferred Depth (Analytical Frustum)', 'FontSize', 14, 'FontWeight', 'bold');
end

%% Helper function for time formatting
function str = formatTime(seconds)
    if seconds < 60
        str = sprintf('%.1fs', seconds);
    elseif seconds < 3600
        mins = floor(seconds / 60);
        secs = mod(seconds, 60);
        str = sprintf('%dm %.0fs', mins, secs);
    else
        hrs = floor(seconds / 3600);
        mins = floor(mod(seconds, 3600) / 60);
        str = sprintf('%dh %dm', hrs, mins);
    end
end


