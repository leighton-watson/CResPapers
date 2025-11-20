% Grid Search Inversion with Analytical Pipe Model
% Uses analytical formula: f = n*c/(4*Leff)
% where Leff = L + (8*a)/(3*pi) with L = depth and a = radius

clear all; clc;

%% ========================================
%% USER-DEFINED PARAMETERS
%% ========================================

% Data processing parameters
gap_time = datetime(2022, 5, 20);  % Date to split data into two sections
downsample_hours = 12;             % Time interval for downsampling (hours)

% Pipe model parameters
c_sound = 387;  % Speed of sound [m/s]

% Geometry search parameters (single radius for cylindrical pipe)
radius_range = 2:10:80;           % Pipe radius [m]

% Depth search parameters
depth_range = 20:10:200;          % Depth to search [m]

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
fprintf('========================================\n');
fprintf('BUILDING FORWARD MODEL LOOKUP TABLE\n');
fprintf('========================================\n');
lookup_table = build_lookup_table_analytical(radius_range, depth_range, c_sound);

%% INVERT SECTION 1
fprintf('\n========================================\n');
fprintf('INVERTING SECTION 1 (PRE-GAP)\n');
fprintf('========================================\n');
[geom1, depths1, misfit1, results1] = invert_pipe_model_fast(section1_sub, lookup_table);

%% INVERT SECTION 2
fprintf('\n========================================\n');
fprintf('INVERTING SECTION 2 (POST-GAP)\n');
fprintf('========================================\n');
[geom2, depths2, misfit2, results2] = invert_pipe_model_fast(section2_sub, lookup_table);

%% SAVE COMBINED RESULTS
combined_results = struct();
combined_results.section1 = results1;
combined_results.section2 = results2;
combined_results.gap_time = gap_time;
combined_results.lookup_table = lookup_table;
combined_results.c_sound = c_sound;

save('inv_results_cylindrical_pipe.mat', 'combined_results');
fprintf('\nCombined results saved to inv_results_cylindrical_pipe.mat\n\n');

%% VISUALIZE BOTH SECTIONS
visualize_split_results(results1, results2, gap_time);

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

%% BUILD LOOKUP TABLE - ANALYTICAL PIPE MODEL
function lookup = build_lookup_table_analytical(radius_range, depth_range, c_sound)
    
    n_radius = length(radius_range);
    n_depth = length(depth_range);
    
    fprintf('Computing analytical pipe model for all parameter combinations...\n');
    fprintf('Radii: %d values\n', n_radius);
    fprintf('Depths: %d values\n', n_depth);
    fprintf('Total evaluations: %d\n', n_radius * n_depth);
    fprintf('----------------------------------------\n\n');
    
    % Initialize lookup table
    % Store as 3D array: (radius_idx, depth_idx, mode_number)
    lookup.frequencies = nan(n_radius, n_depth, 3);
    lookup.radius = radius_range;
    lookup.depth = depth_range;
    lookup.c_sound = c_sound;
    
    tic_total = tic;
    
    % Loop over all radii and depths
    for i_r = 1:n_radius
        for i_d = 1:n_depth
            radius = radius_range(i_r);
            depth = depth_range(i_d);
            
            % Compute frequencies for three modes using analytical formula
            peak_freqs = analytical_pipe_frequencies(radius, depth, c_sound);
            
            % Store in lookup table
            lookup.frequencies(i_r, i_d, :) = peak_freqs;
        end
    end
    
    total_time = toc(tic_total);
    fprintf('\nLookup table complete!\n');
    fprintf('Total time: %s\n', formatTime(total_time));
    fprintf('Total evaluations: %d\n', n_radius * n_depth);
    fprintf('Average time per evaluation: %.6f s\n', total_time / (n_radius * n_depth));
    fprintf('========================================\n');
end

%% ANALYTICAL PIPE MODEL - COMPUTE FREQUENCIES
function peak_freqs = analytical_pipe_frequencies(a, L, c)
    % Analytical formula for pipe resonance
    % f = n*c/(4*Leff) where Leff = L + (8*a)/(3*pi)
    %
    % Inputs:
    %   a - pipe radius [m]
    %   L - pipe depth [m]
    %   c - speed of sound [m/s]
    % Output:
    %   peak_freqs - frequencies for modes 1, 2, 3 [Hz]
    
    % Calculate effective length with end correction
    Leff = L + (8*a)/(3*pi);
    
    % Calculate frequencies for first three modes
    peak_freqs = zeros(1, 3);
    for n = 1:3
        peak_freqs(n) = n * c / (4 * Leff);
    end
end

%% FAST INVERSION USING LOOKUP TABLE
function [best_geom, best_depths, best_misfit, results] = invert_pipe_model_fast(peak_data, lookup)
    
    n_radius = length(lookup.radius);
    n_depth = length(lookup.depth);
    n_times = height(peak_data);
    
    fprintf('Using pre-computed lookup table for inversion...\n');
    fprintf('Time steps: %d\n', n_times);
    fprintf('Total geometries: %d\n', n_radius);
    fprintf('----------------------------------------\n\n');
    
    tic_inversion = tic;
    
    % Initialize results storage
    results = struct();
    results.radius = lookup.radius;
    results.depth = lookup.depth;
    results.best_depths = zeros(n_radius, n_times);
    results.best_misfit_per_radius = inf(n_radius, 1);
    
    % Loop over radii
    for i_r = 1:n_radius
        
        % Track total misfit for this radius across all time steps
        total_misfit_radius = 0;
        
        % Loop over time steps
        for t = 1:n_times
            
            % Get observed frequencies
            f_obs = [peak_data.f1(t), peak_data.f2(t), peak_data.f3(t)];
            valid_idx = ~isnan(f_obs);
            f_obs = f_obs(valid_idx);
            
            if isempty(f_obs)
                results.best_depths(i_r, t) = nan;
                continue;
            end
            
            best_misfit_t = inf;
            best_depth_t = nan;
            
            % Loop over depths - just lookup, no computation!
            for i_d = 1:n_depth
                
                % Get modeled frequencies from lookup table
                f_model = squeeze(lookup.frequencies(i_r, i_d, :))';
                f_model = f_model(valid_idx);
                
                % Compute misfit (relative error)
                mode_weights = [1.0, 1.0, 1.0];
                weights = mode_weights(valid_idx);
                rel_error = ((f_obs - f_model) ./ f_obs).^2;
                misfit = sum(weights .* rel_error) / length(f_obs);
                
                % Track best depth for this time step
                if misfit < best_misfit_t
                    best_misfit_t = misfit;
                    best_depth_t = lookup.depth(i_d);
                end
            end
            
            % Store best depth for this time step
            results.best_depths(i_r, t) = best_depth_t;
            total_misfit_radius = total_misfit_radius + best_misfit_t;
        end
        
        % Average misfit for this radius
        results.best_misfit_per_radius(i_r) = total_misfit_radius / n_times;
    end
    
    inversion_time = toc(tic_inversion);
    
    fprintf('Inversion complete!\n');
    fprintf('Inversion time: %s\n', formatTime(inversion_time));
    
    % Find best radius
    [best_misfit, i_r_best] = min(results.best_misfit_per_radius);
    
    best_radius = lookup.radius(i_r_best);
    best_depths = results.best_depths(i_r_best, :)';
    
    best_geom = struct('radius', best_radius);
    
    fprintf('\nBEST-FIT RESULTS:\n');
    fprintf('  Pipe radius: %.0f m\n', best_radius);
    fprintf('  Average misfit: %.6e\n', best_misfit);
    fprintf('  Depth range: %.1f - %.1f m\n', min(best_depths(~isnan(best_depths))), max(best_depths(~isnan(best_depths))));
    
    % Extract modeled frequencies for best radius at each time step
    fprintf('Extracting modeled frequencies for best-fit radius...\n');
    
    modeled_freqs = nan(n_times, 3);
    for t = 1:n_times
        if ~isnan(best_depths(t))
            i_d = find(lookup.depth == best_depths(t), 1);
            if ~isempty(i_d)
                modeled_freqs(t, :) = squeeze(lookup.frequencies(i_r_best, i_d, :));
            end
        end
    end
    
    % Store best-fit parameters
    results.best_radius = best_radius;
    results.best_depths_timeseries = best_depths;
    results.best_misfit = best_misfit;
    results.peak_data = peak_data;
    results.modeled_freqs = modeled_freqs;
    results.c_sound = lookup.c_sound;
end

%% VISUALIZATION FOR SPLIT SECTIONS
function visualize_split_results(results1, results2, gap_time)
    
    %% FIGURE 1: MISFIT AND GEOMETRY
    figure('Position', [100,100,1200,800]);
    
    % Section 1 misfit vs radius
    subplot(2, 2, 1);
    plot(results1.radius, results1.best_misfit_per_radius, 'b-o', 'LineWidth', 2);
    hold on;
    plot(results1.best_radius, results1.best_misfit, 'r*', 'MarkerSize', 15, 'LineWidth', 2);
    xlabel('Radius (m)');
    ylabel('Average Misfit');
    title('(a) Period 1: Misfit vs Radius');
    grid on;
    
    % Section 1 geometry
    subplot(2, 2, 2);
    depths1 = results1.best_depths_timeseries;
    rectangle('Position', [-results1.best_radius, 0, 2*results1.best_radius, 100], ...
        'EdgeColor', 'b', 'LineWidth', 2);
    hold on;
    % plot([-results1.best_radius, results1.best_radius], [min(depths1), min(depths1)], ...
    %     'r--', 'LineWidth', 1.5);
    % plot([-results1.best_radius, results1.best_radius], [max(depths1), max(depths1)], ...
    %     'b--', 'LineWidth', 1.5);
    set(gca, 'YDir', 'reverse');
    xlabel('Radius (m)');
    ylabel('Depth (m)');
    title(sprintf('(b) Period 1 Geometry\nR = %.0f m', results1.best_radius));
    grid on; box on;
    axis equal;
    ylim([0 100]);
    xlim([-results1.best_radius*1.5, results1.best_radius*1.5]);
    % legend(sprintf('Min depth: %.0fm', min(depths1)), ...
    %     sprintf('Max depth: %.0fm', max(depths1)), 'Location', 'northwest');
    
   
    % Section 2 misfit vs radius
    subplot(2, 2, 3);
    plot(results2.radius, results2.best_misfit_per_radius, 'r-o', 'LineWidth', 2);
    hold on;
    plot(results2.best_radius, results2.best_misfit, 'r*', 'MarkerSize', 15, 'LineWidth', 2);
    xlabel('Radius (m)');
    ylabel('Average Misfit');
    title('(c) Period 2: Misfit vs Radius');
    grid on;
    
    % Section 2 geometry
    subplot(2, 2, 4);
    depths2 = results2.best_depths_timeseries;
    rectangle('Position', [-results2.best_radius, 0, 2*results2.best_radius, 100], ...
        'EdgeColor', 'r', 'LineWidth', 2);
    hold on;
    % plot([-results2.best_radius, results2.best_radius], [min(depths2), min(depths2)], ...
    %     'r--', 'LineWidth', 1.5);
    % plot([-results2.best_radius, results2.best_radius], [max(depths2), max(depths2)], ...
    %     'b--', 'LineWidth', 1.5);
    set(gca, 'YDir', 'reverse');
    xlabel('Radius (m)');
    ylabel('Depth (m)');
    title(sprintf('(d) Period 2 Geometry\nR = %.0f m', results2.best_radius));
    grid on; box on;
    axis equal;
    ylim([0 100]);
    xlim([-results2.best_radius*1.5, results2.best_radius*1.5]);
    % legend(sprintf('Min depth: %.0fm', min(depths2)), ...
    %     sprintf('Max depth: %.0fm', max(depths2)), 'Location', 'northwest');
   
    %% FIGURE 2: FREQUENCIES, MISFIT, AND DEPTH TIME SERIES
    figure('Position', [100,100,1200,1100]);
    
    colors = lines(3);
    mode_weights = [1.0, 1.0, 1.0];
    
    % Get depths and times for plotting
    depths1 = results1.best_depths_timeseries;
    depths2 = results2.best_depths_timeseries;
    t1 = results1.peak_data.time;
    t2 = results2.peak_data.time;
    
    %%% ---------------------- SUBPLOT 1: FREQUENCIES ---------------------- %%%
    subplot(5, 1, [1,2]); hold on;
    
    % Plot Section 1
    for mode = 1:3
        f_obs_1 = results1.peak_data.(['f' num2str(mode)]);
        f_model_1 = results1.modeled_freqs(:, mode);
        
        valid_idx_1 = ~isnan(f_obs_1);
        t_valid_1 = t1(valid_idx_1);
        f_obs_valid_1 = f_obs_1(valid_idx_1);
        f_model_valid_1 = f_model_1(valid_idx_1);
        
        % Observed
        scatter(t_valid_1, f_obs_valid_1, 200, colors(mode,:), ...
            'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.25, ...
            'DisplayName', sprintf('f%d Observed', mode), ...
            'MarkerEdgeColor', 'k');
        
        % Modeled
        scatter(t_valid_1, f_model_valid_1, 100, colors(mode,:), ...
            'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 0.75, ...
            'DisplayName', sprintf('f%d Modeled', mode), ...
            'MarkerEdgeColor', 'k', 'Marker', 's');
    end
    
    % Plot Section 2
    for mode = 1:3
        f_obs_2 = results2.peak_data.(['f' num2str(mode)]);
        f_model_2 = results2.modeled_freqs(:, mode);
        
        valid_idx_2 = ~isnan(f_obs_2);
        t_valid_2 = t2(valid_idx_2);
        f_obs_valid_2 = f_obs_2(valid_idx_2);
        f_model_valid_2 = f_model_2(valid_idx_2);
        
        % Observed
        scatter(t_valid_2, f_obs_valid_2, 200, colors(mode,:), ...
            'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.25, ...
            'HandleVisibility', 'off', ...
            'MarkerEdgeColor', 'k');
        
        % Modeled
        scatter(t_valid_2, f_model_valid_2, 100, colors(mode,:), ...
            'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 0.75, ...
            'HandleVisibility', 'off', ...
            'MarkerEdgeColor', 'k', 'Marker', 's');
    end
    
    ylabel('Frequency (Hz)');
    title('(a) Observed and Modelled Frequencies');
    legend('Location', 'NorthWest');
    grid on; box on;
    
    %%% ---------------------- SUBPLOT 2: MISFIT ---------------------- %%%
    subplot(5, 1, 3);
    
    % Compute combined misfit across modes - Section 1
    combined_misfit1 = nan(size(t1));
    for i = 1:length(t1)
        temp = [];
        for mode = 1:3
            f_obs = results1.peak_data.(['f' num2str(mode)]);
            f_model = results1.modeled_freqs(:, mode);
            if ~isnan(f_obs(i)) && ~isnan(f_model(i))
                temp(end+1) = mode_weights(mode) * ((f_obs(i)-f_model(i))/f_obs(i))^2;
            end
        end
        if ~isempty(temp)
            combined_misfit1(i) = mean(temp);
        end
    end
    
    % Compute combined misfit across modes - Section 2
    combined_misfit2 = nan(size(t2));
    for i = 1:length(t2)
        temp = [];
        for mode = 1:3
            f_obs = results2.peak_data.(['f' num2str(mode)]);
            f_model = results2.modeled_freqs(:, mode);
            if ~isnan(f_obs(i)) && ~isnan(f_model(i))
                temp(end+1) = mode_weights(mode) * ((f_obs(i)-f_model(i))/f_obs(i))^2;
            end
        end
        if ~isempty(temp)
            combined_misfit2(i) = mean(temp);
        end
    end
    
    % Plot misfit
    scatter(t1, combined_misfit1, 200, 'k', 'filled', ...
        'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 1, ...
        'HandleVisibility', 'off', 'MarkerEdgeColor', 'k');
    hold on;
    scatter(t2, combined_misfit2, 200, 'k', 'filled', ...
        'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 1, ...
        'MarkerEdgeColor', 'k', ...
        'DisplayName', '$\sum \left(\frac{f_{\mathrm{obs}} - f_{\mathrm{model}}}{f_{\mathrm{obs}}}\right)^2$');
    
    legend('Interpreter', 'latex', 'Location', 'Best');
    box on; grid on;
    ylim([0, 0.04]);
    set(gca, 'YScale', 'linear');
    ylabel('Misfit');
    title('(b) Misfit');
    

    disp('misfit1')
    mean(combined_misfit1)
    disp('misfit2')
    mean(combined_misfit2)

    
    %%% ---------------------- SUBPLOT 3: DEPTH ---------------------- %%%
    subplot(5, 1, [4,5]);
    scatter(t1, depths1, 200, colors(1,:), 'filled', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 1, ...
        'MarkerEdgeColor', 'k');
    hold on;
    scatter(t2, depths2, 200, colors(1,:), 'filled', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 1, ...
        'MarkerEdgeColor', 'k');
    
    xlabel('Time');
    ylabel('Magma Level Depth (m)');
    title('(c) Inferred Magma Depth Over Time');
    grid on; box on;
    set(gca, 'YDir', 'reverse');
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