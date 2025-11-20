%% plot cylindrical pipe %%

clear all; clc;
set(0,'DefaultAxesFontSize',14)

load('inv_results_cylindrical_pipe.mat');
results1 = combined_results.section1;
results2 = combined_results.section2;

c_min = 0;
c_max = 0.05;

%% ===================== FIGURE 1 =====================
% Misfit vs radius + inferred depth profiles

figure(1); clf;
set(gcf,'Position',[100 100 1200 900]);

% ---- Radius–Misfit Curve: Period 1 ----
subplot(2,2,1)
plot(results1.radius, results1.best_misfit_per_radius, 'b-o','LineWidth',2)
hold on
xline(results1.best_radius,'k--','LineWidth',2)
xlabel('Radius (m)')
ylabel('Misfit')
title('(a) Period 1: Misfit vs Radius')
grid on

% ---- Radius–Misfit Curve: Period 2 ----
subplot(2,2,3)
plot(results2.radius, results2.best_misfit_per_radius, 'r-o','LineWidth',2)
hold on
xline(results2.best_radius,'k--','LineWidth',2)
xlabel('Radius (m)')
ylabel('Misfit')
title('(c) Period 2: Misfit vs Radius')
grid on

% ---- Find uncertainty in radius values ----
% What radii values have misfit values within 1% of the minimum values
relTol = 0.1;      % 1 percent

% period 1
r = results1.radius;
m = results1.best_misfit_per_radius;
[m_min, idx_min] = min(m);
r_best = r(idx_min);
threshold = m_min * (1 + relTol); % Compute threshold
inBand = m <= threshold; % Find r values within the threshold
r_band = r(inBand);
r_lo = min(r_band); % Extract interval
r_hi = max(r_band);
fprintf('Period 1: Best radius: %.6g\n', r_best);
fprintf('Period 1: 1%% interval: [%.6g, %.6g]\n', r_lo, r_hi);

subplot(2,2,1);
fill([r_lo r_hi r_hi r_lo], ...
    [1.5e-3 1.5e-3 threshold threshold],...
    [0.7 0.85 1], ...          % light blue color
     'FaceAlpha', 0.5, ...       % transparency
     'EdgeColor', 'none');       % no border


% period 2
r = results2.radius;
m = results2.best_misfit_per_radius;
[m_min, idx_min] = min(m);
r_best = r(idx_min);
threshold = m_min * (1 + relTol); % Compute threshold
inBand = m <= threshold; % Find r values within the threshold
r_band = r(inBand);
r_lo = min(r_band); % Extract interval
r_hi = max(r_band);
fprintf('Period 2: Best radius: %.6g\n', r_best);
fprintf('Period 2: 1%% interval: [%.6g, %.6g]\n', r_lo, r_hi);

subplot(2,2,3);
fill([r_lo r_hi r_hi r_lo], ...
    [0.02 0.02 threshold threshold],...
    [1 0.8 0.8], ...          % light red color
     'FaceAlpha', 0.5, ...       % transparency
     'EdgeColor', 'none');       % no border


% ---- Crater Geometry: Period 1 ----
subplot(2,2,2);
max_depth1 = 100;
plot([-results1.best_radius,results1.best_radius],[0,0], 'b-','LineWidth',2); hold on;
plot([-results1.best_radius,-results1.best_radius],[0,max_depth1], 'b-','LineWidth',2)
plot([results1.best_radius,results1.best_radius],[0,max_depth1], 'b-','LineWidth',2)
plot([-results1.best_radius,results1.best_radius],[max_depth1,max_depth1], 'b-','LineWidth',2); 
grid on
axis equal
set(gca,'YDir','reverse')
xlabel('Radius (m)')
ylabel('Depth (m)')
% shaded magma region
xmin = -results1.best_radius;
xmax =  results1.best_radius;
ymin = min(results1.best_depths_timeseries);
ymax = max(results1.best_depths_timeseries);
fill([xmin xmax xmax xmin], ...  % x coordinates
     [ymin ymin ymax ymax], ...  % y coordinates
     [0.7 0.85 1], ...          % light blue color
     'FaceAlpha', 0.5, ...       % transparency
     'EdgeColor', 'none');       % no border
title('(b) Period 1: Geometry')


% ---- Crater Geometry: Period 2 ----
subplot(2,2,4);
max_depth2 = 100;
plot([-results2.best_radius,results2.best_radius],[0,0], 'r-','LineWidth',2); hold on;
plot([-results2.best_radius,-results2.best_radius],[0,max_depth2], 'r-','LineWidth',2)
plot([results2.best_radius,results2.best_radius],[0,max_depth2], 'r-','LineWidth',2)
plot([-results2.best_radius,results2.best_radius],[max_depth2,max_depth2], 'r-','LineWidth',2); 
grid on
axis equal
set(gca,'YDir','reverse')
xlabel('Radius (m)')
ylabel('Depth (m)')

% shaded magma region
xmin = -results2.best_radius;
xmax =  results2.best_radius;
ymin = min(results2.best_depths_timeseries);
ymax = max(results2.best_depths_timeseries);
fill([xmin xmax xmax xmin], ...  % x coordinates
     [ymin ymin ymax ymax], ...  % y coordinates
     [1 0.8 0.8], ...          % light red color
     'FaceAlpha', 0.5, ...       % transparency
     'EdgeColor', 'none');       % no border
title('(d) Period 2: Geometry')

%% ===================== FIGURE 2 =====================
% Observed vs modeled freqs, misfit, depth

figure(2); clf;
set(gcf,'Position',[100 100 1200 1000]);

colors = lines(3);
t1 = results1.peak_data.time;
t2 = results2.peak_data.time;

% ---- (a) Observed vs Modeled Freqs ----
subplot(2,1,1); hold on;

for mode = 1:3
    % Period 1
    f_obs_1 = results1.peak_data{:,mode+1};
    f_mod_1 = results1.modeled_freqs(:,mode);
    valid_1 = ~isnan(f_obs_1);
    scatter(t1(valid_1), f_obs_1(valid_1), 200, colors(mode,:), ...
        'filled','MarkerFaceAlpha',0.25,'MarkerEdgeColor','none')
    scatter(t1(valid_1), f_mod_1(valid_1), 120, colors(mode,:), ...
        's','filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.5)

    % Period 2
    f_obs_2 = results2.peak_data{:,mode+1};
    f_mod_2 = results2.modeled_freqs(:,mode);
    valid_2 = ~isnan(f_obs_2);
    scatter(t2(valid_2), f_obs_2(valid_2), 200, colors(mode,:), ...
        'filled','MarkerFaceAlpha',0.25,'HandleVisibility','off','MarkerEdgeColor','none')
    scatter(t2(valid_2), f_mod_2(valid_2), 120, colors(mode,:), ...
        's','filled','MarkerEdgeColor','k','HandleVisibility','off','MarkerFaceAlpha',0.5)
end

ylabel('Frequency (Hz)')
title('(a) Observed and Modeled Frequencies')
grid on, box on
legend('f_1 observed','f_1 model','f_2 observed','f_2 model','f_3 observed','f_3 model',...
    'Location','NorthWest');

% ---- (b) Depth Timeseries ----

subplot(2,1,2); hold on;
scatter(t1, results1.best_depths_timeseries, 120, colors(1,:), 's', ...
        'filled','MarkerFaceAlpha',0.5,'HandleVisibility','off','MarkerEdgeColor','k')
scatter(t2, results2.best_depths_timeseries, 120, colors(1,:), 's',  ...
        'filled','MarkerFaceAlpha',0.5,'HandleVisibility','off','MarkerEdgeColor','k')
set(gca,'YDir','reverse')
ylabel('Depth (m)')
title('(b) Inferred Magma Depth Over Time')
grid on, box on
ylim([15 90])

t1_start = results1.peak_data.time(1);
t1_end = results1.peak_data.time(end);
t2_start = results2.peak_data.time(1);
t2_end = results2.peak_data.time(end);

fill([t1_start t1_end t1_end t1_start], ...
    [15 15 25 25],...
    [0.7 0.85 1], ...          % light blue color
     'FaceAlpha', 0.5, ...       % transparency
     'EdgeColor', 'none');       % no border

fill([t2_start t2_end t2_end t2_start], ...
    [15 15 25 25],...
    [1 0.8 0.8], ...          % light blue color
     'FaceAlpha', 0.5, ...       % transparency
     'EdgeColor', 'none');       % no border

text(mean([t1_start, t1_end]), 20, ...
     sprintf('Period 1'), ...   
     'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'middle', ...
     'FontWeight', 'bold', ...
     'Color', 'b', 'FontSize', 14);

text(mean([t2_start, t2_end]), 20, ...
     sprintf('Period 2'), ...   
     'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'middle', ...
     'FontWeight', 'bold', ...
     'Color', 'r', 'FontSize', 14);
