%% plot misfit %%

clear all; clc;
set(0,'DefaultAxesFontSize',14)

% colors
% c1 = [0.10 0.55 0.70];   % teal-blue
% c2 = [0.80 0.40 0.10];   % warm brown-orange
% c3 = [0.40 0.70 0.20];   % olive green

Paired = [
    166 206 227
    31 120 180
    178 223 138
    51 160 44
    251 154 153
    227 26 28
    253 191 111
    255 127 0
    202 178 214
    106 61 154
    255 255 153
    177 89 40
] / 255;



c1 = Paired(1,:);
c2 = Paired(3,:);
c3 = Paired(7,:);

c1a = Paired(2,:);
c2a = Paired(4,:);
c3a = Paired(8,:);

%% cylinder
load('inv_results_cylindrical_pipe.mat');
cyl1 = combined_results.section1;
cyl2 = combined_results.section2;

figure(1); clf; 

subplot(3,1,1); hold on;

% ---------- PERIOD 1 ----------
R1 = cyl1.peak_data;
% Add model frequencies
R1.f1_model = cyl1.modeled_freqs(:,1);
R1.f2_model = cyl1.modeled_freqs(:,2);
R1.f3_model = cyl1.modeled_freqs(:,3);
% Mark period
R1.period = ones(height(R1),1);
% Observed + model matrices
Fobs1 = [R1.f1, R1.f2, R1.f3];
Fmod1 = [R1.f1_model, R1.f2_model, R1.f3_model];
% Compute misfit only where obs exists
valid1 = ~isnan(Fobs1);
sqerr1 = zeros(size(Fobs1));
sqerr1(valid1) = (Fobs1(valid1) - Fmod1(valid1)).^2;
Nvalid1 = sum(valid1,2);
R1.misfit = sum(sqerr1,2) ./ Nvalid1;

% ---------- PERIOD 2 ----------
R2 = cyl2.peak_data;
% Add model frequencies
R2.f1_model = cyl2.modeled_freqs(:,1);
R2.f2_model = cyl2.modeled_freqs(:,2);
R2.f3_model = cyl2.modeled_freqs(:,3);
% Mark period
R2.period = 2*ones(height(R2),1);
% Observed + model matrices
Fobs2 = [R2.f1, R2.f2, R2.f3];
Fmod2 = [R2.f1_model, R2.f2_model, R2.f3_model];
% Compute misfit only where obs exists
valid2 = ~isnan(Fobs2);
sqerr2 = zeros(size(Fobs2));
sqerr2(valid2) = (Fobs2(valid2) - Fmod2(valid2)).^2;
Nvalid2 = sum(valid2,2);
R2.misfit = sum(sqerr2,2) ./ Nvalid2;

scatter(R1.time, R1.misfit, 200, 'k', 'filled', 'MarkerFaceAlpha', 0.5, ...
    'MarkerEdgeColor','k');
scatter(R2.time, R2.misfit, 200, 'k', 'filled', 'MarkerFaceAlpha', 0.5, ...
    'MarkerEdgeColor','k'); 
set(gca, 'YScale', 'log');   % <-- semilogy effect
grid on
xlabel('Time')
ylabel('Misfit')
title('Misfit vs Time (log scale)')
ylim([10e-6, 10])

% smoothed misfit comparison
figure(2); clf;
set(gcf,'Position',[100 100 1200 600]);
% Convert times to numeric if needed
t1 = datenum(R1.time);
t2 = datenum(R2.time);
% 24-hour window in units of days (since datenum uses days)
window = 24/24;   % = 1 day
% Moving average of misfit values
m1_smooth = movmean(R1.misfit, [12 11]);   % or movmean(R1.misfit, 24)
m2_smooth = movmean(R2.misfit, [12 11]);
% Plot original points (optional, lightly)
scatter(R1.time, R1.misfit, 60, c1, 'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor',c1,'HandleVisibility','off');
hold on
scatter(R2.time, R2.misfit, 60, c1, 'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor',c1);
% Plot smoothed curves
plot(R1.time, m1_smooth, 'Color', c1a, 'LineWidth', 3,'HandleVisibility','off','LineStyle',':');
plot(R2.time, m2_smooth, 'Color', c1a, 'LineWidth', 3,'LineStyle',':');
set(gca,'YScale','log')
ylim([10e-5 10e0]);
box on; grid on;

mean(R1.misfit)
mean(R2.misfit)

%% analytical conical frustum

figure(1);
subplot(3,1,2); hold on;

load('inv_results_analytical_conical_frustum.mat');
num_frustum1 = combined_results.section1;
num_frustum2 = combined_results.section2;

% ---------- PERIOD 1 ----------
R1 = num_frustum1.peak_data;
% Add model frequencies
R1.f1_model = num_frustum1.modeled_freqs(:,1);
R1.f2_model = num_frustum1.modeled_freqs(:,2);
R1.f3_model = num_frustum1.modeled_freqs(:,3);
% Mark period
R1.period = ones(height(R1),1);
% Observed + model matrices
Fobs1 = [R1.f1, R1.f2, R1.f3];
Fmod1 = [R1.f1_model, R1.f2_model, R1.f3_model];
% Compute misfit only where obs exists
valid1 = ~isnan(Fobs1);
sqerr1 = zeros(size(Fobs1));
sqerr1(valid1) = (Fobs1(valid1) - Fmod1(valid1)).^2;
Nvalid1 = sum(valid1,2);
R1.misfit = sum(sqerr1,2) ./ Nvalid1;

% ---------- PERIOD 2 ----------
R2 = num_frustum2.peak_data;
% Add model frequencies
R2.f1_model = num_frustum2.modeled_freqs(:,1);
R2.f2_model = num_frustum2.modeled_freqs(:,2);
R2.f3_model = num_frustum2.modeled_freqs(:,3);
% Mark period
R2.period = 2*ones(height(R2),1);
% Observed + model matrices
Fobs2 = [R2.f1, R2.f2, R2.f3];
Fmod2 = [R2.f1_model, R2.f2_model, R2.f3_model];
% Compute misfit only where obs exists
valid2 = ~isnan(Fobs2);
sqerr2 = zeros(size(Fobs2));
sqerr2(valid2) = (Fobs2(valid2) - Fmod2(valid2)).^2;
Nvalid2 = sum(valid2,2);
R2.misfit = sum(sqerr2,2) ./ Nvalid2;

scatter(R1.time, R1.misfit, 120, 'r', 'filled', 'MarkerFaceAlpha', 0.5, ...
    'MarkerEdgeColor','k');
scatter(R2.time, R2.misfit, 120, 'r', 'filled', 'MarkerFaceAlpha', 0.5, ...
    'MarkerEdgeColor','k'); 
set(gca, 'YScale', 'log');   % <-- semilogy effect
grid on
xlabel('Time')
ylabel('Misfit')
title('Misfit vs Time (log scale)')
ylim([10e-6, 10]);

% smoothed misfit comparison
figure(2); 
% Convert times to numeric if needed
t1 = datenum(R1.time);
t2 = datenum(R2.time);
% 24-hour window in units of days (since datenum uses days)
window = 24/24;   % = 1 day
% Moving average of misfit values
m1_smooth = movmean(R1.misfit, [12 11]);   % or movmean(R1.misfit, 24)
m2_smooth = movmean(R2.misfit, [12 11]);
% Plot original points (optional, lightly)
scatter(R1.time, R1.misfit, 100, c2, 'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor',c2);
hold on
scatter(R2.time, R2.misfit, 100, c2, 'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor',c2,'HandleVisibility','off');
% Plot smoothed curves
plot(R1.time, m1_smooth, 'Color', c2a, 'LineWidth', 3,'LineStyle','-.');
plot(R2.time, m2_smooth, 'Color', c2a, 'LineWidth', 3,'HandleVisibility','off','LineStyle','-.');

mean(R1.misfit)
mean(R2.misfit)

%% numerical conical frustum

figure(1);
subplot(3,1,3); hold on;

load('inv_results_numerical_conical_frustum.mat');
num_frustum1 = combined_results.section1;
num_frustum2 = combined_results.section2;


% ---------- PERIOD 1 ----------
R1 = num_frustum1.peak_data;
% Add model frequencies
R1.f1_model = num_frustum1.modeled_freqs(:,1);
R1.f2_model = num_frustum1.modeled_freqs(:,2);
R1.f3_model = num_frustum1.modeled_freqs(:,3);
% Mark period
R1.period = ones(height(R1),1);
% Observed + model matrices
Fobs1 = [R1.f1, R1.f2, R1.f3];
Fmod1 = [R1.f1_model, R1.f2_model, R1.f3_model];
% Compute misfit only where obs exists
valid1 = ~isnan(Fobs1);
sqerr1 = zeros(size(Fobs1));
sqerr1(valid1) = (Fobs1(valid1) - Fmod1(valid1)).^2;
Nvalid1 = sum(valid1,2);
R1.misfit = sum(sqerr1,2) ./ Nvalid1;

% ---------- PERIOD 2 ----------
R2 = num_frustum2.peak_data;
% Add model frequencies
R2.f1_model = num_frustum2.modeled_freqs(:,1);
R2.f2_model = num_frustum2.modeled_freqs(:,2);
R2.f3_model = num_frustum2.modeled_freqs(:,3);
% Mark period
R2.period = 2*ones(height(R2),1);
% Observed + model matrices
Fobs2 = [R2.f1, R2.f2, R2.f3];
Fmod2 = [R2.f1_model, R2.f2_model, R2.f3_model];
% Compute misfit only where obs exists
valid2 = ~isnan(Fobs2);
sqerr2 = zeros(size(Fobs2));
sqerr2(valid2) = (Fobs2(valid2) - Fmod2(valid2)).^2;
Nvalid2 = sum(valid2,2);
R2.misfit = sum(sqerr2,2) ./ Nvalid2;

scatter(R1.time, R1.misfit, 120, 'b', 'filled', 'MarkerFaceAlpha', 0.5, ...
    'MarkerEdgeColor','k');
scatter(R2.time, R2.misfit, 120, 'b', 'filled', 'MarkerFaceAlpha', 0.5, ...
    'MarkerEdgeColor','k'); 
set(gca, 'YScale', 'log');   % <-- semilogy effect
grid on
xlabel('Time')
ylabel('Misfit')
title('Misfit vs Time (log scale)')
ylim([10e-6, 10]);

% smoothed misfit comparison
figure(2); 
% Convert times to numeric if needed
t1 = datenum(R1.time);
t2 = datenum(R2.time);
% 24-hour window in units of days (since datenum uses days)
window = 24/24;   % = 1 day
% Moving average of misfit values
m1_smooth = movmean(R1.misfit, [12 11]);   % or movmean(R1.misfit, 24)
m2_smooth = movmean(R2.misfit, [12 11]);
% Plot original points (optional, lightly)
scatter(R1.time, R1.misfit, 100, c3, 'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor',c3,'HandleVisibility','off');
hold on
scatter(R2.time, R2.misfit, 100, c3, 'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor',c3);
% Plot smoothed curves
plot(R1.time, m1_smooth, 'Color', c3a, 'LineWidth', 3,'HandleVisibility','off','LineStyle','-');
plot(R2.time, m2_smooth, 'Color', c3a, 'LineWidth', 3,'LineStyle','-');

ylabel('Misfit')
legend('Cylindrical Pipe','Cylindrical Pipe: 24 hour average',...
    'Analytical Conical Frustum','Analytical Conical Frustum: 24 hour average',...
    'Numerical Conical Frustum','Numerical Conical Frustum: 24 hour average',...
    'Location','NorthWest')
mean(R1.misfit)
mean(R2.misfit)