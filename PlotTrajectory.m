%% PlotTrajectory.m
% Plots trajectory and velocity figures for a single bacterium using data
% from SimpleTrackingoutput.mat and pre-computed variables from
% CellTrackAnalysis.m.

% Last updated: 4/11/2026
% By Clara Shin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MAIN PARAMETERS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bacteria_idx = 1;           % column index of the bacterium to plot (integer)

%% =========================================================
%  Load metadata from SimpleTrackingoutput.mat
%  =========================================================

mat_data=       load('SimpleTrackingoutput.mat');   % Load the data calculated from CellTrackAnalysis.m
secondVideo=    10;                                 % video duration in seconds (default)
Nframes=        double(200);                        % total number of frames
time_int=       secondVideo / Nframes;              % time interval


%% =========================================================
%  Extract single-bacterium data from workspace variables
%  =========================================================
% These cell/matrix variables must already be in the workspace from running
% Cell_Track_Analysis_Algorithm.m beforehand.

% Raw x/y for this bacterium (may contain NaN for missing frames)
x_raw = xmat(:, bacteria_idx);
y_raw = ymat(:, bacteria_idx);
duration    = secondVideo;                    % total video duration (s)
frame_count = Nframes;                        % alias kept for clarity

% Frame numbers (1-based row indices) where the bacterium was detected
valid_rows = find(~isnan(x_raw));
frame      = valid_rows;
x          = x_raw(valid_rows);
y          = y_raw(valid_rows);

% t_x and t_y are left empty as requested
t_x = [];
t_y = [];

% Tumble vertex coordinates — one (x,y) point per tumble
t_x_v = tumble_x_mat{bacteria_idx};
t_y_v = tumble_y_mat{bacteria_idx};

% Velocity arrays computed by Cell_Track_Analysis_Algorithm.m
v          = v_mat{bacteria_idx};
w          = w_mat{bacteria_idx};
v_filt     = v_filt_mat{bacteria_idx};
w_filt     = w_filt_mat{bacteria_idx};
w_filt_off = w_filt_off_mat{bacteria_idx};

% Per-bacterium summary metrics
bact_avg_vel      = avg_vel(bacteria_idx);    % average run speed (µm/s) from Cell_Track
bact_avg_tumble_t = avg_tumble_t(bacteria_idx);
bact_avg_run_t    = avg_run_t(bacteria_idx);
bact_avg_angle    = avg_angle(bacteria_idx);

%% =========================================================
%  FIGURE 2(a) — Bacterial trajectory
%  =========================================================

has_tumbles = ~isempty(t_x_v) && ~isempty(t_y_v);

figure('Units','inches', 'Position',[1 1 8 8], 'Color','white');

plot(x, y, 'k-', 'LineWidth', 1);
hold on;

if has_tumbles
    scatter(t_x_v, t_y_v, 30, 'b', 'filled');
end

scatter(x(1), y(1), 200, 'g', 'o', 'LineWidth', 1.5);

xlabel('X (\mum)', 'FontSize', 18);
ylabel('Y (\mum)', 'FontSize', 18);
set(gca, 'YDir','reverse', 'TickDir','in', 'Box','on', 'FontSize', 14);

x_lo = min(x);  x_hi = max(x); x_range = x_hi - x_lo;
y_lo = min(y);  y_hi = max(y); y_range = y_hi - y_lo;
pad       = 0.15;                          
max_range = max(x_range, y_range);

x_mid = (x_lo + x_hi) / 2;
y_mid = (y_lo + y_hi) / 2;
half  = max_range / 2;

xlim([x_mid - half*(1 + pad), x_mid + half*(1 + pad)]);
ylim([y_mid - half*(1 + pad), y_mid + half*(1 + pad)]);
axis square;
hold off;

%% =========================================================
%  FIGURES 2(b)-(d) — Velocity & angular velocity subplots
%  =========================================================

frames     = (0 : length(v)-1)';   % 0-based frame index for x-axis
n_frames   = length(frames);

lightcoral = [1.00, 0.63, 0.63];
darkred    = [0.55, 0.00, 0.00];
grey       = [0.50, 0.50, 0.50];

t_lines = t_tumble_mat{bacteria_idx};

lbl_fs  = 14;
tick_fs = 13;

figure('Units','inches', 'Position',[1 1 9 8], 'Color','white');

% ----- Subplot B: translational velocity -----
ax1 = subplot(3,1,1);
plot(frames, v,      'Color',grey, 'LineWidth',1); hold on;
plot(frames, v_filt, 'k-',         'LineWidth',1);
xlim([0 n_frames]);
v_top = max(max(v), max(v_filt)) * 1.1;
ylim([0, v_top]);
apply_yticks(ax1, v_top);
yl1 = ylabel('Velocity (\mum/s)', 'FontSize',lbl_fs);
set(ax1, 'FontSize',tick_fs);
hold off;

% ----- Subplot C: angular velocity -----
ax2 = subplot(3,1,2);
plot(frames, w,      'Color',lightcoral, 'LineWidth',1); hold on;
plot(frames, w_filt, 'Color',darkred,    'LineWidth',1);
xlim([0 n_frames]);

w_top = max(max(w), max(w_filt)) * 1.1;
ylim([0, w_top]);
apply_yticks(ax2, w_top);
yl2 = ylabel('Angular Velocity (s^{-1})', 'FontSize',lbl_fs);
set(ax2, 'FontSize',tick_fs);
hold off;

% ----- Subplot D: overlay — angular velocity on top (red), velocity below (black)
ax3 = subplot(3,1,3);

v_all = max([v; v_filt], 0);
w_all = max([w; w_filt], 0);

v_min_val = min(v_all);  v_max_val = max(v_all);
w_min_val = min(w_all);  w_max_val = max(w_all);

if v_max_val == v_min_val, v_max_val = v_min_val + 1; end
if w_max_val == w_min_val, w_max_val = w_min_val + 1; end

v_filt_norm = (max(v_filt, 0) - v_min_val) / (v_max_val - v_min_val);
w_filt_norm = (max(w_filt, 0) - w_min_val) / (w_max_val - w_min_val);

gap    = 0.15;
pad_d  = 0.05;
v_plot = v_filt_norm * (0.5 - gap/2 - pad_d) + pad_d;
w_plot = w_filt_norm * (0.5 - gap/2 - pad_d) + (0.5 + gap/2);

hold on;
plot(frames, w_plot, 'Color',darkred, 'LineWidth',1); 
plot(frames, v_plot, 'k-', 'LineWidth',1);

for xi = t_lines'
    xline(xi, 'b-', 'LineWidth',1, 'Alpha',1);
end

xlim([0 n_frames]);
ylim([0 1]);

set(ax3, 'YTick',[], 'FontSize',tick_fs, 'Box','on');
yl3L = ylabel('Velocity (\mum/s)', 'FontSize',lbl_fs, 'Color','k');
xlabel('Frame Number', 'FontSize',lbl_fs);

yyaxis right;
set(ax3, 'YTick',[], 'YColor','k');
yl3R = ylabel('Angular Velocity (s^{-1})', 'FontSize',lbl_fs, 'Color',[0.55 0 0]);
yyaxis left;

hold off;

drawnow;
min_x = min([yl1.Position(1), yl2.Position(1), yl3L.Position(1)]);
yl1.Position(1) = min_x;
yl2.Position(1) = min_x;
yl3L.Position(1) = min_x;

%% =========================================================
%  Helper function — smart y-axis tick marks
%  =========================================================

function apply_yticks(ax, top_val)
    if top_val > 2
        fmt = @(v) sprintf('%d', round(v));
    elseif top_val > 0.15
        fmt = @(v) sprintf('%.1f', v);
    elseif top_val > 0.03
        fmt = @(v) sprintf('%.2f', v);
    else
        fmt = @(v) sprintf('%.1e', v);
    end

    for n_intervals = 4 : -1 : 3
        step   = top_val / n_intervals;
        ticks  = (0 : n_intervals) * step;
        labels = arrayfun(fmt, ticks, 'UniformOutput', false);

        if length(unique(labels)) == length(labels)
            yticks(ax, ticks);
            yticklabels(ax, labels);
            return;
        end
    end

    ticks  = [0, top_val/2, top_val];
    labels = arrayfun(fmt, ticks, 'UniformOutput', false);
    yticks(ax, ticks);
    yticklabels(ax, labels);
end