% CellTrackAnalysis.m
% Reference: Johnson et al., (2024)

% Calculates number of tumbles and velocities from bacteria images tracked 
% from SimpleTracking.m, using a custom algorithm based on the work of 
% Johnson et al. 

% From the tracked x and y coordinates, velocity and angular velocity 
% based on s per frame is calculated and the number, average velocities and
% degrees of tumbles and runs.

% Last updated: 4/14/2026
% By Clara Shin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MAIN PARAMETERS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

secondVideo=        30;                     % Length of Video (sec)

depth_speed_ratio=  0.7;                    % Ratio of Delta v (Depth) and minimum v ratio
v_min_threshold=    0.3;                    % vmin threshold to be included in tumble count

sg_window=          5;                      % Savitzky-Golay filter parameter. Must be odd
sg_poly=            2;                      % Savitzky-Golay filter parameter for smoothing velocities

longRunThreshold=   6;                      % frames; runs of >= this length counted as "long"
minTrackLength=     sg_window + 2;          % Minimum track length to attempt analysis (must be > sg_window)
w_offset=           50;                     % rad/s offset added to w_filt for overlay plot

%% =========================================================
%  Load simpletracking output data
%  =========================================================

data      = load('SimpleTrackingoutput.mat');
objs_link = data.objs_link;

numFrames = double(data.Nframes);

time_int = secondVideo / numFrames;

x_raw     = objs_link(1, :)';
y_raw     = objs_link(2, :)';
frame_raw = objs_link(5, :)';
id_raw    = objs_link(6, :)';

track_ids  = unique(id_raw);
Nbacteria  = length(track_ids);

xmat = NaN(numFrames, Nbacteria);
ymat = NaN(numFrames, Nbacteria);

for b = 1 : Nbacteria
    mask          = (id_raw == track_ids(b));
    frs           = frame_raw(mask);
    xmat(frs, b)  = x_raw(mask);
    ymat(frs, b)  = y_raw(mask);
end

%% =========================================================
%  Allocate output cell arrays and summary metric vectors
%  =========================================================

t_x_mat           = cell(1, Nbacteria);
t_y_mat           = cell(1, Nbacteria);
t_x_v_mat         = cell(1, Nbacteria);
t_y_v_mat         = cell(1, Nbacteria);
v_mat             = cell(1, Nbacteria);
w_mat             = cell(1, Nbacteria);
v_filt_mat        = cell(1, Nbacteria);
w_filt_mat        = cell(1, Nbacteria);
w_filt_off_mat    = cell(1, Nbacteria);
tumble_frames_mat = cell(1, Nbacteria);
t_tumble_mat      = cell(1, Nbacteria);
tumble_x_mat      = cell(1, Nbacteria);
tumble_y_mat      = cell(1, Nbacteria);

num_tumbles       = NaN(1, Nbacteria);
avg_angle         = NaN(1, Nbacteria);
avg_tumble_t      = NaN(1, Nbacteria);
avg_run_t         = NaN(1, Nbacteria);
avg_vel           = NaN(1, Nbacteria);   % avg translational speed during runs (µm/s)
avg_v_tumble      = NaN(1, Nbacteria);   % avg translational speed during tumbles (µm/s)
avg_w_run         = NaN(1, Nbacteria);   % avg angular velocity during runs (rad/s)
avg_w_tumble      = NaN(1, Nbacteria);   % avg angular velocity during tumbles (rad/s)
track_duration_s  = NaN(1, Nbacteria);   % total track duration in seconds

%% =========================================================
%  Bacterium metrics (translational velocity and angular velocity)
%  =========================================================

for b = 1 : Nbacteria

    % ---- Extract valid (non-NaN) rows for this bacterium ----
    x_col = xmat(:, b);
    y_col = ymat(:, b);
    valid = ~isnan(x_col) & ~isnan(y_col);

    if sum(valid) < minTrackLength
        continue
    end

    x_v = x_col(valid);
    y_v = y_col(valid);

    % ---- Velocities ----
    vx = diff(x_v) / time_int;
    vy = diff(y_v) / time_int;
    v_full     = sqrt(vx.^2 + vy.^2);
    phi_full   = acos(clip(vx ./ v_full));
    omega_full = abs(diff(phi_full)) / time_int;

    % Align arrays: drop first 2 points so all have length N-2
    x_out  = x_v(3:end);
    y_out  = y_v(3:end);
    vx_out = vx(2:end);
    vy_out = vy(2:end);
    v_out  = v_full(2:end);
    w_out  = omega_full;

    if length(x_out) < sg_window
        continue
    end

    % ---- Savitzky-Golay filtering ----
    sg_v = sgolayfilt(v_out, sg_poly, sg_window);
    sg_w = sgolayfilt(w_out, sg_poly, sg_window);

    % ---- Store trajectory arrays ----
    t_x_mat{b}        = x_out;
    t_y_mat{b}        = y_out;
    t_x_v_mat{b}      = vx_out;
    t_y_v_mat{b}      = vy_out;
    v_mat{b}          = v_out;
    w_mat{b}          = w_out;
    v_filt_mat{b}     = sg_v;
    w_filt_mat{b}     = sg_w;
    w_filt_off_mat{b} = sg_w + w_offset;

    % ----------------------------------------------------------------
    %  TUMBLE DETECTION  (4 cases)
    % ----------------------------------------------------------------

    N_sg     = length(sg_v);
    min_mask = islocalmin(sg_v, 'FlatSelection', 'first');
    max_mask = islocalmax(sg_v, 'FlatSelection', 'first');
    min_mask([1, N_sg]) = false;
    max_mask([1, N_sg]) = false;
    min_idx  = find(min_mask);
    max_idx  = find(max_mask);
    n_min    = length(min_idx);
    n_max    = length(max_idx);

    if n_min < 2 || n_max < 1
        num_tumbles(b) = 0;
        continue
    end

    v_min_all = sg_v(min_idx);
    v_max_all = sg_v(max_idx);

    delta_v1 = [];  delta_v2 = [];  delta_v = [];
    v_tumble = [];  t_tumble = [];
    v_min_used = [];

    % ---- Case 1: equal counts, min comes first ----
    if n_min == n_max && min_idx(1) < max_idx(1)
        for i = 1 : n_min-1
            delta_v1(end+1) = v_max_all(i+1) - v_min_all(i+1);
            delta_v2(end+1) = v_max_all(i)   - v_min_all(i+1);
        end
        v_min_used = v_min_all(2:end);
        for i = 1 : length(delta_v1)
            delta_v(end+1) = max(delta_v1(i), delta_v2(i));
        end
        v_depth_speed_ratio = delta_v ./ v_min_used';
        for i = 1 : length(v_depth_speed_ratio)
            if v_depth_speed_ratio(i) > depth_speed_ratio
                v_tumble(end+1) = v_min_used(i);
                t_tumble(end+1) = min_idx(i+1);
            end
        end

    % ---- Case 2: equal counts, max comes first ----
    elseif n_min == n_max && min_idx(1) > max_idx(1)
        for i = 1 : n_min-1
            delta_v1(end+1) = v_max_all(i+1) - v_min_all(i);
            delta_v2(end+1) = v_max_all(i)   - v_min_all(i);
        end
        v_min_used = v_min_all(1:end-1);
        for i = 1 : length(delta_v1)
            delta_v(end+1) = max(delta_v1(i), delta_v2(i));
        end
        v_depth_speed_ratio = delta_v ./ v_min_used';
        for i = 1 : length(v_depth_speed_ratio)
            if v_depth_speed_ratio(i) > depth_speed_ratio
                v_tumble(end+1) = v_min_used(i);
                t_tumble(end+1) = min_idx(i);
            end
        end

    % ---- Case 3: more minima than maxima ----
    elseif n_min > n_max
        for i = 1 : n_min-2
            delta_v1(end+1) = v_max_all(i+1) - v_min_all(i+1);
            delta_v2(end+1) = v_max_all(i)   - v_min_all(i+1);
        end
        v_min_used = v_min_all(2:end-1);
        for i = 1 : length(delta_v1)
            delta_v(end+1) = max(delta_v1(i), delta_v2(i));
        end
        v_depth_speed_ratio = delta_v ./ v_min_used';
        for i = 1 : length(v_depth_speed_ratio)
            if v_depth_speed_ratio(i) > depth_speed_ratio
                v_tumble(end+1) = v_min_used(i);
                t_tumble(end+1) = min_idx(i+1);
            end
        end

    % ---- Case 4: more maxima than minima ----
    elseif n_min < n_max
        for i = 1 : n_min
            delta_v1(end+1) = v_max_all(i+1) - v_min_all(i);
            delta_v2(end+1) = v_max_all(i)   - v_min_all(i);
        end
        v_min_used = v_min_all;
        for i = 1 : length(delta_v1)
            delta_v(end+1) = max(delta_v1(i), delta_v2(i));
        end
        v_depth_speed_ratio = delta_v ./ v_min_used';
        for i = 1 : length(v_depth_speed_ratio)
            if abs(v_depth_speed_ratio(i)) > depth_speed_ratio
                v_tumble(end+1) = v_min_used(i);
                t_tumble(end+1) = min_idx(i);
            end
        end
    end

    n_tumble_detected = length(v_tumble);

    if n_tumble_detected == 0
        num_tumbles(b) = 0;
        continue
    end

    % ----------------------------------------------------------------
    %  ADJACENT MAXIMA
    % ----------------------------------------------------------------

    adj_max = [];

    if n_min == n_max && min_idx(1) < max_idx(1)          % Case 1
        for i = 1 : length(t_tumble)
            for j = 1 : n_min
                if t_tumble(i) == min_idx(j)
                    for k = 1 : n_max
                        if max_idx(k) < min_idx(j) && (j == 1 || max_idx(k) > min_idx(j-1))
                            adj_max(end+1) = max_idx(k);
                        end
                        if max_idx(k) > min_idx(j)
                            if n_min == j
                                adj_max(end+1) = max_idx(k);
                            elseif max_idx(k) < min_idx(j+1)
                                adj_max(end+1) = max_idx(k);
                            end
                        end
                    end
                end
            end
        end

    elseif n_min == n_max && min_idx(1) > max_idx(1)      % Case 2
        for i = 1 : length(t_tumble)
            for j = 1 : n_min
                if t_tumble(i) == min_idx(j)
                    for k = 1 : n_max
                        if max_idx(k) < min_idx(j)
                            if j == 1 || max_idx(k) > min_idx(j-1)
                                adj_max(end+1) = max_idx(k);
                            end
                        end
                        if max_idx(k) > min_idx(j) && max_idx(k) < min_idx(j+1)
                            adj_max(end+1) = max_idx(k);
                        end
                    end
                end
            end
        end

    elseif n_min > n_max                                   % Case 3
        for i = 1 : length(t_tumble)
            for j = 1 : n_min
                if t_tumble(i) == min_idx(j)
                    for k = 1 : n_max
                        if max_idx(k) < min_idx(j) && (j == 1 || max_idx(k) > min_idx(j-1))
                            adj_max(end+1) = max_idx(k);
                        end
                        if max_idx(k) > min_idx(j) && max_idx(k) < min_idx(j+1)
                            adj_max(end+1) = max_idx(k);
                        end
                    end
                end
            end
        end

    elseif n_min < n_max                                   % Case 4
        for i = 1 : length(t_tumble)
            for j = 1 : n_min
                if t_tumble(i) == min_idx(j)
                    for k = 1 : n_max
                        if max_idx(k) < min_idx(j)
                            if j == 1 || max_idx(k) > min_idx(j-1)
                                adj_max(end+1) = max_idx(k);
                            end
                        elseif max_idx(k) > min_idx(j)
                            if n_min == j || max_idx(k) < min_idx(j+1)
                                adj_max(end+1) = max_idx(k);
                            end
                        end
                    end
                end
            end
        end
    end

    adj_max = sort(adj_max);

    if length(adj_max) ~= 2 * n_tumble_detected
        num_tumbles(b) = n_tumble_detected;
        continue
    end

    % ----------------------------------------------------------------
    %  TUMBLE FRAMES, TUMBLE TIMES, RUN TIMES
    % ----------------------------------------------------------------

    tumble_time   = zeros(1, n_tumble_detected);
    run_time      = zeros(1, n_tumble_detected + 1);
    tumble_frames = [];
    run_frames    = [];
    N             = length(x_out);

    for i = 1 : n_tumble_detected
        for j = 1 : N
            if j >= adj_max(2*i-1) && j <= adj_max(2*i)
                for k = 1 : length(v_min_used)
                    if v_tumble(i) == v_min_used(k)
                        if sg_v(j) <= v_tumble(i) + v_min_threshold * delta_v(k)
                            if ~ismember(j, tumble_frames)
                                tumble_frames(end+1) = j;
                                tumble_time(i) = time_int + tumble_time(i);
                            end
                        end
                    end
                end
            end
        end
    end

    run_idx = 1;
    for j = 1 : N
        if ~ismember(j, tumble_frames)
            run_frames(end+1) = j;
            run_time(run_idx) = run_time(run_idx) + time_int;
        elseif j > 1 && ~ismember(j-1, tumble_frames)
            run_idx = run_idx + 1;
        else
            for k = 1 : length(adj_max)
                if adj_max(k) == j
                    if ismember(j, tumble_frames) && ...
                       j > 1 && ismember(j-1, tumble_frames) && ...
                       j < N && ismember(j+1, tumble_frames)
                        if run_idx > 1 && run_time(run_idx-1) ~= 0
                            run_time(run_idx) = 0;
                            run_idx = run_idx + 1;
                        end
                    end
                end
            end
        end
    end

    if length(run_time) >= 2
        run_time(1)   = [];
        run_time(end) = [];
    else
        run_time = [];
    end

    if isempty(run_time)
        num_tumbles(b) = n_tumble_detected;
        continue
    end

    % ----------------------------------------------------------------
    %  TUMBLE ANGLES
    % ----------------------------------------------------------------

    dif_ix_arr = [];  dif_iy_arr = [];
    dif_ox_arr = [];  dif_oy_arr = [];
    iphi_arr   = [];  ophi_arr   = [];
    tumble_angles_cell = {};

    for i = 1 : length(tumble_frames)
        h = tumble_frames(i);
        if h > 3 && ~ismember(h-1, tumble_frames)
            dif_ix = x_out(h-1) - x_out(h-3);
            dif_iy = y_out(h-1) - y_out(h-3);
            dif_ix_arr(end+1) = dif_ix;
            dif_iy_arr(end+1) = dif_iy;
            iv_x = dif_ix / time_int;
            iv_y = dif_iy / time_int;
            iv_mag = sqrt(iv_x^2 + iv_y^2);
            iphi_arr(end+1) = acos(clip(iv_x / iv_mag)) * (180/pi);
        end
        if h+2 <= N && ~ismember(h+1, tumble_frames)
            dif_ox = x_out(h+2) - x_out(h);
            dif_oy = y_out(h+2) - y_out(h);
            dif_ox_arr(end+1) = dif_ox;
            dif_oy_arr(end+1) = dif_oy;
            ov_x = dif_ox / time_int;
            ov_y = dif_oy / time_int;
            ov_mag = sqrt(ov_x^2 + ov_y^2);
            ophi_arr(end+1) = acos(clip(ov_x / ov_mag)) * (180/pi);
        end
    end

    n_skip = 0;
    for i = 1 : length(tumble_time)
        if i == 1
            rt_prev = run_time(end);
        elseif (i-1) <= length(run_time)
            rt_prev = run_time(i-1);
        else
            rt_prev = 0;
        end
        idx = i - n_skip;

        if idx < 1 || idx > length(dif_iy_arr) || idx > length(dif_oy_arr)
            continue
        end

        if rt_prev ~= 0
            tumble_angles_cell{end+1} = compute_angle( ...
                dif_iy_arr(idx), dif_oy_arr(idx), iphi_arr(idx), ophi_arr(idx));
        elseif rt_prev == 0 && (i-1) >= 1
            tumble_angles_cell{end+1} = 'See previous entry';
            n_skip = n_skip + 1;
        elseif rt_prev == 0 && (i-1) < 1
            tumble_angles_cell{end+1} = compute_angle( ...
                dif_iy_arr(idx), dif_oy_arr(idx), iphi_arr(idx), ophi_arr(idx));
        end
    end

    % ----------------------------------------------------------------
    %  SUMMARY METRICS
    % ----------------------------------------------------------------

    tumble_time_clean = tumble_time(tumble_time > 0);
    run_time_clean    = run_time(run_time > 0);

    angles_clean = [];
    for i = 1 : length(tumble_angles_cell)
        if isnumeric(tumble_angles_cell{i})
            angles_clean(end+1) = tumble_angles_cell{i};
        end
    end

    num_tumbles(b)  = length(tumble_time_clean);
    avg_angle(b)    = mean(angles_clean,      'omitnan');
    avg_tumble_t(b) = mean(tumble_time_clean, 'omitnan');
    avg_run_t(b)    = mean(run_time_clean,    'omitnan');

    % Average translational and angular speeds during tumble and run frames
    if ~isempty(tumble_frames)
        avg_v_tumble(b) = mean(sg_v(tumble_frames), 'omitnan');
        avg_w_tumble(b) = mean(sg_w(tumble_frames), 'omitnan');
    end
    if ~isempty(run_frames)
        avg_vel(b)    = mean(sg_v(run_frames), 'omitnan');
        avg_w_run(b)  = mean(sg_w(run_frames), 'omitnan');
    end

    % Track duration = number of valid frames * time_int
    track_duration_s(b) = sum(valid) * time_int;

    tumble_frames_mat{b} = tumble_frames;
    t_tumble_mat{b}      = t_tumble;

    valid_tt = t_tumble(t_tumble >= 1 & t_tumble <= length(x_out));
    tumble_x_mat{b} = x_out(valid_tt);
    tumble_y_mat{b} = y_out(valid_tt);

end   % end bacterium loop

%% =========================================================
%  Histograms
%  =========================================================

dodgerblue = [0.118, 0.565, 1.000];
avgline    = [0.85, 0.10, 0.10];

nt_plot   = num_tumbles(~isnan(num_tumbles));
att_plot  = avg_tumble_t(~isnan(avg_tumble_t));
art_plot  = avg_run_t(~isnan(avg_run_t));
aa_plot   = avg_angle(~isnan(avg_angle));
dur_plot  = track_duration_s(~isnan(track_duration_s));
vr_plot   = avg_vel(~isnan(avg_vel));           % translational speed during runs
vt_plot   = avg_v_tumble(~isnan(avg_v_tumble)); % translational speed during tumbles
wr_plot   = avg_w_run(~isnan(avg_w_run));       % angular velocity during runs
wt_plot   = avg_w_tumble(~isnan(avg_w_tumble)); % angular velocity during tumbles

% ---- Population-level summary statistics ----
pop_mean_nt  = mean(nt_plot);    pop_med_nt  = median(nt_plot);
pop_mean_att = mean(att_plot);   pop_med_att = median(att_plot);
pop_mean_art = mean(art_plot);   pop_med_art = median(art_plot);
pop_mean_aa  = mean(aa_plot);    pop_med_aa  = median(aa_plot);
pop_mean_dur = mean(dur_plot);   pop_med_dur = median(dur_plot);
pop_mean_vr  = mean(vr_plot);    pop_med_vr  = median(vr_plot);
pop_mean_vt  = mean(vt_plot);    pop_med_vt  = median(vt_plot);
pop_mean_wr  = mean(wr_plot);    pop_med_wr  = median(wr_plot);
pop_mean_wt  = mean(wt_plot);    pop_med_wt  = median(wt_plot);

% ---- Print summary table ----
fprintf('\n');
fprintf('=====================================================================\n');
fprintf('          Population Summary  (N = %d bacteria analysed)\n', length(nt_plot));
fprintf('=====================================================================\n');
fprintf('%-35s   %10s   %10s\n', 'Metric', 'Mean', 'Median');
fprintf('---------------------------------------------------------------------\n');
fprintf('%-35s   %10.2f   %10.2f\n', 'Number of tumbles', pop_mean_nt, pop_med_nt);
fprintf('%-35s   %10.2f   %10.2f\n', 'Tumble duration (s)', pop_mean_att, pop_med_att);
fprintf('%-35s   %10.2f   %10.2f\n', 'Run duration (s)', pop_mean_art, pop_med_art);
fprintf('%-35s   %10.2f   %10.2f\n', 'Tumble angle (deg)', pop_mean_aa, pop_med_aa);
fprintf('%-35s   %10.2f   %10.2f\n', 'Tumble linear speed (µm/s)', pop_mean_vt, pop_med_vt);
fprintf('%-35s   %10.2f   %10.2f\n', 'Run linear speed (µm/s)', pop_mean_vr, pop_med_vr);
fprintf('=====================================================================\n\n');

% ---- Histogram panels ----
figure('Units','inches', 'Position',[1 1 15 8], 'Color','white');

subplot(2,3,1); smart_histogram(nt_plot,  pop_mean_nt,  'Number of Tumbles',    'Tumble Count',        dodgerblue, avgline, true);
subplot(2,3,2); smart_histogram(att_plot, pop_mean_att, 'Tumble Duration (s)',   'Tumble Duration',     dodgerblue, avgline, false);
subplot(2,3,3); smart_histogram(art_plot, pop_mean_art, 'Run Duration (s)',      'Run Duration',        dodgerblue, avgline, false);
subplot(2,3,4); smart_histogram(aa_plot,  pop_mean_aa,  'Tumble Angle (deg)',    'Tumble Angle',        dodgerblue, avgline, false);
subplot(2,3,5); smart_histogram(vt_plot,  pop_mean_vt,  'Tumble Speed (\mum/s)', 'Tumble Linear Speed', dodgerblue, avgline, false);
subplot(2,3,6); smart_histogram(vr_plot,  pop_mean_vr,  'Run Speed (\mum/s)',    'Run Linear Speed',    dodgerblue, avgline, false);


%% =========================================================
%  Helper functions
%  =========================================================

function smart_histogram(data, mean_val, xlbl, ttl, dodgerblue, avgline, is_integer)
% Plots a histogram with auto-selected bin width, padded axes, and a mean line.

    if isempty(data) || all(isnan(data))
        title(ttl, 'FontSize',14);
        xlabel(xlbl, 'FontSize',14);
        return;
    end

    data = data(~isnan(data));
    data_min = min(data);
    data_max = max(data);
    data_range = data_max - data_min;

    % Guard: all values identical
    if data_range == 0
        data_range = max(1, abs(data_min));
    end

    if is_integer
        bw = max(1, round(data_range / 20));
        h  = histogram(data, 'BinWidth',bw, 'FaceColor',dodgerblue, 'EdgeColor','k');
    else
        h  = histogram(data, 15, 'FaceColor',dodgerblue, 'EdgeColor','k');
    end

    hold on;
    xline(mean_val, '--', 'Color',avgline, 'LineWidth',1.5, ...
        'Label', sprintf('Mean=%.2f', mean_val), ...
        'LabelOrientation','horizontal', 'FontSize',10);
    hold off;

    x_pad = data_range * 0.05;
    xlim([data_min - x_pad,  data_max + x_pad]);

    if ~isempty(h.Values) && max(h.Values) > 0
        ylim([0, max(h.Values) * 1.10]);
    end

    xlabel(xlbl,       'FontSize',14);
    ylabel('Frequency','FontSize',14);
    title(ttl,         'FontSize',14);
    set(gca, 'FontSize',12, 'TickDir','in', 'Box','on');
end

function val = clip(x)
% Clamp x to [-1, 1] to guard against floating-point errors in acos
    val = max(-1, min(1, x));
end

function angle = compute_angle(diy, doy, iphi, ophi)
% Compute tumble angle (deg) from in/out heading components
    if diy > 0 && doy > 0
        d = abs(ophi - iphi);
    elseif diy > 0 && doy < 0
        d = abs(-ophi - iphi);
    elseif diy < 0 && doy > 0
        d = abs(ophi + iphi);
    else
        d = abs(-ophi + iphi);
    end
    if d <= 180
        angle = round(d, 2);
    else
        angle = round(abs(360 - d), 2);
    end
end