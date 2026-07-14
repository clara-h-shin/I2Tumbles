%% CellTrackAnalysis.m
% Reference: Johnson et al., (2024)

% Calculates number of tumbles and velocities from bacteria images tracked 
% from SimpleTracking.m, using a custom algorithm based on the work of 
% Johnson et al. 

% From the tracked x and y coordinates, velocity and angular velocity 
% based on s per frame is calculated and the number, average velocities and
% degrees of tumbles and runs.

% < IMPORTANT: Please change the file names after running CellTrackAnalysis.m >
% ^ If you want to perform statistical testing. 

% Last updated: 7/13/2026
% By Clara Shin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MAIN PARAMETERS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

secondVideo=        10;                     % Length of Video (sec)

% --- Pixel-to-micron conversion ---
convert_to_um =     true;                   
px_to_um =          0.33;                   % µm per pixel

% --- Savitzky-Golay smoothing ---
fps_scale_sg_window= true;                  % true = auto-scale | false = use fixed value
sg_window=          5;                      % Odd numbers only. Used when fps_scale_sg_window = false
sg_poly=            2;                      % polynomial order (keep at 2)

minTrackLength=     sg_window + 2;          % minimum track length; updated below if scaling
w_offset=           50;                     % rad/s offset added to w_filt for overlay plot

% --- Tumble detection thresholds ---
depth_speed_ratio=  0.7;                    % Value: 0.5 - 0.999 (The higher, more conservative)
v_min_threshold=    0.3;                    % tumble frame width threshold

% --- Angular speed confirmation ---
use_angular_confirmation = true;            
w_confirm_factor  = 0.20;            

% --- Minimum speed at tumble dip ---
min_tumble_speed  = 3.0;                    

% --- Scatter plot: show bacteria IDs next to each dot ---
show_scatter_labels = false;               % true = annotate each dot with its TrackID

% --- Display toggle: set to false to skip all histograms and scatter plot ---
show_plots = true;

%% =========================================================
%  Load simpletracking output data
%  =========================================================

data      = load('SimpleTrackingoutput.mat');
objs_link = data.objs_link;

numFrames = double(data.Nframes);

time_int = secondVideo / numFrames;

% --- Scale sg_window to match the paper's ~200 ms smoothing window ---
if fps_scale_sg_window
    fps_data    = numFrames / secondVideo;
    fps_paper   = 25;
    paper_window_ms = 5 / fps_paper * 1000;           % = 200 ms
    equiv_frames    = round(paper_window_ms / 1000 * fps_data);
    if mod(equiv_frames, 2) == 0, equiv_frames = equiv_frames + 1; end  % must be odd
    sg_window   = max(5, equiv_frames);                % never go below paper's value
    fprintf('FPS=%.1f: sg_window auto-scaled to %d (%.0f ms window)\n', ...
        fps_data, sg_window, sg_window / fps_data * 1000);
end
minTrackLength = sg_window + 2;

x_raw     = objs_link(1, :)';
y_raw     = objs_link(2, :)';
frame_raw = objs_link(5, :)';
id_raw    = objs_link(6, :)';

track_ids  = unique(id_raw)';   % row vector — must match orientation of tumble_freq etc.
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
firstsec_disp     = NaN(1, Nbacteria);   % first-second displacement (µm); 0 if track < 0.5s
total_migration   = NaN(1, Nbacteria);   % total displacement start→end (µm)
first_sec_speed   = NaN(1, Nbacteria);   % average speed over first second (µm/s)
avg_vel_total     = NaN(1, Nbacteria);   % mean speed over all frames (run + tumble) (µm/s)

%% =========================================================
%  Bacterium metrics (translational velocity and angular velocity)
%  =========================================================

for b = 1 : Nbacteria

    % ---- Extract valid (non-NaN) rows for this bacterium ----
    x_col = xmat(:, b);
    y_col = ymat(:, b);
    valid = ~isnan(x_col) & ~isnan(y_col);

    % ---- First-second displacement, first-second speed, total migration ----
    % All computed before any continue so every bacterium gets a value.
    % Uses NaN-safe moving median (window=3) to smooth centroid jitter,
    % matching migration_v3.m's x_use/y_use approach.
    nFrames1s_mig   = max(1, round(numFrames / secondVideo));
    nFrames_half_s  = round(0.5 * numFrames / secondVideo);  % frames in 0.5 s
    x_sm = movmedian_nan_col(x_col, 3);
    y_sm = movmedian_nan_col(y_col, 3);
    valid_idx_mig = find(~isnan(x_sm) & ~isnan(y_sm));

    if numel(valid_idx_mig) >= 2
        idx_start_mig  = valid_idx_mig(1);
        idx_end_target = idx_start_mig + nFrames1s_mig;
        idx_within_mig = valid_idx_mig(valid_idx_mig <= idx_end_target);
        n_within       = numel(idx_within_mig);
        elapsed_frames = idx_within_mig(end) - idx_within_mig(1);  % frame span

        % First-second displacement: only if track spans >= 0.5 s
        if n_within >= 2 && elapsed_frames >= nFrames_half_s
            xm = x_sm(idx_within_mig(end)) - x_sm(idx_within_mig(1));
            ym = y_sm(idx_within_mig(end)) - y_sm(idx_within_mig(1));
            if convert_to_um, xm = xm * px_to_um; ym = ym * px_to_um; end
            firstsec_disp(b) = sqrt(xm^2 + ym^2);

            % First-second speed: total path length / elapsed time
            xs_win = x_sm(idx_within_mig);
            ys_win = y_sm(idx_within_mig);
            if convert_to_um
                xs_win = xs_win * px_to_um;
                ys_win = ys_win * px_to_um;
            end
            path_len  = sum(sqrt(diff(xs_win).^2 + diff(ys_win).^2));
            elapsed_s = elapsed_frames * time_int;
            if elapsed_s > 0
                first_sec_speed(b) = path_len / elapsed_s;
            end
        end
        % tracks < 0.5 s: firstsec_disp stays NaN (leave as 0 per spec)
        if isnan(firstsec_disp(b)) && numel(valid_idx_mig) >= 2
            firstsec_disp(b) = 0;
        end

        % Total migration: straight-line displacement across entire track
        x_all = x_sm(valid_idx_mig);
        y_all = y_sm(valid_idx_mig);
        if convert_to_um
            x_all = x_all * px_to_um;
            y_all = y_all * px_to_um;
        end
        dx_tot = x_all(end) - x_all(1);
        dy_tot = y_all(end) - y_all(1);
        total_migration(b) = sqrt(dx_tot^2 + dy_tot^2);
    end

    if sum(valid) < minTrackLength
        continue
    end

    x_v = x_col(valid);
    y_v = y_col(valid);

    % ---- Pixel-to-micron conversion (applied before velocity — cascades to all outputs) ----
    if convert_to_um
        x_v = x_v * px_to_um;
        y_v = y_v * px_to_um;
    end

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

    % Total linear speed: mean of sg_v over ALL frames (run + tumble + no-tumble).
    % Computed here so it is populated even for bacteria that exit tumble
    % detection early (0 tumbles, angular confirmation failure, etc.).
    avg_vel_total(b) = mean(sg_v, 'omitnan');

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

    % NOTE: delta_v and v_min_used are appended in the SAME loop iteration
    % as v_tumble and t_tumble, so all four arrays are always the same length.
    % This is required for the angular confirmation filter below to work correctly.

    % ---- Case 1: equal counts, min comes first ----
    if n_min == n_max && min_idx(1) < max_idx(1)
        v_min_cand = v_min_all(2:end);
        for i = 1 : n_min-1
            dv = max(v_max_all(i+1) - v_min_cand(i), v_max_all(i) - v_min_cand(i));
            if dv / v_min_cand(i) > depth_speed_ratio
                v_tumble(end+1)   = v_min_cand(i);
                t_tumble(end+1)   = min_idx(i+1);
                delta_v(end+1)    = dv;
                v_min_used(end+1) = v_min_cand(i);
            end
        end

    % ---- Case 2: equal counts, max comes first ----
    elseif n_min == n_max && min_idx(1) > max_idx(1)
        v_min_cand = v_min_all(1:end-1);
        for i = 1 : n_min-1
            dv = max(v_max_all(i+1) - v_min_cand(i), v_max_all(i) - v_min_cand(i));
            if dv / v_min_cand(i) > depth_speed_ratio
                v_tumble(end+1)   = v_min_cand(i);
                t_tumble(end+1)   = min_idx(i);
                delta_v(end+1)    = dv;
                v_min_used(end+1) = v_min_cand(i);
            end
        end

    % ---- Case 3: more minima than maxima ----
    elseif n_min > n_max
        v_min_cand = v_min_all(2:end-1);
        for i = 1 : n_min-2
            dv = max(v_max_all(i+1) - v_min_cand(i), v_max_all(i) - v_min_cand(i));
            if dv / v_min_cand(i) > depth_speed_ratio
                v_tumble(end+1)   = v_min_cand(i);
                t_tumble(end+1)   = min_idx(i+1);
                delta_v(end+1)    = dv;
                v_min_used(end+1) = v_min_cand(i);
            end
        end

    % ---- Case 4: more maxima than minima ----
    elseif n_min < n_max
        v_min_cand = v_min_all;
        for i = 1 : n_min
            dv = max(v_max_all(i+1) - v_min_cand(i), v_max_all(i) - v_min_cand(i));
            if abs(dv / v_min_cand(i)) > depth_speed_ratio
                v_tumble(end+1)   = v_min_cand(i);
                t_tumble(end+1)   = min_idx(i);
                delta_v(end+1)    = dv;
                v_min_used(end+1) = v_min_cand(i);
            end
        end
    end

    n_tumble_detected = length(v_tumble);

    if n_tumble_detected == 0
        num_tumbles(b) = 0;
        continue
    end

    % ----------------------------------------------------------------
    %  ANGULAR SPEED CONFIRMATION  (Paper Appendix B, Criterion 2)
    %  Requires angular speed at the tumble frame (±1) to exceed the
    %  track-mean angular speed by w_confirm_factor (paper value: 0.20).
    %  All four arrays — v_tumble, t_tumble, delta_v, v_min_used — are
    %  the same length and are trimmed together with the same keep mask.
    % ----------------------------------------------------------------
    if use_angular_confirmation
        w_baseline  = mean(sg_w, 'omitnan');
        w_threshold = w_baseline * (1 + w_confirm_factor);
        keep = false(1, n_tumble_detected);
        for i = 1 : n_tumble_detected
            check = t_tumble(i) + (-1:1);
            check = check(check >= 1 & check <= N_sg);
            if any(sg_w(check) >= w_threshold)
                keep(i) = true;
            end
        end
        v_tumble          = v_tumble(keep);
        t_tumble          = t_tumble(keep);
        delta_v           = delta_v(keep);
        v_min_used        = v_min_used(keep);
        n_tumble_detected = length(v_tumble);
    end

    if n_tumble_detected == 0
        num_tumbles(b) = 0;
        continue
    end

    % ----------------------------------------------------------------
    %  MINIMUM SPEED GATE
    %  Reject tumble candidates where the speed at the dip is already
    %  near zero. These arise when the bacterium stops moving (wall,
    %  focal-plane exit) and only angular noise remains. A real tumble
    %  must start from an active run, so v_min must exceed the threshold.
    % ----------------------------------------------------------------
    if min_tumble_speed > 0
        speed_keep        = v_min_used > min_tumble_speed;
        v_tumble          = v_tumble(speed_keep);
        t_tumble          = t_tumble(speed_keep);
        delta_v           = delta_v(speed_keep);
        v_min_used        = v_min_used(speed_keep);
        n_tumble_detected = length(v_tumble);
    end

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
        num_tumbles(b)    = n_tumble_detected;
        tumble_x_mat{b}   = [];
        tumble_y_mat{b}   = [];
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
    elseif length(run_time) == 1
        % Only one tumble: keep the single interior run segment as-is.
        % (The leading/trailing boundary trimming doesn't apply here.)
    else
        run_time = [];
    end

    if isempty(run_time)
        % run_time was stripped because there was only 1 tumble and the
        % trimming logic removed both the first and last boundary.
        % We can still compute tumble-level metrics — fall through to the
        % summary block below rather than skipping with continue.
        % run_time stays empty; avg_run_t will remain NaN for this track
        % (the 0-tumble fill in the export section handles it correctly).
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
            if isempty(run_time)
                rt_prev = 0;
            else
                rt_prev = run_time(end);
            end
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
%  Export results to Excel
%  =========================================================

outFile = 'CellTrackAnalysis_results.xlsx';

% Delete existing file 
if exist(outFile, 'file')
    delete(outFile);
end

% ---- TrackDuration_s: ensure every bacterium with any valid frames gets a value ----
% track_duration_s is set inside the bacterium loop only when a track reaches
% the summary-metrics block.  Bacteria that exit early (short track, no tumbles
% found after filtering, etc.) are still NaN here.  Fill those gaps now from
% the raw valid-frame count so every row has a duration.
for b = 1 : Nbacteria
    if isnan(track_duration_s(b))
        valid_b = ~isnan(xmat(:, b)) & ~isnan(ymat(:, b));
        n_valid = sum(valid_b);
        if n_valid >= 1
            track_duration_s(b) = n_valid * time_int;
        end
    end
end

% ---- TumbleFrequency_per_s: TumbleCount / TrackDuration_s ----
tumble_freq = nan(1, Nbacteria);
valid_dur   = ~isnan(track_duration_s) & track_duration_s > 0 & ~isnan(num_tumbles);
tumble_freq(valid_dur) = num_tumbles(valid_dur) ./ track_duration_s(valid_dur);

% ---- MeanRunDuration_s: fill NaN cases ----
% 0 tumbles: whole track is a single run → use TrackDuration_s
% 1 tumble:  run time can't be measured from boundaries → TrackDuration_s - MeanTumbleDuration_s
mean_run_dur = avg_run_t(:)';
zero_tumble  = (~isnan(num_tumbles)) & (num_tumbles == 0);
one_tumble   = (~isnan(num_tumbles)) & (num_tumbles == 1);
mean_run_dur(zero_tumble) = track_duration_s(zero_tumble);
mean_run_dur(one_tumble & isnan(mean_run_dur)) = ...
    track_duration_s(one_tumble & isnan(mean_run_dur)) - avg_tumble_t(one_tumble & isnan(mean_run_dur));

% ---- MeanRunSpeed_um_per_s: if 0 tumbles use the whole-track average speed ----
% avg_vel is only computed over labelled run-frames; for 0-tumble bacteria
% every frame is a run frame, but avg_vel may still be NaN if the track
% exited the summary block early.  Compute a fallback from sg_v directly.
for b = 1 : Nbacteria
    if zero_tumble(b) && isnan(avg_vel(b)) && ~isempty(v_filt_mat{b})
        avg_vel(b) = mean(v_filt_mat{b}, 'omitnan');
    end
end

% ---- MeanTumbleDuration_s, MeanTumbleAngle_deg, MeanTumbleSpeed_um_per_s:
%      NaN when num_tumbles == 0 (already the default — no change needed).
%      When num_tumbles >= 1, the values were computed in the loop.         ----

% ---- Shared filtered vectors ----
% Defined once here and reused identically for T_PerTrack, Metadata,
% printout, and histograms — guaranteeing all outputs are consistent.
nt_plot   = num_tumbles(~isnan(num_tumbles));
tf_plot   = tumble_freq(~isnan(tumble_freq));
att_plot  = avg_tumble_t(~isnan(avg_tumble_t));
art_plot  = mean_run_dur(~isnan(mean_run_dur));
aa_plot   = avg_angle(~isnan(avg_angle));
vt_plot   = avg_v_tumble(~isnan(avg_v_tumble));
vr_plot   = avg_vel(~isnan(avg_vel));
vall_plot = avg_vel_total(~isnan(avg_vel_total));
wr_plot   = avg_w_run(~isnan(avg_w_run));
wt_plot   = avg_w_tumble(~isnan(avg_w_tumble));
dur_plot  = track_duration_s(~isnan(track_duration_s));
fsd_plot  = firstsec_disp(firstsec_disp > 0 & ~isnan(firstsec_disp));
tmig_plot = total_migration(~isnan(total_migration));
spd_plot  = first_sec_speed(~isnan(first_sec_speed));

% ---- Build per-track table ----
T_PerTrack = table( ...
    (1:Nbacteria)', ...
    num_tumbles(:), ...
    tumble_freq(:), ...
    avg_tumble_t(:), ...
    mean_run_dur(:), ...
    avg_angle(:), ...
    avg_v_tumble(:), ...
    avg_vel(:), ...
    avg_vel_total(:), ...
    total_migration(:), ...
    track_duration_s(:), ...
    firstsec_disp(:), ...
    first_sec_speed(:), ...
    'VariableNames', { ...
        'TrackID', ...
        'TumbleCount', ...
        'Tumbles_per_s', ...
        'MeanTumbleDuration_s', ...
        'MeanRunDuration_s', ...
        'MeanTumbleAngle_deg', ...
        'MeanTumbleSpeed_um_per_s', ...
        'MeanRunSpeed_um_per_s', ...
        'TotalLinearSpeed_um_per_s', ...
        'Migration_um', ...
        'TrackDuration_s', ...
        'FirstSecondDisplacement_um', ...
        'FirstSecondSpeed_um_per_s' ...
    });

% ---- Build Metadata summary-statistics table ----
% Uses the same _plot vectors defined above → guaranteed consistent with
% printout and histograms.
vecs = {nt_plot, tf_plot, att_plot, art_plot, aa_plot, vt_plot, vr_plot, vall_plot, tmig_plot, dur_plot, fsd_plot, spd_plot};

stat_fns  = {@mean, @median, @std, @var};
stat_names = {'Mean', 'Median', 'StdDev', 'Variance', 'SEM'};

Ncols = 12;
meta_vals = zeros(5, Ncols);
for s = 1 : 4
    for c = 1 : Ncols
        if isempty(vecs{c})
            meta_vals(s, c) = NaN;
        else
            meta_vals(s, c) = stat_fns{s}(vecs{c});
        end
    end
end
% SEM = std / sqrt(N)
for c = 1 : Ncols
    if isempty(vecs{c}) || numel(vecs{c}) < 2
        meta_vals(5, c) = NaN;
    else
        meta_vals(5, c) = std(vecs{c}) / sqrt(numel(vecs{c}));
    end
end

T_Meta = array2table(meta_vals, ...
    'RowNames', stat_names, ...
    'VariableNames', { ...
        'NumberOfTumbles', ...
        'TumblesPerSecond_1_per_s', ...
        'TumbleDuration_s', ...
        'RunDuration_s', ...
        'TumbleAngle_deg', ...
        'TumbleLinearSpeed_um_per_s', ...
        'RunLinearSpeed_um_per_s', ...
        'TotalLinearSpeed_um_per_s', ...
        'Migration_um', ...
        'TrackDuration_s', ...
        'FirstSecondDisplacement_um', ...
        'FirstSecondSpeed_um_per_s' ...
    });

% ---- Write to Excel ----
writetable(T_PerTrack, outFile, 'Sheet', 'PerTrack_Metrics', 'WriteMode', 'overwrite');
writetable(T_Meta,     outFile, 'Sheet', 'Metadata',         'WriteMode', 'overwrite', 'WriteRowNames', true);

fprintf('Excel saved → %s\n', outFile);

%% =========================================================
%  Histograms
%  =========================================================

dodgerblue = [0.118, 0.565, 1.000];
avgline    = [0.85, 0.10, 0.10];

% (shared _plot vectors already defined above, before Metadata)

% ---- Population-level summary statistics ----
pop_mean_nt   = mean(nt_plot);    pop_med_nt   = median(nt_plot);
pop_mean_att  = mean(att_plot);   pop_med_att  = median(att_plot);
pop_mean_art  = mean(art_plot);   pop_med_art  = median(art_plot);
pop_mean_aa   = mean(aa_plot);    pop_med_aa   = median(aa_plot);
pop_mean_dur  = mean(dur_plot);   pop_med_dur  = median(dur_plot);
pop_mean_vr   = mean(vr_plot);    pop_med_vr   = median(vr_plot);
pop_mean_vall = mean(vall_plot);  pop_med_vall = median(vall_plot);
pop_mean_vt   = mean(vt_plot);    pop_med_vt   = median(vt_plot);
pop_mean_wr   = mean(wr_plot);    pop_med_wr   = median(wr_plot);
pop_mean_wt   = mean(wt_plot);    pop_med_wt   = median(wt_plot);
pop_mean_tf   = mean(tf_plot);    pop_med_tf   = median(tf_plot);
pop_mean_tmig = mean(tmig_plot);  pop_med_tmig = median(tmig_plot);

% ---- Print summary table ----
fprintf('\n');
fprintf('=====================================================================\n');
fprintf('          Population Summary  (N = %d bacteria analysed)\n', length(nt_plot));
fprintf('=====================================================================\n');
fprintf('%-35s   %10s   %10s\n', 'Metric', 'Mean', 'Median');
fprintf('---------------------------------------------------------------------\n');
fprintf('%-35s   %10.2f   %10.2f\n', 'Number of tumbles',          pop_mean_nt,   pop_med_nt);
fprintf('%-35s   %10.2f   %10.2f\n', 'Tumbles per second (1/s)',   pop_mean_tf,   pop_med_tf);
fprintf('%-35s   %10.2f   %10.2f\n', 'Tumble duration (s)',        pop_mean_att,  pop_med_att);
fprintf('%-35s   %10.2f   %10.2f\n', 'Run duration (s)',           pop_mean_art,  pop_med_art);
fprintf('%-35s   %10.2f   %10.2f\n', 'Tumble angle (deg)',         pop_mean_aa,   pop_med_aa);
fprintf('%-35s   %10.2f   %10.2f\n', 'Tumble linear speed (µm/s)', pop_mean_vt,   pop_med_vt);
fprintf('%-35s   %10.2f   %10.2f\n', 'Run linear speed (µm/s)',    pop_mean_vr,   pop_med_vr);
fprintf('%-35s   %10.2f   %10.2f\n', 'Total linear speed (µm/s)',  pop_mean_vall, pop_med_vall);
fprintf('%-35s   %10.2f   %10.2f\n', 'Migration (µm)',             pop_mean_tmig, pop_med_tmig);
fprintf('%-35s   %10.2f   %10.2f\n', 'Track duration (s)',         pop_mean_dur,  pop_med_dur);
fprintf('=====================================================================\n\n');

if show_plots

% ---- Histogram panels (2 rows × 4 cols = 8 panels) ----
figure('Units','inches', 'Position',[1 1 22 9], 'Color','white');

subplot(2,4,1); smart_histogram(nt_plot,   pop_mean_nt,   'Number of Tumbles',     'Tumble Count',        dodgerblue, avgline, true);
subplot(2,4,2); smart_histogram(att_plot,  pop_mean_att,  'Tumble Duration (s)',    'Tumble Duration',     dodgerblue, avgline, false);
subplot(2,4,3); smart_histogram(art_plot,  pop_mean_art,  'Run Duration (s)',       'Run Duration',        dodgerblue, avgline, false);
subplot(2,4,4); smart_histogram(aa_plot,   pop_mean_aa,   'Tumble Angle (deg)',     'Tumble Angle',        dodgerblue, avgline, false);
subplot(2,4,5); smart_histogram(vt_plot,   pop_mean_vt,   'Tumble Speed (\mum/s)', 'Tumble Linear Speed', dodgerblue, avgline, false);
subplot(2,4,6); smart_histogram(vr_plot,   pop_mean_vr,   'Run Speed (\mum/s)',    'Run Linear Speed',    dodgerblue, avgline, false);
subplot(2,4,7); smart_histogram(tmig_plot, pop_mean_tmig, 'Migration (\mum)',      'Migration',           dodgerblue, avgline, false);
subplot(2,4,8); smart_histogram(dur_plot,  pop_mean_dur,  'Track Duration (s)',    'Track Duration',      dodgerblue, avgline, false);

% ---- Tumble frequency scatter plot ----
scatter_mask = ~isnan(tumble_freq) & ~isnan(track_duration_s);
x_sc    = tumble_freq(scatter_mask);
y_sc    = track_duration_s(scatter_mask);
ids_sc  = find(scatter_mask);

figure('Units','inches', 'Position',[1 1 7 6], 'Color','white');
scatter(x_sc, y_sc, 40, dodgerblue, 'filled', 'MarkerFaceAlpha', 0.6);
hold on;
xline(pop_mean_tf, '--', 'Color', avgline, 'LineWidth', 1.5, ...
    'Label', sprintf('Mean=%.4f', pop_mean_tf), ...
    'LabelOrientation', 'horizontal', 'FontSize', 10);

if show_scatter_labels
    x_pad_lbl = (max(x_sc) - min(x_sc)) * 0.012;
    for i = 1 : numel(x_sc)
        text(x_sc(i) + x_pad_lbl, y_sc(i), num2str(ids_sc(i)), ...
            'FontSize', 8, 'Color', [0.2 0.2 0.2], ...
            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
    end
end

hold off;
xlabel('Tumbles per Second (1/s)', 'FontSize', 14);
ylabel('Track Duration (s)',       'FontSize', 14);
title('Tumble Frequency vs Track Duration', 'FontSize', 14);
set(gca, 'FontSize', 12, 'TickDir', 'in', 'Box', 'on');

x_pad = (max(x_sc) - min(x_sc)) * 0.05;
y_pad = (max(y_sc) - min(y_sc)) * 0.05;
x_left_pad = max(x_pad, max(x_sc) * 0.08);
xlim([min(x_sc) - x_left_pad,  max(x_sc) + x_pad]);
ylim([max(0, min(y_sc) - y_pad),  max(y_sc) + y_pad]);

end   % end show_plots


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

function y = movmedian_nan_col(x, w)
% NaN-safe moving median over a column vector, window width w.
    n  = numel(x);
    y  = nan(size(x));
    hw = floor(w / 2);
    for i = 1 : n
        a  = max(1, i - hw);
        b  = min(n, i + hw);
        xi = x(a:b);
        xi = xi(~isnan(xi));
        if ~isempty(xi)
            y(i) = median(xi);
        end
    end
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