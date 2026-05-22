% === migration_v3.m ===

%% ================= USER SETTINGS =================
matFile   = 'WT1.mat';                  
outFile   = 'cheYH2O21.xlsx';                 % <-- NEW name recommended

totalDur  = 30;                         % seconds
um_per_px = 0.33;                       % µm/px

% --- XY smoothing (reduces centroid jitter) ---
xy_smooth_frames = 3;                   % odd: 1,3,5

% --- k-frame displacement/turn-rate ---
k_disp_frames = 5;                      % 3–7 recommended; larger = less jitter sensitivity

% --- Gate: require enough displacement over k frames to define heading ---
min_step_um = 1.0;                      % µm over k frames (start 0.8–1.5)

% --- Tumble definition (angle-based) ---
% With k=5 at ~90 fps: k*dt ≈ 0.056 s.
% A 90° reorientation over 0.056 s corresponds to ~1600 deg/s.
turnrate_thresh_deg_s = 900;           % start 1200–2500, tune with WT vs cheY
min_tumble_s = 0.03;                    % seconds

% --- Track QC ---
min_valid_tr_frames = 30;               % require at least this many valid turn-rate samples

% --- Debug prints ---
debug_print_firstN = 5;                 % print diagnostics for first N tracks (0 to disable)

%% ================= 0) DELETE OUTPUT FILE (prevents stale/cached sheets) =================
if exist(outFile,'file')
    delete(outFile);
end

%% ================= 1) LOAD + POSITIONS =================
load(matFile);                          % must contain objs_link
[xmat, ymat] = objs2pos(objs_link);     % [frames × tracks]

%% ================= 2) TIMING =================
frameNums  = unique(objs_link(5,:));
nFrames    = numel(frameNums);
fps        = nFrames / totalDur;
dt         = 1 / fps;

fprintf('Frames=%d | totalDur=%.3f s | fps=%.3f | dt=%.6f s\n', nFrames, totalDur, fps, dt);

%% ================= 3) ALIGN MATRICES =================
xmat = xmat(1:min(end,nFrames),:);
ymat = ymat(1:min(end,nFrames),:);

if size(xmat,1) < nFrames
    nT = size(xmat,2);
    xmat(end+1:nFrames,1:nT) = NaN;
    ymat(end+1:nFrames,1:nT) = NaN;
end

nTracks  = size(xmat,2);
time_vec = (0:nFrames-1)' * dt;

%% ================= 4) SMOOTH X,Y (NaN-safe moving median) =================
x_use = xmat;
y_use = ymat;

if xy_smooth_frames > 1
    for t = 1:nTracks
        x_use(:,t) = movmedian_nan(xmat(:,t), xy_smooth_frames);
        y_use(:,t) = movmedian_nan(ymat(:,t), xy_smooth_frames);
    end
end

%% ================= 5) OPTIONAL QC: FIRST-SECOND METRICS =================
nFrames1s = max(1, round(fps));

dist1s_um  = nan(nTracks,1);
speed1s_um = nan(nTracks,1);

for t = 1:nTracks
    xi = x_use(:,t);
    yi = y_use(:,t);

    valid = find(~isnan(xi) & ~isnan(yi));
    if numel(valid) < 2
        continue
    end

    idx_start = valid(1);
    idx_end_target = idx_start + nFrames1s;
    idx_within = valid(valid <= idx_end_target);

    if numel(idx_within) < 2
        continue
    end

    idx_end = idx_within(end);

    dx = xi(idx_end) - xi(idx_start);
    dy = yi(idx_end) - yi(idx_start);
    dist1s_um(t) = sqrt(dx^2 + dy^2) * um_per_px;

    dx_step = diff(xi(idx_within));
    dy_step = diff(yi(idx_within));
    path_len_px = sum(sqrt(dx_step.^2 + dy_step.^2));

    duration_s = (numel(idx_within)-1) * dt;
    speed1s_um(t) = (path_len_px * um_per_px) / duration_s;
end

TrackNames = arrayfun(@(k) sprintf('Track %d',k),(1:nTracks).','Uni',0);
T_FirstSecond = table(TrackNames, dist1s_um, speed1s_um, ...
    'VariableNames', {'Track','FirstSecondDisplacement_um','MeanSpeed_FirstSecond_um_per_s'});

%% ================= 6) TURN RATE (deg/s) FROM k-FRAME HEADING =================
k = max(1, round(k_disp_frames));

turn_rate_deg_s = nan(nFrames, nTracks);
valid_tr_mask   = false(nFrames, nTracks);

% per-track diagnostics
validTR_count = zeros(nTracks,1);
maxTR_deg_s   = nan(nTracks,1);

for t = 1:nTracks
    xi = x_use(:,t);
    yi = y_use(:,t);

    % k-frame displacement vectors assigned to frame i
    dxk = nan(nFrames,1);
    dyk = nan(nFrames,1);
    dxk(1+k:end) = xi(1+k:end) - xi(1:end-k);
    dyk(1+k:end) = yi(1+k:end) - yi(1:end-k);

    stepk_um = sqrt(dxk.^2 + dyk.^2) * um_per_px;

    % heading from k-frame vectors
    headingk = atan2(dyk, dxk);
    headingk = unwrap_nan_segments(headingk);

    % gate: undefined heading if too little displacement
    headingk(stepk_um < min_step_um) = NaN;

    % turn angle over k frames
    dh = nan(nFrames,1);
    dh(1+k:end) = headingk(1+k:end) - headingk(1:end-k);
    dh = atan2(sin(dh), cos(dh));  % wrap to [-pi, pi]

    tr = abs(rad2deg(dh)) / (k*dt); % deg/s

    turn_rate_deg_s(:,t) = tr;

    vmask = ~isnan(tr) & isfinite(tr);
    valid_tr_mask(:,t) = vmask;

    validTR_count(t) = nnz(vmask);
    if validTR_count(t) > 0
        maxTR_deg_s(t) = max(tr(vmask));
    end

    if debug_print_firstN > 0 && t <= debug_print_firstN
        fprintf('Track %d: validTR=%d, maxTR=%.1f deg/s\n', t, validTR_count(t), maxTR_deg_s(t));
    end
end

%% ================= 7) TUMBLE EVENTS FROM TURN RATE =================
min_frames = max(1, round(min_tumble_s * fps));

tumble_count       = nan(nTracks,1);
tumble_freq_per_s  = nan(nTracks,1);
track_duration_s   = nan(nTracks,1);
frac_time_highturn = nan(nTracks,1);

for t = 1:nTracks
    tr = turn_rate_deg_s(:,t);
    valid = valid_tr_mask(:,t);

    if nnz(valid) < min_valid_tr_frames
        tumble_count(t) = 0;
        tumble_freq_per_s(t) = NaN;
        track_duration_s(t) = 0;
        frac_time_highturn(t) = NaN;
        continue
    end

    isHighTurn = (tr >= turnrate_thresh_deg_s) & valid;

    d = diff([0; isHighTurn; 0]);
    startIdx = find(d == 1);
    endIdx   = find(d == -1) - 1;

    durations = endIdx - startIdx + 1;
    keep = durations >= min_frames;

    tumble_count(t) = sum(keep);

    track_duration_s(t) = nnz(valid) * dt;
    tumble_freq_per_s(t) = tumble_count(t) / track_duration_s(t);

    frac_time_highturn(t) = nnz(isHighTurn) / nnz(valid);
end

T_Tumble = table(TrackNames, tumble_count, tumble_freq_per_s, track_duration_s, ...
                 validTR_count, maxTR_deg_s, frac_time_highturn, ...
    'VariableNames', {'Track','TumbleCount','TumbleFrequency_per_s','TrackDuration_s', ...
                      'ValidTurnRateFrames','MaxTurnRate_deg_per_s','FracTimeHighTurn'});

%% ================= 8) EXPORT COORDINATES (QC) =================
varX = ['Time_s', arrayfun(@(k) sprintf('Track%d_X_px',k),1:nTracks,'Uni',0)];
varY = ['Time_s', arrayfun(@(k) sprintf('Track%d_Y_px',k),1:nTracks,'Uni',0)];
T_X = array2table([time_vec, xmat], 'VariableNames', varX);
T_Y = array2table([time_vec, ymat], 'VariableNames', varY);

%% ================= 9) METADATA =================
T_Meta = table( ...
    string(matFile), string(outFile), totalDur, fps, dt, um_per_px, ...
    xy_smooth_frames, k_disp_frames, min_step_um, ...
    turnrate_thresh_deg_s, min_tumble_s, min_frames, min_valid_tr_frames, ...
    'VariableNames', {'MatFile','OutFile','TotalDur_s','FPS','dt_s','um_per_px', ...
                      'xy_smooth_frames','k_disp_frames','min_step_um', ...
                      'turnrate_thresh_deg_s','min_tumble_s','min_frames','min_valid_tr_frames'});

%% ================= 10) WRITE EXCEL =================
writetable(T_X, outFile, 'Sheet','X_Coordinates_px',           'WriteMode','overwrite');
writetable(T_Y, outFile, 'Sheet','Y_Coordinates_px',           'WriteMode','overwrite');
writetable(T_FirstSecond, outFile, 'Sheet','FirstSecond_Metrics','WriteMode','overwrite');
writetable(T_Tumble, outFile, 'Sheet','Tumble_Metrics',         'WriteMode','overwrite');
writetable(T_Meta, outFile, 'Sheet','Metadata',                 'WriteMode','overwrite');

fprintf('✅ Done. Per-track tumbles/s: Tumble_Metrics -> TumbleFrequency_per_s\n');
fprintf('Params: k=%d, min_step_um=%.2f, thr=%.1f deg/s, min_event=%.3f s (%d frames)\n', ...
    k, min_step_um, turnrate_thresh_deg_s, min_tumble_s, min_frames);

%% ================= LOCAL FUNCTIONS =================
function y = movmedian_nan(x, w)
% NaN-safe moving median without toolboxes.
    n = numel(x);
    y = nan(size(x));
    if w <= 1
        y = x;
        return
    end
    hw = floor(w/2);
    for i = 1:n
        a = max(1, i-hw);
        b = min(n, i+hw);
        xi = x(a:b);
        xi = xi(~isnan(xi));
        if isempty(xi)
            y(i) = NaN;
        else
            y(i) = median(xi);
        end
    end
end

function h_out = unwrap_nan_segments(h_in)
% Unwrap within contiguous non-NaN segments only.
    h_out = h_in;
    isn = isnan(h_in);
    n = numel(h_in);
    i = 1;
    while i <= n
        while i <= n && isn(i)
            i = i + 1;
        end
        if i > n, break; end
        j = i;
        while j <= n && ~isn(j)
            j = j + 1;
        end
        seg = h_in(i:j-1);
        seg = unwrap(seg);
        h_out(i:j-1) = seg;
        i = j;
    end
end