% === process_tracking_data.m ===
% Reads C272A HOCl 3.mat, exports full-length coordinates,
% and computes first-second displacement per track.

%% 1) Load your data
load('TrackGUIoutput.mat');              % contains objs_link

%% 2) Convert to position matrices
[xmat, ymat] = objs2pos(objs_link);    % [frames×tracks]

%% 3) Build time vector (0 → 30 s over nFrames)
frameNums  = unique(objs_link(5,:));
nFrames    = numel(frameNums);
totalDur   = 30;                       % total seconds
dt         = totalDur / nFrames;
time_vec   = (0:nFrames-1)' * dt;      % [nFrames×1]

%% 4) Align xmat/ymat to exactly nFrames rows
xmat = xmat(1:min(end,nFrames),:);
ymat = ymat(1:min(end,nFrames),:);
if size(xmat,1)<nFrames
    nT = size(xmat,2);
    xmat(end+1:nFrames,1:nT) = NaN;
    ymat(end+1:nFrames,1:nT) = NaN;
end

%% 5) Export full-duration X_Coordinates
nTracks = size(xmat,2);
varX    = ['Frame', arrayfun(@(k) sprintf('Track%d_X',k),1:nTracks,'Uni',0)];
T_X     = array2table([time_vec, xmat], 'VariableNames', varX);

%% 6) Export full-duration Y_Coordinates
varY    = ['Frame', arrayfun(@(k) sprintf('Track%d_Y',k),1:nTracks,'Uni',0)];
T_Y     = array2table([time_vec, ymat], 'VariableNames', varY);

%% 7) Compute each track's displacement in the first 1 second
dist1s  = nan(nTracks,1);
for t = 1:nTracks
    xi = xmat(:,t); yi = ymat(:,t);
    valid = find(~isnan(xi) & ~isnan(yi));
    if numel(valid)<2, continue; end

    % start at first appearance
    idx_start = valid(1);
    % restrict to indices where time_vec<=1 s
    idx_within1 = valid(time_vec(valid)<=1);
    if numel(idx_within1)<2, continue; end
    idx_end = idx_within1(end);

    dx = xi(idx_end)-xi(idx_start);
    dy = yi(idx_end)-yi(idx_start);
    dist1s(t) = sqrt(dx^2 + dy^2);
end

TrackNames = arrayfun(@(k) sprintf('Track %d',k),(1:nTracks).','Uni',0);
T_D        = table(TrackNames, dist1s, 'VariableNames',{'Track','FirstSecondDisplacement'});

%% 8) Write to Excel
outFile = 'tracked_positions_firstSecond.xlsx';
writetable(T_X, outFile, 'Sheet','X_Coordinates',   'WriteMode','overwrite');
writetable(T_Y, outFile, 'Sheet','Y_Coordinates',   'WriteMode','overwrite');
writetable(T_D, outFile, 'Sheet','DistanceResults', 'WriteMode','overwrite');

fprintf('✅ Exported %d tracks: full coords + first-1s displacement to %s\n', nTracks, outFile);