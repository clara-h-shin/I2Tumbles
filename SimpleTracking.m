
% SimpleTracking - Version 2.0

% This script was written to utilize TrackingGUI without the GUI and
% pre-set settings that are frequently used in Baylink Lab. 

% Set the parameters before running. 
% Please refer to SimpleTracking_Instructions.pdf for detailed instructions.

% Baylink Lab
% www.baylink-lab.com

% Version 1.0 Written Fall 2025
% Version 2.0 Written Spring 2026

% Last updated: 5/22/2026
% By Clara Shin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MAIN PARAMETERS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bpfiltsize=         20;                     % bpfiltsize        [A]
hthreshopt=         2;                      % Threshold Option  [B] 2. Std. Thresh.
defaultthresh2=     4.0;                    % thr (>1)          [C]

stdevMinThreshold=  2;                      % Cull: Std. Min. Threshold
numFrames=          200;                    % Cull: # of Video Image Frames
secondVideo=        10;                     % Cull: Length of Video (sec)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Neighborhood Parameters
processopt=         'spatialfilter';        % Process Option
nsize=              bpfiltsize;             % nsize (Same as bpfiltsize)
lockobjsize=        true;                   % Lock

grobjsize=          7;                      % grobjsiz
graddiroptval=      1;         % dir. 1: Both dir, 2. Positive, 3. Negative
defaultgrthresh=    0.99;                   % grthresh

defaultthresh1=     0.9900;                 % thr (0-1)
defaultthresh3=     3;                      % N (>=1)

% Display Processed Images
showprocess=        false;                  % DispProcess
showthreshprocess=  false;                  % DispThreshProcess

% Tracking
ctrfitstr=          'radial';               % Localiz. method
orientstr=          'none';                 % Orient. method
sizestr=            'none';                 % Size method

husenhoodctrs=      false;                  % Use prev. ctrs.
h1pernhood=         false;                  % 1/nhood

linkstep=           1000;                   % Link MaxStep^2
linkmem=            0;                      % Link Memory

% Parameters for displayimage.m
hdispcircles=       true;                   % Display circles
hdispIDs=           false;                  % Display object IDs
hdisptracks=        false;                  % Display tracks

% Cull Options
cullOptionStrings=  {'StdDev', 'TrackLength'};
cullMaxThreshold=  Inf;                    % Cull: Std. Max. Threshold

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% BEGINNNING OF CODE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get Filenames and images info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

programdir = cd; % present directory.

if ~exist('im', 'var') || isempty(im)
    % images were not input; get file name info
    imagesloaded = false;
    %    Can be multipage TIFF
    [fbase, frmin, frmax, formatstr, FileName1, ~, PathName1, ext, ismultipage] = ...
        getnumfilelist;
    % If there's only one image, getnumfilelist returns empty arrays for
    %   frmin, frmax, formatstr.  Replace frmin and frmax with 1s.
    %   fbase is the complete filename
    if isempty(frmin)
        frmin = 1;
        frmax = 1;
    end
    if frmax < frmin
        % User may have selected the ending image as the start; reverse
        warndlg('"Start" and "End" image files appear to be reversed. Flipping these...')
        pause(1)
        temp_fr = frmax;
        frmax = frmin;
        frmin = temp_fr;
    end
else
    % images were input
    imagesloaded = true;
    frmin = 1;
    frmax = size(im,3);
    FileName1 = '[input 3D array]';
end
Nframes = frmax - frmin + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create all variables here

% Outputs
processparam = [];
thresh = defaultthresh1;
threshopt = [];

% Set processparam, graddir based on process option
% Call setProcessparam.m
[processparam, graddir]=setProcessparam(processopt, processparam, ...
                                        bpfiltsize, nsize, graddiroptval, ...
                                        grobjsize, defaultgrthresh);

objs = [];
objs_link = [];
objs_original = []; % to revert completely -- first objs matrix that contains all frames.
objs_link_original = []; % to revert completely -- first objs_link matrix
savedobjs = [];   % to revert if un-culling objects
savedobjs_link = []; % to revert if un-culling tracks
objs_dedrift = []; % de-drifted linked object (objs_link) array, to save (not display);
savedobjs_dedrift = []; % to revert if un-culling tracks

% Other
fitstr = {ctrfitstr; orientstr; sizestr};
lsqoptions = [];
colchannel = 0;
outfile = [];

if exist('optimset', 'file') && exist('lsqnonlin', 'file')
    % use try / catch, since this still seems to lead to problems:
    try
        lsqoptions = optimset('lsqnonlin');
    catch
        warning('Problem assigning lsqoptions in TrackingGUI_rp.m; leave empty')
        lsqoptions = []; % fine; see below
    end
else
    % Optimization toolbox probably doesn't exist. Send fo5_rp an empty
    % options variable, which will only cause problems if non-linear least
    % squares fitting is used.
    lsqoptions = [];
end

% Load one image, or use the first frame of the 'im' array
% In general, will use 'A' for the displayed image
if imagesloaded
    A = im(:,:,1);
    isColor = false;
else
    cd(PathName1) %Go to the directory of the images
    if ~ismultipage
        A=imread(FileName1);    % load the first file image
    else
        % multipage TIFF
        A = imread(FileName1, 1);
    end
    isColor = (ndims(A)==3);  % TRUE if the image is color (3 layers)
    if isColor
        prompt = {'If color: which channel to use?  (1, 2, 3): '};
        dlg_title = 'Color option'; num_lines= 1;
        % guess at which channel to use as default, based on total brightness
        b = zeros(1,3);
        for ch = 1:3
            b(ch) = sum(sum(A(:,:,ch)));
        end
        [~, ic] = max(b);
        def     = {num2str(ic)};  % default values
        answer  = inputdlg(prompt,dlg_title,num_lines,def);
        colchannel = round(str2double(answer(1)));
        A = A(:,:,colchannel);  % select the color channel to examine
        pause(0.5);  % seems to be necessary to keep MATLAB from crashing!
    end
end

% Determine the bit depth of the image, for adjustable display range
if isa(A, 'uint8') % is A a 'uint8' variable?
    bitDepth = 8;
elseif isa(A, 'uint16') % is z a 'uint16' variable?
    bitDepth = 16;
else
    bitDepth = [];
end

displayRange = []; % 0-1, Define here to use later. Max caxis will be displayRange*2^bitDepth-1)
currframe = frmin;  % # of the current (primary) frame

% Show filtered, thresholded image
processedA = [];  % filtered image
threshprocessA = [];  % thresholded, filtered image

cd(programdir); % Go back to the directory from which the GUI was called

% Handles to display images
im1 = [];  % handle to primary image display

% Status variables
trackthisdone = false(Nframes,1);  % true for frames that have been segmented
islinkdone = false;  % is linkage of objects into tracks done?
isdeDriftDone = false; % is there a de-drifted object array?
colormap('gray');

% make the filtering and neighborhood object sizes the same
if lockobjsize
    nsize = round(bpfiltsize);
    % update the filtered images, if these are being displayed
    [processedA, threshprocessA]=updateprocessandthresh(processedA, threshprocessA, ...
        showprocess, showthreshprocess);
end

thresh=whichthresh(thresh, hthreshopt, ...
    defaultthresh1, defaultthresh2, defaultthresh3);

[processedA, threshprocessA]=updateprocessandthresh(processedA, threshprocessA, ...
        showprocess, showthreshprocess);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Track This Frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trackthisdone = false;                  % Track this frame - Done
trackalldone = false;                   % Track all frames - Done
islinkdone = false;                     % Link Objs -> Tracks - Done

% Track a single frame (find objects)
tmpobj = fo5_rp(A, processopt, processparam, thresh, fitstr, ...
            h1pernhood, [], lsqoptions);
tmpobj(5,:) = currframe-frmin+1;
objs = [objs tmpobj];
trackthisdone(currframe-frmin+1) = true;  % note that tracking has been done
savedobjs = objs;  % saved, in case culling and un-culling are done.

% Make an image pop up then ask to proceed
% Call displayImage.m

[im1,objsthisframe]=displayImage(A, im1, ...
                        showthreshprocess, showprocess, ...
                        threshprocessA, processedA, displayRange, bitDepth, ...
                        trackthisdone, islinkdone, ...
                        currframe, frmin, objs_link, objs, ...
                        hdispcircles, hdispIDs, hdisptracks, ...
                        orientstr, sizestr);

% Ask if the image is highlighting correct objects:
askProceed = questdlg(strcat('Does the image highlight correct cell objects?'), 'show objs', 'Yes', 'No', 'Yes');  % last item is default
answerProceed = strcmpi(askProceed, 'yes');

if ~answerProceed
    % Proceed if 'Yes', otherwise close the program
    close;
    return;
end
close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Track All Frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

progtitle = sprintf('Tracking All Frames : ');
progbar = waitbar(0, progtitle);  % will display progress
if (Nframes > 1)
    progbar = waitbar(0, progbar, progtitle);  % will display progress
end
objs = [];
oldA = A;
oldcurrframe = currframe;

for j = frmin:frmax
    % Loop through all frames
    % similar code as in im2obj_rp.m
    if husenhoodctrs && (j>frmin)
        % use x, y positions from
        % previous frame for neighborhood centers of this frame
        nhoodctrs = [tmpobj(1,:)' tmpobj(2,:)'];  
    else
        nhoodctrs = [];
    end
    if imagesloaded
        % all images already loaded
        tmpobj = fo5_rp(im(:,:,j), processopt, processparam, ...
                        thresh, fitstr, get(h1pernhood, 'value'), ...
                        nhoodctrs, lsqoptions);
    else
        % read from file; % use the variablea "A" and "currframe"
        currframe = j;
        [A, processedA, threshprocessA]=loadImages(PathName1, ...
                        ismultipage, Nframes, A, fbase, ...
                        formatstr, currframe, ext, programdir, isColor, colchannel, ...
                        showprocess, showthreshprocess, ...
                        processedA, threshprocessA, imagesloaded);
        tmpobj = fo5_rp(A, processopt, processparam, ...
                        thresh, fitstr, h1pernhood, ...
                        nhoodctrs, lsqoptions);
    end
    if ~isempty(tmpobj)
        tmpobj(5,:) = j-frmin+1;
    end
    objs = [objs tmpobj];
    % show progress -- not called if just one frame
    if mod(j-frmin+1,10)==0
        waitbar((j-frmin+1)/Nframes, progbar, ...
            strcat(progtitle, sprintf(' Frame %d of %d', (j-frmin+1), Nframes)), 'Interpreter','none');
    end
    trackthisdone(j-frmin+1) = true;  % note that tracking has been done
end
A = oldA;
currframe = oldcurrframe;
if Nframes>1
    close(progbar)
end
savedobjs = objs;  % saved, in case culling and un-culling are done.

if sum(trackthisdone)==length(trackthisdone)
    % All frames have been tracked
    if isempty(objs_original)
        objs_original = objs; % This is the first objs matrix; save if we want to revert
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Link Objs -> Tracks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hdisptracks = true;                  % Display tracks

% Linking objects into trajectories, using nnlink_rp
if sum(trackthisdone)==length(trackthisdone)
    objs_link = nnlink_rp(objs, linkstep, linkmem, true);
    if isempty(objs_link_original)
        objs_link_original = objs_link; % the first linked object matrix, in case we want to revert
    end

    islinkdone = true;
    savedobjs_link = objs_link;  % saved, in case culling and un-culling are done.

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cull
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Perform Culling of either objects or tracks.
% Options (from Menu, hcullOpt):
% 3: Cull tracks from objs_link, based on criteria such as standard
%    deviation, track length, etc. 
%    Does not change objs, only objs_link
%    Also saves the dedrifted & culled matrix (objs_dedrift), which
%    can be output.
%    Uses cullTracks_function.m, and doesn't use dedrift_rp.m
%    since we don't need dedrift function.

% Cull Tracks (Option 3 above)
if islinkdone
    for i=1:length(cullOptionStrings)
        savedobjs = objs;
        savedobjs_link = objs_link; 
        objs_toCull = objs_link;

        switch cullOptionStrings{1,i}
            case 'StdDev'
                cullThrash = [stdevMinThreshold cullMaxThreshold];
            case 'TrackLength'
                frameRate = numFrames * 1.0 / secondVideo;
                cullThrash = [frameRate cullMaxThreshold];
        end
        objs_culled = cullTracks_f(objs_toCull, cullOptionStrings{1,i}, cullThrash);
        objs_link = objs_culled;
    end

    disp('Number of Tracks After Cull:');
    disp(length(unique(objs_link(6,:))));

    [im1,objsthisframe]=displayImage(A, im1, ...
                        showthreshprocess, showprocess, ...
                        threshprocessA, processedA, displayRange, bitDepth, ...
                        trackthisdone, islinkdone, ...
                        currframe, frmin, objs_link, objs, ...
                        hdispcircles, hdispIDs, hdisptracks, ...
                        orientstr, sizestr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save results in a MAT file
outfile = 'SimpleTrackingoutput.mat';
if ~isempty(outfile)
    save(outfile, 'Nframes', 'objs', 'objs_link', 'objs_dedrift',  ...
        'threshopt', 'thresh', 'processopt', 'processparam', ...
        'trackthisdone', 'islinkdone', 'isdeDriftDone');  % Save these variables
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kailie's MATLAB Coding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('SimpleTrackingoutput.mat');
uniqueFrames=unique(objs_link(5,:));
objs_in_frame=zeros(size(uniqueFrames));
for j=1:length(uniqueFrames)
    objs_in_frame(j)=sum(objs_link(5,:)==uniqueFrames(j));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% migration.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% === process_tracking_data.m ===
% Reads tracking data, exports coordinates,
% computes first-second displacement & speed,
% and instantaneous per-frame speeds (µm/s)

%% 1) Load data
load('SimpleTrackingoutput.mat');      % contains objs_link

%% 2) Convert to position matrices
[xmat, ymat] = objs2pos(objs_link);        % [frames × tracks]

%% 3) Video timing + calibration
frameNums  = unique(objs_link(5,:));
nFrames    = numel(frameNums);
totalDur   = 30;                           % seconds
fps        = nFrames / totalDur;
dt         = 1 / fps;

um_per_px  = 0.33;                         % spatial calibration

%% 4) Align matrices to nFrames
xmat = xmat(1:min(end,nFrames),:);
ymat = ymat(1:min(end,nFrames),:);

if size(xmat,1) < nFrames
    nT = size(xmat,2);
    xmat(end+1:nFrames,1:nT) = NaN;
    ymat(end+1:nFrames,1:nT) = NaN;
end

nTracks = size(xmat,2);

%% 5) Export full-duration coordinates
time_vec = (0:nFrames-1)' * dt;

varX = ['Time_s', arrayfun(@(k) sprintf('Track%d_X_px',k),1:nTracks,'Uni',0)];
varY = ['Time_s', arrayfun(@(k) sprintf('Track%d_Y_px',k),1:nTracks,'Uni',0)];

T_X = array2table([time_vec, xmat], 'VariableNames', varX);
T_Y = array2table([time_vec, ymat], 'VariableNames', varY);

%% 6) First-second displacement & mean speed (per track)
nFrames1s = round(fps);

dist1s_um  = nan(nTracks,1);
speed1s_um = nan(nTracks,1);

for t = 1:nTracks
    xi = xmat(:,t);
    yi = ymat(:,t);

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

    % Straight-line displacement
    dx = xi(idx_end) - xi(idx_start);
    dy = yi(idx_end) - yi(idx_start);
    dist1s_um(t) = sqrt(dx^2 + dy^2) * um_per_px;

    % Path length for mean speed
    dx_step = diff(xi(idx_within));
    dy_step = diff(yi(idx_within));
    path_len_px = sum(sqrt(dx_step.^2 + dy_step.^2));

    duration_s = (numel(idx_within)-1) * dt;
    speed1s_um(t) = (path_len_px * um_per_px) / duration_s;
end

%% 7) Instantaneous speed per frame (µm/s)
% One column per track, NaN where undefined

inst_speed_um_s = nan(nFrames, nTracks);

for t = 1:nTracks
    xi = xmat(:,t);
    yi = ymat(:,t);

    valid = find(~isnan(xi) & ~isnan(yi));
    if numel(valid) < 2
        continue
    end

    dx = diff(xi(valid));
    dy = diff(yi(valid));
    step_dist_um = sqrt(dx.^2 + dy.^2) * um_per_px;

    inst_speed = step_dist_um / dt;    % µm/s

    inst_speed_um_s(valid(2:end), t) = inst_speed;
end

%% 8) Build results tables
TrackNames = arrayfun(@(k) sprintf('Track %d',k),(1:nTracks).','Uni',0);

T_D = table(TrackNames, dist1s_um, speed1s_um, ...
    'VariableNames', ...
    {'Track','FirstSecondDisplacement_um','MeanSpeed_FirstSecond_um_per_s'});

varS = ['Time_s', arrayfun(@(k) sprintf('Track%d_Speed_um_per_s',k), ...
        1:nTracks,'Uni',0)];

T_S = array2table([time_vec, inst_speed_um_s], 'VariableNames', varS);

%% 9) Write to Excel
outFile = 'tracked_positions.xlsx';

writetable(T_X, outFile, 'Sheet','X_Coordinates_px',      'WriteMode','overwrite');
writetable(T_Y, outFile, 'Sheet','Y_Coordinates_px',      'WriteMode','overwrite');
writetable(T_D, outFile, 'Sheet','FirstSecond_Metrics',   'WriteMode','overwrite');
writetable(T_S, outFile, 'Sheet','Instantaneous_Speed',   'WriteMode','overwrite');

fprintf('✅ Exported %d tracks with calibrated displacement and speed to %s\n', ...
        nTracks, outFile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Necessary Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A, processedA, threshprocessA]=loadImages(PathName1, ...
                        ismultipage, Nframes, A, fbase, ...
                        formatstr, currframe, ext, programdir, isColor, colchannel, ...
                        showprocess, showthreshprocess, ...
                        processedA, threshprocessA, imagesloaded)

    if imagesloaded
        % all images are in the 'im' array
        A = im(:,:,currframe);
    else
        % load from file
        cd(PathName1)     %Go to the directory of the images
        % load the image into "A"
        % Don't do any calculations or updates
        % Note that all file names were determined previously
        if ~ismultipage
            if Nframes==1
                A = imread(fbase);
            else
                framestr = sprintf(formatstr, currframe);
                A  = imread(strcat(fbase, framestr, ext));
            end
        else
            % multipage TIFF
            A  = imread(strcat(fbase, ext), currframe);
        end
        cd(programdir) %go back to the original directory
    
        if isColor
            % select the color channel to examine
            A = A(:,:,colchannel);
        end
    end
    
    % if the boxes are checked, calculate the filtered and thresholded
    % images
    if showprocess || showthreshprocess
        processedA = calcprocessimg;
    end
    if showthreshprocess
        threshprocessA = calcthreshprocess;
    end
end


function [processedA, threshprocessA]=updateprocessandthresh(processedA, threshprocessA, ...
        showprocess, showthreshprocess)
    % calls functions to recalculate filtered and post-threshold
    % images, and display
    if showprocess || showthreshprocess
        processedA = calcprocessimg;
    end
    if showthreshprocess
        threshprocessA = calcthreshprocess;
    end
end

function thresh=whichthresh(thresh, hthreshopt, ...
    defaultthresh1, defaultthresh2, defaultthresh3)
        % determine which threshold value to use, based on thresholding
        % option
        switch hthreshopt
            case 1
                thresh = defaultthresh1;
            case 2
                thresh = -1.0*defaultthresh2;
                % fo5_rp.m interprets negative thresholds as "option 2"
                % inputs.
            case 3
                thresh = round(defaultthresh3);
        end
    end

function processedA = calcprocessimg
    % calculate the 'processed' image, for neighborhood finding
    % either spatially filtered, or gradient voting option
    % Note that this simply copies the calculation in fo5_rp.m -- 
    %   recalculated when fo5_rp.m is called.
    switch processopt
        case 'spatialfilter'
            if processparam(1)>0
                processedA = bpass(A,1,processparam(1));
            else
                processedA = A;  % if filter parameter = 0, don't filter.
            end
        case 'gradientvote'
            processedA = gradientvote(A, processparam(1), processparam(2), processparam(3));
        case 'none'
            processedA = A;
    end
end
    
function threshprocessA = calcthreshprocess
    % calculate the dilated image of local maxima that pass the
    % threshold.  
    % Calculate processed image even if previously calculated, in case
    % processing parameters have changed.
    processedA = calcprocessimg;  % uses bpfiltsize for filtering size
    % neighborhod size
    switch processopt
        case 'spatialfilter'
            nsize = processparam(2);
        case 'gradientvote'
            nsize = processparam(1);
        case 'none'
            nsize = processparam(2);
    end
    [y, x] = calcthreshpts(processedA, hthreshopt, thresh, nsize);
    threshprocessAmask = false(size(processedA));
    threshprocessAmask(sub2ind(size(processedA), round(y), round(x)))=true;
    threshprocessAmask = imdilate(threshprocessAmask, strel('disk', floor(nsize/2)));
    threshprocessA = processedA.*threshprocessAmask;
end