% cullTracks_f.m
% Modified by Clara Shin Nov 21, 2025

function objs_out = cullTracks_f(objs, cullString, cullThrash)
    
    switch cullString
        case {'StdDev'}
            % Calculate standard deviation per frame (px)
            track_property_function = @get_track_StdDev; % handle to function
            
        case {'TrackLength'}
            % Calculate length of each track (number of frames)
            track_property_function = @get_trackLength; % handle to function
    end
    
    %% Evaluate track properties -- call various functions
    [track_IDs, ~, iUtrk] = unique(objs(6,:)); % Unique Track IDs, and indexing
    track_properties = track_property_function(objs, track_IDs, iUtrk);
    
    %% Cull, based on parameters
    objs_out = cullTracks(objs, track_properties, cullThrash);
    % [track_IDs, ~, iUtrk] = unique(objs_out(6,:)); % Unique Track IDs, and indexing
    % track_properties = track_property_function(objs_out, track_IDs, iUtrk);
    % track_IDs = unique(objs_out(6,:));
    %%
    
    %% Functions for evaluating properties
    
    function trk_StdDevPerSqrtFrame = get_track_StdDev(objs, utrk, iUtrk)
    % Calculate standard deviation / Number of frames for each track
    % Do this by incrementing sum(x) and sum(x^2) arrays, using this to
    % calculate variance. This is *far* faster (~100x) than extracting the
    % columns corresponding to each track and calculating the variance of each
    % set of positions. See "deleted_from_trackedit.m" file -- Mar 11, 2020
        Nframes = zeros(1,length(utrk));
        sum_x = zeros(1,length(utrk));
        sum_x2 = zeros(1,length(utrk));
        sum_y = zeros(1,length(utrk));
        sum_y2 = zeros(1,length(utrk));
        progtitle = sprintf('cullTracks_function: Calculating variance...  ');
        progbar = waitbar(0, progtitle);  % will display progress
        for j=1:size(objs,2)
            % loop through each column
            % Increment sum(x), sum(x2), Nframes for the track this column is
            % part of.
            thisTrackNo = iUtrk(j); % the "unique" track index corresponding to this object
            % Increment positions; don't keep track of NaNs
            sum_x(thisTrackNo) = sum_x(thisTrackNo) + objs(1,j);
            sum_x2(thisTrackNo) = sum_x2(thisTrackNo) + (objs(1,j))^2;
            sum_y(thisTrackNo) = sum_y(thisTrackNo) + objs(2,j);
            sum_y2(thisTrackNo) = sum_y2(thisTrackNo) + (objs(2,j))^2;
            Nframes(thisTrackNo) = Nframes(thisTrackNo) + 1;
            if (mod(j,round(size(objs,2)/25))==0)
               waitbar(j/size(objs,2), progbar, strcat(progtitle, sprintf('Object %d of %d', j, size(objs,2))));
            end
        end
        close(progbar)
        % Calculate variance (Normalize by N rather than N-1, just to be
        % concise.)
        vx = sum_x2./Nframes - (sum_x./Nframes).^2;
        vy = sum_y2./Nframes - (sum_y./Nframes).^2;
        trk_StdDevPerSqrtFrame = sqrt((vx+vy)./Nframes); 
    end
    
    function trk_Length = get_trackLength(objs, utrk, iUtrk)
        % Determine the length of each track. 
        % Avoid "Nframes = sum(objs(6,:)==j);" -- very slow for large arrays
        trk_Length = zeros(1,length(utrk));
        for j=1:size(objs,2)
            thisTrackNo = iUtrk(j); % the "unique" track index corresponding to this object
            trk_Length(thisTrackNo) = trk_Length(thisTrackNo) + 1;
        end
    end
    
    function objs_out = cullTracks(objs, trackProperties, params)
    % Remove tracks with parameter values outside the params range (min, max)
    % If max < min, assume that these are Angles, with a desired range to keep
    % that crosses 360 degrees.
        utrk = unique(objs(6,:));  % all the unique track ids
        objs_out = zeros(size(objs));  % the largest it could possibly be
    
        nc = 1;  % number of columns, for re-sizing the array
        progtitle = sprintf('Culling tracks...  ');
        progbar = waitbar(0, progtitle);  % will display progress
        for j=1:length(utrk)
            if params(2) >= params(1)
                % Usual [min max]
                keep_track_condition = (trackProperties(j) >= params(1)) && (trackProperties(j) <= params(2));
            else
                % max < min, so assume need to wrap angles around 360
                keep_track_condition = (trackProperties(j) <= params(2)) || (trackProperties(j) >= params(1));
            end
            if keep_track_condition
                % keep this track
                trtmp = objs(:, objs(6,:)==utrk(j));  % objects that are part of track j
                objs_out(:,nc:(nc+size(trtmp,2))-1) = trtmp; % keep these
                %k = k+1;
                nc = nc + size(trtmp,2);
            end
            if (mod(j,round(length(utrk)/25))==0)
                waitbar(j/length(utrk), progbar, strcat(progtitle, sprintf('track %d of %d', j, length(utrk))));
            end
        end        
        close(progbar)
        nc = nc-1;
        objs_out = objs_out(:,1:nc);  % "re-sizing" the array
    end

end