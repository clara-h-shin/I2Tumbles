% Program:  objs2pos.m
%
% Summary:  Extract matrices of object positions froms objs matrix,
%           output of TrackingGUI
%
% Inputs:   objs = objs matrix, output from TrackingGUI
%
% Outputs:  xmat, ymat = matrices of positions, time down, track # to right
%
% Author:   Brandon Schlomann
%
% Date:     5/12/16 - first written
%

function [xmat,ymat] = objs2pos(objs)

% get ids
ids = sort(unique(objs(6,:)));

% number of uniqe objects
nobjs = length(ids);

% preallocate NaN arrays with the biggest possible size
biglength  = size(objs,2);

xmat = NaN(biglength,nobjs);
ymat = NaN(biglength,nobjs);

% use a counter to keep track of the longgest track.  Will
% use this to trim NaN padding.
longesttrack = 0;

% Loop through ids, assemble matrix, update longgestrack counter
for i = 1:nobjs
    xtmp = objs(1,objs(6,:)==ids(i));
    ytmp = objs(2,objs(6,:)==ids(i));
    xmat(1:length(xtmp),i) = xtmp';
    ymat(1:length(ytmp),i) = ytmp';
    
    longtmp = sum(~isnan(xmat(:,i)));
    if longtmp > longesttrack
        longesttrack = longtmp;
    end
end

% trim NaN pads
xmat = xmat(1:longesttrack,:);
ymat = ymat(1:longesttrack,:);

end

