%% StatTesting.m
% Performs Kruskal-Wallis statistical testing across multiple groups.
% Each group corresponds to one CellTrackAnalysis_results.xlsx file.
% < IMPORTANT: Please change the file names after running CellTrackAnalysis.m >
% Reads a single column from the PerTrack_Metrics sheet of each file.

% Last updated: 6/2/2026
% By Clara Shin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%% % %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MAIN PARAMETERS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

group_count=        3;              % number of groups (= number of Excel files to select)
                                    % must be a positive integer

excel_column_name=  'Tumbles_per_s';% column to read from the PerTrack_Metrics sheet
                                    % e.g. 'TumbleCount', 'Tumbles_per_s',
                                    %      'MeanRunDuration_s', 'TrackDuration_s'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% =========================================================
%  Validate parameters
%  =========================================================

if ~isnumeric(group_count) || group_count < 2 || floor(group_count) ~= group_count
    error('group_count must be an integer >= 2. Got: %g', group_count);
end
group_count = int32(group_count);

%% =========================================================
%  Select files
%  =========================================================

file_paths  = cell(1, group_count);
group_names = cell(1, group_count);

fprintf('Please select %d Excel file(s), one per group.\n', group_count);

for g = 1 : group_count
    [fname, fpath] = uigetfile('*.xlsx', ...
        sprintf('Select Excel file for Group %d of %d', g, group_count));

    if isequal(fname, 0)
        error('File selection cancelled for Group %d. Aborting.', g);
    end

    file_paths{g}  = fullfile(fpath, fname);
    group_names{g} = sprintf('Group%d', g);   % Group1, Group2, Group3, ...

    fprintf('  Group %d → %s\n', g, fname);
end

%% =========================================================
%  Read column from each file
%  =========================================================

groups = cell(1, group_count);   % groups{g} = column vector of values

for g = 1 : group_count
    try
        T = readtable(file_paths{g}, 'Sheet', 'PerTrack_Metrics');
    catch ME
        error('Could not read sheet "PerTrack_Metrics" from file:\n  %s\nError: %s', ...
            file_paths{g}, ME.message);
    end

    if ~ismember(excel_column_name, T.Properties.VariableNames)
        error('Column "%s" not found in PerTrack_Metrics sheet of:\n  %s\nAvailable columns: %s', ...
            excel_column_name, file_paths{g}, ...
            strjoin(T.Properties.VariableNames, ', '));
    end

    col_data = T.(excel_column_name);

    % Keep only finite numeric values (drop NaN / Inf)
    col_data = col_data(isfinite(col_data));

    if isempty(col_data)
        error('Column "%s" in Group %d contains no finite values after removing NaNs.', ...
            excel_column_name, g);
    end

    groups{g} = col_data(:);   % ensure column vector
    fprintf('  Group %d: %d valid values read.\n', g, numel(col_data));
end

%% =========================================================
%  Data matrix for testing
%  =========================================================

% Assign group1, group2, ... as individual workspace variables
for g = 1 : group_count
    assignin('base', group_names{g}, groups{g});
end

% kruskalwallis expects a matrix where each column is one group.
% Pad shorter columns with NaN so all columns have equal length.
max_n      = max(cellfun(@numel, groups));
dataMatrix = NaN(max_n, group_count);

for g = 1 : group_count
    n = numel(groups{g});
    dataMatrix(1:n, g) = groups{g};
end

fprintf('\nData matrix built: %d rows × %d groups.\n', max_n, group_count);

%% =========================================================
%  Kruskal-Wallis test
%  =========================================================

fprintf('\n--- Kruskal-Wallis Test ---\n');
fprintf('Column tested : %s\n', excel_column_name);
fprintf('Groups        : %s\n\n', strjoin(group_names, ', '));

[p, tbl, stats] = kruskalwallis(dataMatrix, group_names, 'on');

fprintf('p-value       : %.6f\n', p);

if p < 0.001
    sig_str = '*** (p < 0.001)';
elseif p < 0.01
    sig_str = '** (p < 0.01)';
elseif p < 0.05
    sig_str = '* (p < 0.05)';
else
    sig_str = 'n.s. (p >= 0.05)';
end

fprintf('Significance  : %s\n\n', sig_str);

%% =========================================================
%  Post-hoc: Dunn-Sidak multiple comparisons (if p < 0.05)
%  =========================================================

if p < 0.05 && group_count > 2
    fprintf('--- Post-hoc Multiple Comparisons (multcompare) ---\n');
    figure('Units','inches', 'Position',[1 1 8 5], 'Color','white');
    [results, means] = multcompare(stats, 'CriticalValueType', 'dunn-sidak');
    title(sprintf('Multiple Comparisons — %s', excel_column_name), ...
        'FontSize', 13, 'FontWeight', 'bold');
    fprintf('Pairwise comparison results: \n\n');
end

%% =========================================================
%  Summary statistics per group
%  =========================================================

fprintf('--- Summary Statistics: %s ---\n', excel_column_name);
fprintf('%-12s  %6s  %10s  %10s  %10s  %10s\n', ...
    'Group', 'N', 'Mean', 'Median', 'Std', 'IQR');
fprintf('%s\n', repmat('-', 1, 66));

for g = 1 : group_count
    d = groups{g};
    fprintf('%-12s  %6d  %10.4f  %10.4f  %10.4f  %10.4f\n', ...
        group_names{g}, numel(d), mean(d), median(d), std(d), iqr(d));
end
fprintf('\n');