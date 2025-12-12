function [processparam, graddir]=setProcessparam(processopt, processparam, ...
                                        bpfiltsize, nsize, graddiroptval, ...
                                        grobjsize, defaultgrthresh)

% Set processparam, graddiroptval, graddir based on process option
graddir = 0;
switch processopt
    case 'spatialfilter'
        % Spatial filtering
        processparam = [bpfiltsize nsize];
    case 'gradientvote'
        % Gradient voting
         switch graddiroptval
            case 1
                % both directions
                graddir = 0;
            case 2
                % positive gradients
                graddir = 1;
            case 3
                % negative gradients
                graddir = -1;
            otherwise
                errordlg('Error -- bad graddiroptval')
        end
        processparam = [grobjsize graddir defaultgrthresh];
    case 'none'
        % none, but set region size as first (only) element of processparam
        processparam = bpfiltsize;
end