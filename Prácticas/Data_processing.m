% ***********************************************************************
%                   DATA PROCESSING OF DYNAMIC TESTING
% ***********************************************************************

% Collects all the numeric data from the beam's testing and places them
% into different structures to automate the post-processing

clc;clear;close all

%% BEAM'S PROPERTIES

beam.input = importdata('data\geomdata.txt',';',2);
beam.m = beam.input.data(:,1);
beam.L = beam.input.data(:,2);
beam.b = beam.input.data(:,3);
beam.t = beam.input.data(:,4);
beam.Ixx = beam.b.*beam.t.^3./12;                      
beam.rho = beam.m./beam.L./beam.b./beam.t;             
beam.rhom = beam.rho.*beam.b.*beam.t;                 

%% CREATE TREE STRUCTURE

% 3-D printing parameters

beam.seed{1} = {'F','V','E'};                                       % Printing orientation (Flat, Vertical, On Edge)
beam.seed{2} = {'99','90' '85','70','60','50','40'};                % Infill percentage (%) 
beam.seed{3} = {'00','15','30','45','60','90'};                     % Raster angle (deg)
beam.seed{4} = {'1','2','3','4'};                                   % Layer thickness (1E-4 m)
beam.seed{5} = {'ZZ','T'};                                          % Infill pattern (Zig-Zag, Triangular)   
beam.seed{6} = {'a','b','c'};                                       % Sets of the same beam (3 models of each one)

casestoread = {'E99001ZZa','E90451ZZa','E85001ZZa','E70001Ta','E60001Ta','E50001Ta','E40001Ta','E99002ZZa','E99003ZZa','E99004ZZa',...
               'E99001ZZb','E90451ZZb','E85001ZZb','E70001Tb','E60001Tb','E50001Tb','E99002ZZb','E99003ZZb','E99004ZZb',...
               'E99001ZZc','E90451ZZc','E85001ZZc','E70001Tc','E60001Tc','E50001Tc','E99002ZZc','E99003ZZc','E99004ZZc'};    % Specify the files that are going to be read here

% Generate all possible combinations of beams

cont = 1;
for i=1:length(beam.seed{1})
    for j=1:length(beam.seed{2})
        for k=1:length(beam.seed{3})
            for l=1:length(beam.seed{4})
                for m=1:length(beam.seed{5})
                    for n=1:length(beam.seed{6})
                        beam.name{cont,1} = [beam.seed{1}{i} beam.seed{2}{j} beam.seed{3}{k} beam.seed{4}{l} beam.seed{5}{m} beam.seed{6}{n}];
                        cont = cont + 1;
                    end
                end
            end
        end
    end
end

for s=1:length(casestoread)
    filecontent = readmatrix(['data\' casestoread{s} '.csv'],'Delimiter',';','NumHeaderLines',20,'DecimalSeparator',',');
    test.(casestoread{s}).freq = filecontent(:,1);
    test.(casestoread{s}).real = filecontent(:,2);
    test.(casestoread{s}).imag = filecontent(:,3);
    
    test.(casestoread{s}).amp = sqrt(filecontent(:,2).^2 + filecontent(:,3).^2);
    
    % Findpeaks

    tolook_f = test.(casestoread{s}).freq(10:end);
    tolook_v = log(test.(casestoread{s}).amp(10:end));
    [test.(casestoread{s}).peak, test.(casestoread{s}).rf] = findpeaks(tolook_v,tolook_f,'MinPeakDistance',800,'MinPeakProminence',.1,'NPeaks',3);
    tolook_f = test.(casestoread{s}).freq(10:end);
    tolook_v = log(1./test.(casestoread{s}).amp(10:end));
    [test.(casestoread{s}).min, test.(casestoread{s}).af] = findpeaks(tolook_v,tolook_f,'MinPeakDistance',500,'MinPeakProminence',.78,'NPeaks',3);                                           % Find the antiresonance frequencies
    test.(casestoread{s}).min = 1./test.(casestoread{s}).min;

    % Compute the elastic module with three different methods

    test.(casestoread{s}).E_iso = diag(E_ISO(beam.L(s),beam.rhom(s),beam.Ixx(s),test.(casestoread{s}).af));
    test.(casestoread{s}).E_af = diag(E_AFF(beam.L(s),beam.rhom(s),beam.Ixx(s),test.(casestoread{s}).af));
    test.(casestoread{s}).E_rf = diag(E_RFF(beam.L(s),beam.rhom(s),beam.Ixx(s),test.(casestoread{s}).rf));

end

%% FIGURES
figure()
tiledlayout(7,4)
for t=1:length(casestoread)
    nexttile
    hold on
    plot(test.(casestoread{t}).freq,log(test.(casestoread{t}).amp))
    plot(test.(casestoread{t}).rf,test.(casestoread{t}).peak,'ro')
    plot(test.(casestoread{t}).af,test.(casestoread{t}).min,'bo')
    title(casestoread{t})

end

%% EXPORT RESULTS TO .CSV FILE

for i=1:length(casestoread)
    T(i,:) = table(casestoread(i),beam.m(i),beam.L(i),beam.b(i),beam.t(i),test.(casestoread{i}).rf',test.(casestoread{i}).af', ...
        test.(casestoread{i}).E_iso',test.(casestoread{i}).E_af',test.(casestoread{i}).E_rf');
end
T = renamevars(T,["Var1","Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10"], ...
    ["Beam","Mass (kg)","Length (m)","Width (m)","Thickness (m)","Resonance Frequencies (Hz)","Antirresonance Frequencies (Hz)","E_ISO (GPa)","E_AF (GPa)","E_RF (GPa)"]);

filename = 'Set1.xlsx';              % Rename the file name when changing the data files
writetable(T,filename,"WriteMode","append","PreserveFormat",true)
%% FUNCTIONS

function E_I = E_ISO(L,rhom,Ixx,f)                           % Young's modulus with ISO-16940 method
    lambda = [1.87510 4.69410 7.85476];
    E_I = rhom/Ixx*(2*pi*(L/2)^2*f./lambda.^2).^2*1E-9;
end

function E_A = E_AFF(L,rhom,Ixx,f)                           % Young's modulus with antiresonance frequencies of a free-free beam method
    lambda = [3.75038 9.39740 15.73438];
    E_A = rhom/Ixx*(2*pi*L^2*f./lambda.^2).^2*1E-9;
end

function E_R = E_RFF(L,rhom,Ixx,f)                           % Young's modulus with resonance frequencies of a free-free beam method
    lambda = [4.73004 10.99561 17.27876];
    E_R = rhom/Ixx*(2*pi*L^2*f./lambda.^2).^2*1E-9;
end

% function varargout = peakfinder(x0, sel, thresh, extrema, includeEndpoints, interpolate)
% %PEAKFINDER Noise tolerant fast peak finding algorithm
% %   INPUTS:
% %       x0 - A real vector from the maxima will be found (required)
% %       sel - The amount above surrounding data for a peak to be,
% %           identified (default = (max(x0)-min(x0))/4). Larger values mean
% %           the algorithm is more selective in finding peaks.
% %       thresh - A threshold value which peaks must be larger than to be
% %           maxima or smaller than to be minima.
% %       extrema - 1 if maxima are desired, -1 if minima are desired
% %           (default = maxima, 1)
% %       includeEndpoints - If true the endpoints will be included as
% %           possible extrema otherwise they will not be included
% %           (default = true)
% %       interpolate - If true quadratic interpolation will be performed
% %           around each extrema to estimate the magnitude and the
% %           position of the peak in terms of fractional indicies. Note that
% %           unlike the rest of this function interpolation assumes the
% %           input is equally spaced. To recover the x_values of the input
% %           rather than the fractional indicies you can do:
% %           peakX = x0 + (peakLoc - 1) * dx
% %           where x0 is the first x value and dx is the spacing of the
% %           vector. Output peakMag to recover interpolated magnitudes.
% %           See example 2 for more information.
% %           (default = false)
% %
% %   OUTPUTS:
% %       peakLoc - The indicies of the identified peaks in x0
% %       peakMag - The magnitude of the identified peaks
% %
% %   [peakLoc] = peakfinder(x0) returns the indicies of local maxima that
% %       are at least 1/4 the range of the data above surrounding data.
% %
% %   [peakLoc] = peakfinder(x0,sel) returns the indicies of local maxima
% %       that are at least sel above surrounding data.
% %
% %   [peakLoc] = peakfinder(x0,sel,thresh) returns the indicies of local
% %       maxima that are at least sel above surrounding data and larger
% %       (smaller) than thresh if you are finding maxima (minima).
% %
% %   [peakLoc] = peakfinder(x0,sel,thresh,extrema) returns the maxima of the
% %       data if extrema > 0 and the minima of the data if extrema < 0
% %
% %   [peakLoc] = peakfinder(x0,sel,thresh,extrema, includeEndpoints)
% %       returns the endpoints as possible extrema if includeEndpoints is
% %       considered true in a boolean sense
% %
% %   [peakLoc, peakMag] = peakfinder(x0,sel,thresh,extrema,interpolate)
% %       returns the results of results of quadratic interpolate around each
% %       extrema if interpolate is considered to be true in a boolean sense
% %
% %   [peakLoc, peakMag] = peakfinder(x0,...) returns the indicies of the
% %       local maxima as well as the magnitudes of those maxima
% %
% %   If called with no output the identified maxima will be plotted along
% %       with the input data.
% %
% %   Note: If repeated values are found the first is identified as the peak
% %
% % Example 1:
% % t = 0:.0001:10;
% % x = 12*sin(10*2*pi*t)-3*sin(.1*2*pi*t)+randn(1,numel(t));
% % x(1250:1255) = max(x);
% % peakfinder(x)
% %
% % Example 2:
% % ds = 100;  % Downsample factor
% % dt = .001; % Time step
% % ds_dt = ds*dt; % Time delta after downsampling
% % t0 = 1;
% % t = t0:dt:5 + t0;
% % x = 0.2-sin(0.01*2*pi*t)+3*cos(7/13*2*pi*t+.1)-2*cos((1+pi/10)*2*pi*t+0.2)-0.2*t;
% % x(end) = min(x);
% % x_ds = x(1:ds:end); % Downsample to test interpolation
% % [minLoc, minMag] = peakfinder(x_ds, .8, 0, -1, false, true);
% % minT = t0 + (minLoc - 1) * ds_dt; % Take into account 1 based indexing
% % p = plot(t,x,'-',t(1:ds:end),x_ds,'o',minT,minMag,'rv');
% % set(p(2:end), 'linewidth', 2); % Show the markers more clearly
% % legend('Actual Data', 'Input Data', 'Estimated Peaks');
% % Copyright Nathanael C. Yoder 2015 (nyoder@gmail.com)
% % Perform error checking and set defaults if not passed in
% narginchk(1, 6);
% nargoutchk(0, 2);
% s = size(x0);
% flipData =  s(1) < s(2);
% len0 = numel(x0);
% if len0 ~= s(1) && len0 ~= s(2)
%     error('PEAKFINDER:Input','The input data must be a vector')
% elseif isempty(x0)
%     varargout = {[],[]};
%     return;
% end
% if ~isreal(x0)
%     warning('PEAKFINDER:NotReal','Absolute value of data will be used')
%     x0 = abs(x0);
% end
% if nargin < 2 || isempty(sel)
%     sel = (max(x0)-min(x0))/4;
% elseif ~isnumeric(sel) || ~isreal(sel)
%     sel = (max(x0)-min(x0))/4;
%     warning('PEAKFINDER:InvalidSel',...
%         'The selectivity must be a real scalar.  A selectivity of %.4g will be used',sel)
% elseif numel(sel) > 1
%     warning('PEAKFINDER:InvalidSel',...
%         'The selectivity must be a scalar.  The first selectivity value in the vector will be used.')
%     sel = sel(1);
% end
% if nargin < 3 || isempty(thresh)
%     thresh = [];
% elseif ~isnumeric(thresh) || ~isreal(thresh)
%     thresh = [];
%     warning('PEAKFINDER:InvalidThreshold',...
%         'The threshold must be a real scalar. No threshold will be used.')
% elseif numel(thresh) > 1
%     thresh = thresh(1);
%     warning('PEAKFINDER:InvalidThreshold',...
%         'The threshold must be a scalar.  The first threshold value in the vector will be used.')
% end
% if nargin < 4 || isempty(extrema)
%     extrema = 1;
% else
%     extrema = sign(extrema(1)); % Should only be 1 or -1 but make sure
%     if extrema == 0
%         error('PEAKFINDER:ZeroMaxima','Either 1 (for maxima) or -1 (for minima) must be input for extrema');
%     end
% end
% if nargin < 5 || isempty(includeEndpoints)
%     includeEndpoints = true;
% end
% if nargin < 6 || isempty(interpolate)
%     interpolate = false;
% end
% x0 = extrema*x0(:); % Make it so we are finding maxima regardless
% thresh = thresh*extrema; % Adjust threshold according to extrema.
% dx0 = diff(x0); % Find derivative
% dx0(dx0 == 0) = -eps; % This is so we find the first of repeated values
% ind = find(dx0(1:end-1).*dx0(2:end) < 0)+1; % Find where the derivative changes sign
% % Include endpoints in potential peaks and valleys as desired
% if includeEndpoints
%     x = [x0(1);x0(ind);x0(end)];
%     ind = [1;ind;len0];
%     minMag = min(x);
%     leftMin = minMag;
% else
%     x = x0(ind);
%     minMag = min(x);
%     leftMin = min(x(1), x0(1));
% end
% % x only has the peaks, valleys, and possibly endpoints
% len = numel(x);
% if len > 2 % Function with peaks and valleys
%     % Set initial parameters for loop
%     tempMag = minMag;
%     foundPeak = false;
%     if includeEndpoints
%         % Deal with first point a little differently since tacked it on
%         % Calculate the sign of the derivative since we tacked the first
%         %  point on it does not neccessarily alternate like the rest.
%         signDx = sign(diff(x(1:3)));
%         if signDx(1) <= 0 % The first point is larger or equal to the second
%             if signDx(1) == signDx(2) % Want alternating signs
%                 x(2) = [];
%                 ind(2) = [];
%                 len = len-1;
%             end
%         else % First point is smaller than the second
%             if signDx(1) == signDx(2) % Want alternating signs
%                 x(1) = [];
%                 ind(1) = [];
%                 len = len-1;
%             end
%         end
%     end
%     % Skip the first point if it is smaller so we always start on a
%     %   maxima
%     if x(1) >= x(2)
%         ii = 0;
%     else
%         ii = 1;
%     end
%     % Preallocate max number of maxima
%     maxPeaks = ceil(len/2);
%     peakLoc = zeros(maxPeaks,1);
%     peakMag = zeros(maxPeaks,1);
%     cInd = 1;
%     % Loop through extrema which should be peaks and then valleys
%     while ii < len
%         ii = ii+1; % This is a peak
%         % Reset peak finding if we had a peak and the next peak is bigger
%         %   than the last or the left min was small enough to reset.
%         if foundPeak
%             tempMag = minMag;
%             foundPeak = false;
%         end
%         % Found new peak that was lager than temp mag and selectivity larger
%         %   than the minimum to its left.
%         if x(ii) > tempMag && x(ii) > leftMin + sel
%             tempLoc = ii;
%             tempMag = x(ii);
%         end
%         % Make sure we don't iterate past the length of our vector
%         if ii == len
%             break; % We assign the last point differently out of the loop
%         end
%         ii = ii+1; % Move onto the valley
%         % Come down at least sel from peak
%         if ~foundPeak && tempMag > sel + x(ii)
%             foundPeak = true; % We have found a peak
%             leftMin = x(ii);
%             peakLoc(cInd) = tempLoc; % Add peak to index
%             peakMag(cInd) = tempMag;
%             cInd = cInd+1;
%         elseif x(ii) < leftMin % New left minima
%             leftMin = x(ii);
%         end
%     end
%     % Check end point
%     if includeEndpoints
%         if x(end) > tempMag && x(end) > leftMin + sel
%             peakLoc(cInd) = len;
%             peakMag(cInd) = x(end);
%             cInd = cInd + 1;
%         elseif ~foundPeak && tempMag > minMag % Check if we still need to add the last point
%             peakLoc(cInd) = tempLoc;
%             peakMag(cInd) = tempMag;
%             cInd = cInd + 1;
%         end
%     elseif ~foundPeak
%         if x(end) > tempMag && x(end) > leftMin + sel
%             peakLoc(cInd) = len;
%             peakMag(cInd) = x(end);
%             cInd = cInd + 1;
%         elseif tempMag > min(x0(end), x(end)) + sel
%             peakLoc(cInd) = tempLoc;
%             peakMag(cInd) = tempMag;
%             cInd = cInd + 1;
%         end
%     end
%     % Create output
%     if cInd > 1
%         peakInds = ind(peakLoc(1:cInd-1));
%         peakMags = peakMag(1:cInd-1);
%     else
%         peakInds = [];
%         peakMags = [];
%     end
% else % This is a monotone function where an endpoint is the only peak
%     [peakMags,xInd] = max(x);
%     if includeEndpoints && peakMags > minMag + sel
%         peakInds = ind(xInd);
%     else
%         peakMags = [];
%         peakInds = [];
%     end
% end
% % Apply threshold value.  Since always finding maxima it will always be
% %   larger than the thresh.
% if ~isempty(thresh)
%     m = peakMags>thresh;
%     peakInds = peakInds(m);
%     peakMags = peakMags(m);
% end
% if interpolate && ~isempty(peakMags)
%     middleMask = (peakInds > 1) & (peakInds < len0);
%     noEnds = peakInds(middleMask);
%     magDiff = x0(noEnds + 1) - x0(noEnds - 1);
%     magSum = x0(noEnds - 1) + x0(noEnds + 1)  - 2 * x0(noEnds);
%     magRatio = magDiff ./ magSum;
%     peakInds(middleMask) = peakInds(middleMask) - magRatio/2;
%     peakMags(middleMask) = peakMags(middleMask) - magRatio .* magDiff/8;
% end
% % Rotate data if needed
% if flipData
%     peakMags = peakMags.';
%     peakInds = peakInds.';
% end
% % Change sign of data if was finding minima
% if extrema < 0
%     peakMags = -peakMags;
%     x0 = -x0;
% end
% % Plot if no output desired
% if nargout == 0
%     if isempty(peakInds)
%         disp('No significant peaks found')
%     else
%         figure;
%         plot(1:len0,x0,'.-',peakInds,peakMags,'ro','linewidth',2);
%     end
% else
%     varargout = {peakInds,peakMags};
% end
% end