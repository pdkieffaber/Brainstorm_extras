function varargout = process_interp( varargin )
% PROCESS_interp: Interpolate Bad Electrodes 
%
% USAGE:                sProcess = process_interp('GetDescription')
%                         sInput = process_interp('Run', sProcess, sInput, method=[])
%        [x, FiltSpec, Messages] = process_interp('Compute', x, method, badchans)
%                              x = process_bandpass('Compute', x, method, badchans)

% @=============================================================================
% This function is part of the CPL Toolbox for Brainstorm Software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c)2020 College of William & Mary
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Paul Kieffaber (I authored the brainstorm process, but
% subfunctions that do all the real work were authored by Arnaud Delorme

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Interpolate Channels';
sProcess.FileTag     = @GetFileTag;
sProcess.Category    = 'Filter';
sProcess.SubGroup    = {'CPL','Interpolate Channels'};
sProcess.Index       = 66;
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data', 'results', 'raw', 'matrix'};
sProcess.OutputTypes = {'data', 'results', 'raw', 'matrix'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
sProcess.processDim  = 1;   % Process channel by channel
sProcess.Description = 'https://neuroimage.usc.edu/brainstorm/Tutorials/ArtifactsFilter#What_filters_to_apply.3F';
% Definition of the options
% === Sensor types
sProcess.options.sensortypes.Comment = 'Sensor types or names (empty=all): ';
sProcess.options.sensortypes.Type    = 'text';
sProcess.options.sensortypes.Value   = 'MEG, EEG';
sProcess.options.sensortypes.InputTypes = {'data', 'raw'};
% ==== Parameters
sProcess.options.label1.Comment = '<BR><U><B>Interpolation Parameters</B></U>:';
sProcess.options.label1.Type    = 'label';
% === Interp method
sProcess.options.method.Comment = {'Spherical Spline','Nearest Neighbor'};
sProcess.options.method.Type    = 'radio';
sProcess.options.method.Value   = 1;

end

%% ===== FORMAT COMMENT =====
function [Comment, fileTag] = FormatComment(sProcess)
    % Get options
    [method] = GetOptions(sProcess);
    % Format comment
    if ~isempty(method)
        Comment = ['Interp:spherical'];
        fileTag = 'interp';
    else
        Comment = '';
    end
end

%% ===== GET FILE TAG =====
function fileTag = GetFileTag(sProcess)
    [Comment, fileTag] = FormatComment(sProcess);
end

%% ===== GET OPTIONS =====
function [method] = GetOptions(sProcess)
% Method selection
switch (sProcess.options.method.Value)
    case 1
        method = 'spherical';
    case 2
        disp('Not yet implemented')
        return
end
end

%% ===== RUN =====
function sInput = Run(sProcess, sInput) %#ok<DEFNU>
% Get options
[method] = GetOptions(sProcess);

%Get Polar Channel Coords
disp(['load(' sInput.ChannelFile ',Channel)'])
eval(['load(''/Volumes/Data/Brainstorm_DB/SM_2020/data/' sInput.ChannelFile ''',''Channel'')']);
badchans=find(sInput.ChannelFlag==-1);
% Filter signals
% sfreq = 1 ./ (sInput.TimeVector(2) - sInput.TimeVector(1));
% [sInput.A, FiltSpec, Messages] = Compute(sInput.A, sfreq, HighPass, LowPass, Method, isMirror, isRelax, TranBand);

%Interpolate Channels
[sInput.A, Messages] = Compute(sInput.A, method, Channel, badchans);
% % Process warnings
% if ~isempty(Messages)
%     bst_report('Warning', sProcess, sInput, Messages);
% end
% 
% % File comment
% if ~isempty(HighPass) && ~isempty(LowPass)
%     filterComment = ['band(' num2str(HighPass) '-' num2str(LowPass) 'Hz)'];
% elseif ~isempty(HighPass)
%     filterComment = ['high(' num2str(HighPass) 'Hz)'];
% elseif ~isempty(LowPass)
%     filterComment = ['low(', num2str(LowPass) 'Hz)'];
% else
%     filterComment = '';
% end
% sInput.CommentTag = filterComment;
% % Do not keep the Std field in the output
% if isfield(sInput, 'Std') && ~isempty(sInput.Std)
%     sInput.Std = [];
% end
end


%% ===== EXTERNAL CALL =====
% USAGE: [x, FiltSpec, Messages] = process_interp('Compute', x, method, badchans)
%
function [x, Messages] = Compute(x, method, Channel, badchans)
Messages=['I did some interpolation'];
%Interpolate Electrodes
orix = x;
% find non-empty good channels
% ----------------------------
goodchans = setdiff(1:length(Channel),badchans);
%goodchans = goodchans( sort(indgood) );
datachans = getdatachans(goodchans,badchans);
%badchans  = intersect_bc(badchans, nonemptychans);
if isempty(badchans), return; end

%convert chanlocs for use with EEGlab functions
for i=1:length(Channel)
    elec(i).X=Channel(i).Loc(1);
    elec(i).Y=Channel(i).Loc(2);
    elec(i).Z=Channel(i).Loc(3);
end

% scan data points
% ----------------
if strcmpi(method, 'spherical')
    % get theta, rad of electrodes
    % ----------------------------
    tmpgoodlocs = elec(goodchans);
    xelec = [ tmpgoodlocs.X ];
    yelec = [ tmpgoodlocs.Y ];
    zelec = [ tmpgoodlocs.Z ];
    rad = sqrt(xelec.^2+yelec.^2+zelec.^2);
    xelec = xelec./rad;
    yelec = yelec./rad;
    zelec = zelec./rad;
    tmpbadlocs = elec(badchans);
    xbad = [ tmpbadlocs.X ];
    ybad = [ tmpbadlocs.Y ];
    zbad = [ tmpbadlocs.Z ];
    rad = sqrt(xbad.^2+ybad.^2+zbad.^2);
    xbad = xbad./rad;
    ybad = ybad./rad;
    zbad = zbad./rad;
    
    
    %[tmp1 tmp2 tmp3 tmpchans] = spheric_spline_old( xelec, yelec, zelec, EEG.data(goodchans,1));
    %max(tmpchans(:,1)), std(tmpchans(:,1)),
    %[tmp1 tmp2 tmp3 EEG.data(badchans,:)] = spheric_spline( xelec, yelec, zelec, xbad, ybad, zbad, EEG.data(goodchans,:));
    [tmp1 tmp2 tmp3 badchansdata] = spheric_spline( xelec, yelec, zelec, xbad, ybad, zbad, x(datachans,:));
    %max(EEG.data(goodchans,1)), std(EEG.data(goodchans,1))
    %max(EEG.data(badchans,1)), std(EEG.data(badchans,1))
    
elseif strcmpi(method, 'spacetime') % 3D interpolation, works but x10 times slower
    disp('Warning: if processing epoch data, epoch boundary are ignored...');
    disp('3-D interpolation, this can take a long (long) time...');
    tmpgoodlocs = EEG.chanlocs(goodchans);
    tmpbadlocs = EEG.chanlocs(badchans);
    [xbad ,ybad]  = pol2cart([tmpbadlocs.theta],[tmpbadlocs.radius]);
    [xgood,ygood] = pol2cart([tmpgoodlocs.theta],[tmpgoodlocs.radius]);
    pnts = size(EEG.data,2)*size(EEG.data,3);
    zgood = [1:pnts];
    zgood = repmat(zgood, [length(xgood) 1]);
    zgood = reshape(zgood,prod(size(zgood)),1);
    xgood = repmat(xgood, [1 pnts]); xgood = reshape(xgood,prod(size(xgood)),1);
    ygood = repmat(ygood, [1 pnts]); ygood = reshape(ygood,prod(size(ygood)),1);
    tmpdata = reshape(EEG.data, prod(size(EEG.data)),1);
    zbad = 1:pnts;
    zbad = repmat(zbad, [length(xbad) 1]);
    zbad = reshape(zbad,prod(size(zbad)),1);
    xbad = repmat(xbad, [1 pnts]); xbad = reshape(xbad,prod(size(xbad)),1);
    ybad = repmat(ybad, [1 pnts]); ybad = reshape(ybad,prod(size(ybad)),1);
    badchansdata = griddata3(ygood, xgood, zgood, tmpdata,...
        ybad, xbad, zbad, 'nearest'); % interpolate data
else
    % get theta, rad of electrodes
    % ----------------------------
    tmpchanlocs = EEG.chanlocs;
    [xbad ,ybad]  = pol2cart([tmpchanlocs( badchans).theta],[tmpchanlocs( badchans).radius]);
    [xgood,ygood] = pol2cart([tmpchanlocs(goodchans).theta],[tmpchanlocs(goodchans).radius]);
    
    fprintf('Points (/%d):', size(EEG.data,2)*size(EEG.data,3));
    badchansdata = zeros(length(badchans), size(EEG.data,2)*size(EEG.data,3));
    
    for t=1:(size(EEG.data,2)*size(EEG.data,3)) % scan data points
        if mod(t,100) == 0, fprintf('%d ', t); end
        if mod(t,1000) == 0, fprintf('\n'); end;
        
        %for c = 1:length(badchans)
        %   [h EEG.data(badchans(c),t)]= topoplot(EEG.data(goodchans,t),EEG.chanlocs(goodchans),'noplot', ...
        %        [EEG.chanlocs( badchans(c)).radius EEG.chanlocs( badchans(c)).theta]);
        %end
        tmpdata = reshape(EEG.data, size(EEG.data,1), size(EEG.data,2)*size(EEG.data,3) );
        if strcmpi(method, 'invdist'), method = 'v4'; end
        [Xi,Yi,badchansdata(:,t)] = griddata(ygood, xgood , double(tmpdata(datachans,t)'),...
            ybad, xbad, method); % interpolate data
    end
    fprintf('\n');
end

x(badchans,:)=badchansdata;
% tmpdata = zeros(length(badchans), size(sInput.A,1));
% tmpdata(origoodchans, :,:) = EEG.data;
% %if input data are epoched reshape badchansdata for Octave compatibility...
% if length(size(tmpdata))==3
%     badchansdata = reshape(badchansdata,length(badchans),size(tmpdata,2),size(tmpdata,3));
% end
% tmpdata(badchans,:,:) = badchansdata;
% EEG.data = tmpdata;
% EEG.nbchan = size(EEG.data,1);
% EEG = eeg_checkset(EEG);

% get data channels
% -----------------
    function datachans = getdatachans(goodchans, badchans);
        datachans = goodchans;
        badchans  = sort(badchans);
        for index = length(badchans):-1:1
            datachans(find(datachans > badchans(index))) = datachans(find(datachans > badchans(index)))-1;
        end
    end
% -----------------
% spherical splines
% -----------------
    function [x, y, z, Res] = spheric_spline_old( xelec, yelec, zelec, values);
        
        SPHERERES = 20;
        [x,y,z] = sphere(SPHERERES);
        x(1:(length(x)-1)/2,:) = []; x = [ x(:)' ];
        y(1:(length(y)-1)/2,:) = []; y = [ y(:)' ];
        z(1:(length(z)-1)/2,:) = []; z = [ z(:)' ];
        
        Gelec = computeg(xelec,yelec,zelec,xelec,yelec,zelec);
        Gsph  = computeg(x,y,z,xelec,yelec,zelec);
        
        % equations are
        % Gelec*C + C0  = Potential (C unknow)
        % Sum(c_i) = 0
        % so
        %             [c_1]
        %      *      [c_2]
        %             [c_ ]
        %    xelec    [c_n]
        % [x x x x x]         [potential_1]
        % [x x x x x]         [potential_ ]
        % [x x x x x]       = [potential_ ]
        % [x x x x x]         [potential_4]
        % [1 1 1 1 1]         [0]
        
        % compute solution for parameters C
        % ---------------------------------
        meanvalues = mean(values);
        values = values - meanvalues; % make mean zero
        C = pinv([Gelec;ones(1,length(Gelec))]) * [values(:);0];
        
        % apply results
        % -------------
        Res = zeros(1,size(Gsph,1));
        for j = 1:size(Gsph,1)
            Res(j) = sum(C .* Gsph(j,:)');
        end
        Res = Res + meanvalues;
        Res = reshape(Res, length(x(:)),1);
    end

    function [xbad, ybad, zbad, allres] = spheric_spline( xelec, yelec, zelec, xbad, ybad, zbad, values);
        
        newchans = length(xbad);
        numpoints = size(values,2);
        
        %SPHERERES = 20;
        %[x,y,z] = sphere(SPHERERES);
        %x(1:(length(x)-1)/2,:) = []; xbad = [ x(:)'];
        %y(1:(length(x)-1)/2,:) = []; ybad = [ y(:)'];
        %z(1:(length(x)-1)/2,:) = []; zbad = [ z(:)'];
        
        Gelec = computeg(xelec,yelec,zelec,xelec,yelec,zelec);
        Gsph  = computeg(xbad,ybad,zbad,xelec,yelec,zelec);
        
        % compute solution for parameters C
        % ---------------------------------
        meanvalues = mean(values);
        values = values - repmat(meanvalues, [size(values,1) 1]); % make mean zero
        
        values = [values;zeros(1,numpoints)];
        C = pinv([Gelec;ones(1,length(Gelec))]) * values;
        clear values;
        allres = zeros(newchans, numpoints);
        
        % apply results
        % -------------
        for j = 1:size(Gsph,1)
            allres(j,:) = sum(C .* repmat(Gsph(j,:)', [1 size(C,2)]));
        end
        allres = allres + repmat(meanvalues, [size(allres,1) 1]);
        
        % compute G function
        % ------------------
        function g = computeg(x,y,z,xelec,yelec,zelec)
            
            unitmat = ones(length(x(:)),length(xelec));
            EI = unitmat - sqrt((repmat(x(:),1,length(xelec)) - repmat(xelec,length(x(:)),1)).^2 +...
                (repmat(y(:),1,length(xelec)) - repmat(yelec,length(x(:)),1)).^2 +...
                (repmat(z(:),1,length(xelec)) - repmat(zelec,length(x(:)),1)).^2);
            
            g = zeros(length(x(:)),length(xelec));
            %dsafds
            m = 4; % 3 is linear, 4 is best according to Perrin's curve
            for n = 1:7
                if ismatlab
                    L = legendre(n,EI);
                else % Octave legendre function cannot process 2-D matrices
                    for icol = 1:size(EI,2)
                        tmpL = legendre(n,EI(:,icol));
                        if icol == 1, L = zeros([ size(tmpL) size(EI,2)]); end
                        L(:,:,icol) = tmpL;
                    end
                end
                g = g + ((2*n+1)/(n^m*(n+1)^m))*squeeze(L(1,:,:));
            end
            g = g/(4*pi);
        end
    end
end