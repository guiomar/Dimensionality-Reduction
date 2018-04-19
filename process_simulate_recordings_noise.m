function varargout = process_simulate_recordings_noise( varargin )
% PROCESS_SIMULATE_SOURCES_NOISE: Simulate source files based on some scouts and add noise
% at the level of the sources (random noise, SNR1) or at the level of the sensors 
% (from a noise covariance matrix, SNR2).
%
% USAGE:  OutputFiles = process_simulate_recordings_noise('Run', sProcess, sInputA)
 
% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2014 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
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
% Authors: Guiomar Niso, Francois Tadel 2014
%
% It calls: get_noise_signals.m

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Simulate recordings from scouts (noise)';
    sProcess.FileTag     = '';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = 'Simulate'; 
    sProcess.Index       = 915; 
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'matrix'};
    sProcess.OutputTypes = {'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    % === SCOUTS
    sProcess.options.scouts.Comment = '';
    sProcess.options.scouts.Type    = 'scout';
    sProcess.options.scouts.Value   = {};
    % === SAVE SOURCES
    sProcess.options.savesources.Comment = 'Save full sources';
    sProcess.options.savesources.Type    = 'checkbox';
    sProcess.options.savesources.Value   = 0;
    % === LEVEL OF NOISE (SNR1)
    sProcess.options.noise1.Comment = 'Level of random noise (SNR1):';
    sProcess.options.noise1.Type    = 'value';
    sProcess.options.noise1.Value   = {0, '', 2};
    % === LEVEL OF SENSOR NOISE (SNR2)
    sProcess.options.noise2.Comment = 'Level of sensor noise (SNR2):';
    sProcess.options.noise2.Type    = 'value';
    sProcess.options.noise2.Value   = {0, '', 2};

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput) %#ok<DEFNU>
    OutputFiles = {};

    % === GET OPTIONS ===
    
    % Get other optinos
    SaveSources = sProcess.options.savesources.Value;
    SNR1 = sProcess.options.noise1.Value{1};
    SNR2 = sProcess.options.noise2.Value{1};

    % Get scouts
    AtlasList = sProcess.options.scouts.Value;
    if isempty(AtlasList)
        bst_report('Error', sProcess, [], 'No scouts selected.');
        return;
    end
    
    % === LOAD CHANNEL FILE / HEAD MODEL===
    % Get condition
    sStudy = bst_get('Study', sInput.iStudy);
    % Get channel file
    [sChannel, iStudyChannel] = bst_get('ChannelForStudy', sInput.iStudy);
    if isempty(sChannel)
        bst_report('Error', sProcess, [], ['No channel file available.' 10 'Please import a channel file in this study before running simulations.']);
        return;
    end
    % Get study channel
    sStudyChannel = bst_get('Study', iStudyChannel);
    % Check head model
    if isempty(sStudyChannel.iHeadModel)
        bst_report('Error', sProcess, [], ['No head model file available.' 10 'Please calculate a head model before running simulations.']);
        return;
    end
    % Check noisecov
    if isempty(sStudyChannel.NoiseCov)
        bst_report('Error', sProcess, [], ['No noise covariance file available.' 10 'Please import a noise covariance before running simulations.']);
        return;
    end
    % Load channel file
    ChannelMat = in_bst_channel(sChannel.FileName);
%     Nchannels = sum(strcmp({ChannelMat.Channel(:).Type},'MEG')); % MEG channels
    % Load the noise covariance matrix
    NoiseCovMat = load(file_fullpath(sStudyChannel.NoiseCov.FileName));
    % Load head model
    HeadModelFile = sStudyChannel.HeadModel(sStudyChannel.iHeadModel).FileName;
    HeadModelMat = in_headmodel_bst(HeadModelFile);
    % If no orientations: error
    if isempty(HeadModelMat.GridOrient)
        bst_report('Error', sProcess, [], 'No source orientations available in this head model.');
        return;
    end
    % Apply the fixed orientation to the Gain matrix (normal to the cortex)
    HeadModelMat.Gain = bst_gain_orient(HeadModelMat.Gain, HeadModelMat.GridOrient);
    % Get all the MEG/EEG channels
    Modalities = {};
    if ~isempty(HeadModelMat.MEGMethod)
        Modalities{end+1} = 'MEG';
    end
    if ~isempty(HeadModelMat.EEGMethod)
        Modalities{end+1} = 'EEG';
    end
    if ~isempty(HeadModelMat.SEEGMethod)
        Modalities{end+1} = 'SEEG';
    end
    if ~isempty(HeadModelMat.ECOGMethod)
        Modalities{end+1} = 'ECOG';
    end
    iChannels = channel_find(ChannelMat.Channel, Modalities);
    
    
    % === LOAD CORTEX ===
    % Get surface from the head model
    SurfaceFile = HeadModelMat.SurfaceFile;
    % Load surface
    sSurface = in_tess_bst(SurfaceFile);
    % Get scout structures
    sScouts = process_extract_scout('GetScoutsInfo', sProcess, sInput, SurfaceFile, AtlasList);
    if isempty(sScouts)
        return;
    end
    
    % === LOAD INPUT FILE ===
    % Read input file
    sMatrix = in_bst_matrix(sInput.FileName);
    % Check dimensions
    if (length(sScouts) ~= size(sMatrix.Value,1))
        bst_report('Error', sProcess, [], sprintf('The number of selected scouts (%d) does not match the number of signals (%d).', length(sScouts), size(sMatrix.Value,1)));
        return;
    end

    % === GENERATE SOURCE MATRIX ===
    nSources = length(sSurface.Vertices);
    nTime = size(sMatrix.Value,2);
    % Initialize space matrix
    ImageGridAmp = sparse([],[],[],nSources, nTime, length([sScouts.Vertices])*nTime);

    
    % Fill matrix
    for i = 1:length(sScouts)
        ImageGridAmp(sScouts(i).Vertices,:) = repmat(sMatrix.Value(i,:), length(sScouts(i).Vertices), 1);
    end
    
    %%% NOISE ON THE SOURCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Add noise SNR1 (random noise on the sources) 
    ImageGridAmp = ImageGridAmp + SNR1.*(rand(size(ImageGridAmp))-0.5).*max(max(abs(ImageGridAmp)));
    % Set unit range to pAm
    ImageGridAmp = 1e-9 .* ImageGridAmp;
    
    % === SAVE RECORDINGS ===
    % Generate data matrix
    F = zeros(length(ChannelMat.Channel), nTime);
    F(iChannels,:) = HeadModelMat.Gain(iChannels,:) * ImageGridAmp;
    
    %%% NOISE ON THE SENSORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Add noise SNR2 (sensor noise) 
    xn = get_noise_signals (NoiseCovMat.NoiseCov(iChannels,iChannels), nTime);
    xnn = xn./max(max(xn)); % Noise signal between 0 and 1
    xns = xnn.*max(max(F(iChannels,:))); % Make the noise of similar amplitude than the signal
    F(iChannels,:) = F(iChannels,:) + SNR2*xns; % Apply the SNR2

    
    
    % Create a new data file structure
    DataMat = db_template('datamat');
    DataMat.F           = F;
    DataMat.Comment     = [sMatrix.Comment,',Nsn:',num2str(SNR1),',Nsc:',num2str(SNR2)];
    DataMat.ChannelFlag = strcmp({ChannelMat.Channel(:).Type},'MEG'); % ones(length(ChannelMat.Channel), 1);
    DataMat.Time        = sMatrix.Time;
    DataMat.DataType    = 'recordings';
    DataMat.Device      = 'simulation';
    DataMat.nAvg        = 1;
    DataMat.Events      = [];
    % Add history entry
    DataMat = bst_history('add', DataMat, 'simulate', ['Simulated from file: ' sInput.FileName]);
    % Output filename
    DataFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_sim');
    % Save on disk
    bst_save(DataFile, DataMat, 'v6');
    % Register in database
    db_add_data(sInput.iStudy, DataFile, DataMat);
    % Return data file
    OutputFiles{1} = DataFile;
    
    % === SAVE SOURCE FILE ===
    if SaveSources
        % Create a new source file structure
        ResultsMat = db_template('resultsmat');
        ResultsMat.ImagingKernel = [];
        ResultsMat.ImageGridAmp  = full(ImageGridAmp);
        ResultsMat.nComponents   = 1;
        ResultsMat.Comment       = sMatrix.Comment;
        ResultsMat.Function      = 'Simulation';
        ResultsMat.Time          = sMatrix.Time;
        ResultsMat.DataFile      = file_short(DataFile);
        ResultsMat.HeadModelFile = HeadModelFile;
        ResultsMat.HeadModelType = HeadModelMat.HeadModelType;
        ResultsMat.ChannelFlag   = [];
        ResultsMat.GoodChannel   = iChannels;
        ResultsMat.SurfaceFile   = SurfaceFile;
        ResultsMat.GridLoc       = [];
        % Add history entry
        ResultsMat = bst_history('add', ResultsMat, 'simulate', ['Simulated from file: ' sInput.FileName]);
        % Output filename
        ResultsFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'results_sim');
        % Save on disk
        bst_save(ResultsFile, ResultsMat, 'v6');
        % Register in database
        db_add_data(sInput.iStudy, ResultsFile, ResultsMat);
    end
end




