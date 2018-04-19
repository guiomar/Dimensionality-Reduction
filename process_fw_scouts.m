function varargout = process_fw_scouts( varargin )
% PROCESS_SIMULATE_SOURCES: Simulate source files based on some scouts.
%
% USAGE:  OutputFiles = process_fw_scouts('Run', sProcess, sInputA)
 
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
% Authors: Guiomar Niso 2014

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'FW model scouts';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom'; 
    sProcess.SubGroup    = 'Simulate'; 
    sProcess.Index       = 917; 
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1; % No imput file needed
    sProcess.isSeparator = 0;

    % === GROUND TRUTH
    sProcess.options.atlas.Comment = 'Original Scouts to threshold:';
    sProcess.options.atlas.Type    = 'atlas';
    sProcess.options.atlas.Value   = [];
    % === ALL
    sProcess.options.all.Comment = 'Compute all the scouts at the same time';
    sProcess.options.all.Type    = 'checkbox';
    sProcess.options.all.Value   = 0;

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput) %#ok<DEFNU>
    OutputFiles = {};

    % === GET OPTIONS ===    
   
    all = sProcess.options.all.Value;

    % Load the surface filename from results file
    sResults = in_bst_results(sInput.FileName, 0);

    
    % === LOAD HEAD MODEL===
    
    % Get condition
%     sStudy = bst_get('Study', sInput.iStudy);
    % Get channel file
    [sChannel, iStudyChannel] = bst_get('ChannelForStudy', sInput.iStudy);
    if isempty(sChannel)
        bst_report('Error', sProcess, [], ['No channel file available.' 10 'Please import a channel file in this study before running simulations.']);
        return;
    end
    % Get study channel
    sStudyChannel = bst_get('Study', iStudyChannel);
    % Check head model and noisecov
    if isempty(sStudyChannel.iHeadModel)
        bst_report('Error', sProcess, [], ['No head model file available.' 10 'Please calculate a head model before running simulations.']);
        return;
    end
    if isempty(sStudyChannel.NoiseCov)
        bst_report('Error', sProcess, [], ['No noise covariance file available.' 10 'Please import a noise covariance before running simulations.']);
        return;
    end
    % Load channel file
    ChannelMat = in_bst_channel(sChannel.FileName);
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
    
    % Get the original atlas
    iAtlas = [];
    if ~isempty(sProcess.options.atlas.Value)
        iAtlas = find(strcmpi({sSurface.Atlas.Name}, sProcess.options.atlas.Value));
    end
    if isempty(iAtlas),   iAtlas = sSurface.iAtlas;    end
    if isempty(iAtlas),   iAtlas = 1;    end

    Nvertices = size(sSurface.Vertices,1);
    ATLAS  = zeros(Nvertices, 1);
    
    % Mask of the scouts
    for i=1:numel(sSurface.Atlas(iAtlas).Scouts)
        ATLAS( sSurface.Atlas(iAtlas).Scouts(i).Vertices ) = i; 
    end    
 
    
    
    % === GENERATE SOURCE MATRIX ===
         
    if all % == process all scouts at the same time
        
        %%%%% IGA-IK(Nsources x Nsensors) * HM(Nsensors x Nsources) * ATLAS(Nsources x 1)
        % (1e-9) Set unit range to pAm
        
        %%% Si esta el LINK a las sources
        THRESMAP = sResults.ImagingKernel * HeadModelMat.Gain(iChannels,:) * (ATLAS>0); 
        %%% Si estan las sources tal cual de la simulacion
        % THRESMAP = sResults.ImageGridAmp(:,iChannels) * HeadModelMat.Gain(iChannels,:) * (1e-9 .*ATLAS);
        
        % Automatic threshold
%         Labels = abs(THRESMAP) > 0.2*max(abs(THRESMAP));

        % Same number of vertices than the original scout
        NverticesON = sum(ATLAS>0);
        sorThres = sort(abs(THRESMAP),'descend');
        Labels = abs(THRESMAP) >= sorThres(NverticesON);

    else  % == process scouts one by one
        
        Nscouts = numel(sSurface.Atlas(iAtlas).Scouts);
        Labels = zeros (Nvertices, Nscouts);
        for i=1:Nscouts
            THRESMAP = sResults.ImagingKernel * HeadModelMat.Gain(iChannels,:) * (ATLAS==i); 
            
%             % Same number of vertices than the original scout
            NverticesON = sum(ATLAS==i);
            sorThres = sort(abs(THRESMAP),'descend');
            Labels(abs(THRESMAP) >= sorThres(NverticesON),i)  = i;

            % Automatic threshold
%             Labels(:,i) = i .* (abs(THRESMAP) > 0.4*max(abs(THRESMAP)) );
        end
        
        Labels = sum(Labels,2); % becareful with overlaps, find a better way to keep the colors

    end
     
    %  figure; hist(abs(THRESMAP))
    %  hold on; plot(0.3*max(abs(THRESMAP)),0,'*r')
    %  plot(sorThres(NverticesON),0,'*g')


    % ===== FILL NEW ATLAS =====
    
    LabelFunction = 'mean';
    
    sSurface = fill_new_atlas ( Labels, sSurface, 'FWM', LabelFunction );
    % Save the modifications back to the surface file
    bst_save(file_fullpath(sResults.SurfaceFile), sSurface, 'v6');
    
   
   
    % ===== OUTPUT =====
    OutputFiles{1} = sInput.FileName;

end




