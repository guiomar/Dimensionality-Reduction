function varargout = process_cluster_ik( varargin )
% PROCESS_SOURCE_ATLAS: Project a source file on an atlas (one time series per scout).

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
% Authors: Guiomar Niso, 2014

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % ===== PROCESS =====
    % Description the process
    sProcess.Comment     = 'Cluster IK';
    sProcess.FileTag     = '';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = 'Connectivity';
    sProcess.Index       = 650;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 1;
    
    
    % === SCOUTS
    sProcess.options.scouts.Comment = '';
    sProcess.options.scouts.Type    = 'scout';
    sProcess.options.scouts.Value   = {};
    % === SCOUT FUNCTION ===
    sProcess.options.scoutfunc.Comment = {'Mean', 'Max', 'PCA', 'Std', 'Scout function:'};
    sProcess.options.scoutfunc.Type    = 'radio_line';
    sProcess.options.scoutfunc.Value   = 1;
%     % === N ITERATIONS
% %     sProcess.options.niter.Comment = 'Number of iterations for cluster IK algorithm';
% %     sProcess.options.niter.Type    = 'value';
% %     sProcess.options.niter.Value   = {1, '', 0};
    % === THRESHOLD FOR CORRELATION
    sProcess.options.threshold.Comment = 'Threshold value for correlation to cluster IK elements';
    sProcess.options.threshold.Type    = 'value';
    sProcess.options.threshold.Value   = {1, '', 3};
    % === NORM XYZ
    sProcess.options.isnorm.Comment = 'Unconstrained sources: Norm of the three orientations (x,y,z)';
    sProcess.options.isnorm.Type    = 'checkbox';
    sProcess.options.isnorm.Value   = 0;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput) %#ok<DEFNU>

% == INPUTS ====

% % niter = sProcess.options.niter.Value{1};
LIM   = sProcess.options.threshold.Value{1};

% Get scouts
AtlasList = sProcess.options.scouts.Value;

% Convert from older structure (keep for backward compatibility)
if isstruct(AtlasList) && ~isempty(AtlasList)
    AtlasList = {'User scouts', {AtlasList.Label}};
end
% No scouts selected: exit
if isempty(AtlasList) || ~iscell(AtlasList) || (size(AtlasList,2) < 2) || isempty(AtlasList{1,2})
    bst_report('Error', sProcess, [], 'No scout selected.');
    return;
end

% Override scouts function
if ~isempty(sProcess.options.scoutfunc.Value)
    switch lower(sProcess.options.scoutfunc.Value)
        case {1, 'mean'}, ScoutFunc = 'mean';
        case {2, 'max'},  ScoutFunc = 'max';
        case {3, 'pca'},  ScoutFunc = 'pca';
        case {4, 'std'},  ScoutFunc = 'std';
        otherwise,  bst_report('Error', sProcess, [], 'Invalid scout function.');  return;
    end
else
    ScoutFunc = [];
end

% Unconstrained function
isNorm = sProcess.options.isnorm.Value;
if isNorm,    xyzFunction = 'norm';
else          xyzFunction = 'none';     end

% ---------------------------------------
nComponents = 1;
%%%% isFlip?
%  ---------------------------------------

% == LOAD INFO ====

% Load the surface filename from results file
sResults = in_bst_results(sInput.FileName, 0);
% sResultsFull = in_bst_results(sInput.FileName, 1);  %%%%%%% Only for ImageGridAmp !!

% SourceData = sResultsFull.ImageGridAmp; 
% clear sResultsFull

% Get selected scouts and surface
[sScouts, AllAtlasNames, sSurface] = process_extract_scout('GetScoutsInfo', sProcess, sInput, sResults.SurfaceFile, AtlasList);

% nSamples   = size(SourceData,2);
nchannels  = size(sResults.ImagingKernel,2);
nVertices  = size(sSurface.Vertices,1);
Keep       = [sScouts.Vertices];


% COMPUTE CLUSTER IK ======================================================

% % for iter = 1:niter
    
% Get labels for clustering Max IK
% % if iter==1
    % Initialize for 1st iteration
    [Labels1,err] = tess_cluster_ik( sSurface, sResults.ImagingKernel, Keep, LIM);
    Labels2 = zeros(1,nVertices);
    Labels2(Keep) = Labels1;
    
% % else
% %     [Labels1,err] = tess_cluster_ik( sSurface2, ImagingKernel2, Keep, LIM);
% %     % UPDATE Labels 
% %     Labels2 = updateLabels (Labels2, Labels1);
% % end

uLabels2 = unique(Labels2);
nLabels2 = length(uLabels2);

% Downsize to a new surface acording to labels (only if more than 1 iteration)
% % if niter>1
% %     sSurface2 = tess_downsize_labels(sSurface, Labels2); 
% % end

% Average data acording to labels obtained
% DataR2   = zeros(nLabels2, nSamples);
ImagingKernel2 = zeros(nLabels2, nchannels);

for i=1:nLabels2
    % Get the scout vertices
    iV = find( Labels2 == uLabels2(i) );
    % Get the scout orientation
    ScoutOrient = sSurface.VertNormals(iV,:);
    % Get scout value
%     DataR2(i,:) = bst_scout_value( SourceData(iV,:), ScoutFunc, ScoutOrient, nComponents, xyzFunction );
    % Update IK for next iteration (MEAN, other?)
    ImagingKernel2(i,:) = mean(sResults.ImagingKernel(iV,:));

end

% ===== FILL NEW ATLAS =====
sSurface = fill_new_atlas ( Labels2, sSurface, ['MLF:',num2str(LIM),'thres'], ScoutFunc );
% Save the modifications back to the surface file
bst_save(file_fullpath(sResults.SurfaceFile), sSurface, 'v6');

% Update for next iteration
% % if niter>1
% %     Keep = 1:nLabels2;
% % end

% % end

% ===== OUTPUT =====
OutputFiles{1} = sInput.FileName;

end
