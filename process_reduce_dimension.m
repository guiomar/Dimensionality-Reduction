function varargout = process_reduce_dimension( varargin )

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
% Authors: Guiomar Niso, Francois Tadel, 2014

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % ===== PROCESS =====
    % Description the process
    sProcess.Comment     = 'Reduce Dimensionality';
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
    % === N INITIAL SCOUTS
    sProcess.options.ninic.Comment = 'Number of initial regions';
    sProcess.options.ninic.Type    = 'value';
    sProcess.options.ninic.Value   = {50, '', 0};
    % === N ITERATIONS
    sProcess.options.niter.Comment = 'Number of iterations for MRA algorithm';
    sProcess.options.niter.Type    = 'value';
    sProcess.options.niter.Value   = {1, '', 0};
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
   
% == INPUTS ===

nIter  = sProcess.options.niter.Value{1};

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
        case {1, 'mean'}, ScoutFunc = 'mean'; % Most likely, 'mean'
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
%%%% isFlip??

% ---------------------------------------
nComponents = 1;
SplitFactor = 2;
%  ---------------------------------------
    
% ===== LOAD ALL INFO ===== 

% Load the surface filename from results file
sResults = in_bst_results(sInput.FileName, 0);
sResultsFull = in_bst_results(sInput.FileName, 1);  %%%%%%% Only for ImageGridAmp !!

SourceData = sResultsFull.ImageGridAmp; 
clear sResultsFull

% Get selected scouts and surface
[sScouts, AllAtlasNames, sSurface] = process_extract_scout('GetScoutsInfo', sProcess, sInput, sResults.SurfaceFile, AtlasList);

% Options for MRA algorithm
options.xyzFunction   = xyzFunction;
options.nComponents   = nComponents;
options.ScoutFunction = ScoutFunc;
options.SplitFactor   = SplitFactor;
options.nInicial      = sProcess.options.ninic.Value{1};

% MRA algorithm
nVertices  = size(sSurface.Vertices,1);
Labels1    = zeros(nVertices,1);
for i=1:nIter
    fprintf('BST> MRA iter #%i/%i \n', i, nIter);
    Labels1 = Labels1 + reduce_dimension_mra( SourceData, sSurface, sScouts, options );
end
LabelsF = round(Labels1 / nIter);
% uLabelsF = unique(LabelsF)'; % clusters
% nLabelsF = length(uLabelsF); % nClusters

% ===== FILL NEW ATLAS =====
sSurface = fill_new_atlas ( LabelsF, sSurface, ['MRA:',num2str(nIter),'iter'], options.ScoutFunction );

% Save the modifications back to the surface file
bst_save(file_fullpath(sResults.SurfaceFile), sSurface, 'v6');


% ===== OUTPUT =====
OutputFiles{1} = sInput.FileName;

end



