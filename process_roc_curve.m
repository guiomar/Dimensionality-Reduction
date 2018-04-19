function varargout = process_roc_curve( varargin )
% process_roc_curve: Generates ROC curve
%
% USAGE:   
 
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
    % Description the process
    sProcess.Comment     = 'ROC curve';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Simulate'; 
    sProcess.Index       = 908; 
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;

    % === SUBJECT NAME
    sProcess.options.subjectname.Comment = 'Subject name:';
    sProcess.options.subjectname.Type    = 'subjectname';
    sProcess.options.subjectname.Value   = 'Test';
    % === CONDITION NAME
    sProcess.options.condition.Comment = 'Condition name:';
    sProcess.options.condition.Type    = 'text';
    sProcess.options.condition.Value   = 'Simulation';
    % === GROUND TRUTH
    sProcess.options.GT.Comment = 'Ground Truth:';
    sProcess.options.GT.Type    = 'atlas';
    sProcess.options.GT.Value   = [];
    % === RESULTS TO COMPARE 1
    sProcess.options.RES1.Comment = 'Results1 to compare:';
    sProcess.options.RES1.Type    = 'atlas';
    sProcess.options.RES1.Value   = [];
%     % === RESULTS TO COMPARE 2
%     sProcess.options.RES2.Comment = 'Results2 to compare:';
%     sProcess.options.RES2.Type    = 'atlas';
%     sProcess.options.RES2.Value   = [];
    % === NUMBER OF POINTS FOR THE ROC CURVE
%     sProcess.options.Nt.Comment = 'Number of point for the ROC curve:';
%     sProcess.options.Nt.Type    = 'value';
%     sProcess.options.Nt.Value   = {10,[],0}; %
    
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput) %#ok<DEFNU>
    
    OutputFiles = {};
    
    % === GET OPTIONS ===
    % Get subject name
    SubjectName = file_standardize(sProcess.options.subjectname.Value);
    if isempty(SubjectName)
        bst_report('Error', sProcess, sInputs, 'Subject name is empty.');
        return
    end
    % Get condition name
    Condition = file_standardize(sProcess.options.condition.Value);
    % Get number of vertices for the ROC curve
%     Nt = sProcess.options.Nt.Value{1};
    
   
    % Load the surface filename from results file
    ResultsMat = in_bst_results(sInput.FileName, 0);
    % Load surface
    SurfaceMat = in_tess_bst(ResultsMat.SurfaceFile);
    
    % Get the atlas to use as GT
    iAtlasGT = [];
    if ~isempty(sProcess.options.GT.Value)
        iAtlasGT = find(strcmpi({SurfaceMat.Atlas.Name}, sProcess.options.GT.Value));
    end
    if isempty(iAtlasGT),   iAtlasGT = SurfaceMat.iAtlas;    end
    if isempty(iAtlasGT),   iAtlasGT = 1;    end
    
    % Get the atlas to use as RESULTS 1
    iAtlasRES1 = [];
    if ~isempty(sProcess.options.RES1.Value)
        iAtlasRES1 = find(strcmpi({SurfaceMat.Atlas.Name}, sProcess.options.RES1.Value));
    end
    if isempty(iAtlasRES1),  iAtlasRES1 = SurfaceMat.iAtlas;    end
    if isempty(iAtlasRES1),  iAtlasRES1 = 1;     end

    % Get the atlas to use as RESULTS 2
%     iAtlasRES2 = [];
%     if ~isempty(sProcess.options.RES2.Value)
%         iAtlasRES2 = find(strcmpi({SurfaceMat.Atlas.Name}, sProcess.options.RES2.Value));
%     end
%     if isempty(iAtlasRES2),  iAtlasRES2 = SurfaceMat.iAtlas;    end
%     if isempty(iAtlasRES2),  iAtlasRES2 = 1;     end

    % Convert to matrices
    Nvertices = size(SurfaceMat.Vertices,1);
    GT   = zeros(1, Nvertices);
    RES1 = zeros(1, Nvertices);
%     RES2 = zeros(1, Nvertices);

    for i=1:numel(SurfaceMat.Atlas(iAtlasGT).Scouts)
        GT( 1, SurfaceMat.Atlas(iAtlasGT).Scouts(i).Vertices ) = 1; % i: si queremos conservar el grupo
    end    
    for i=1:numel(SurfaceMat.Atlas(iAtlasRES1).Scouts)
        RES1( 1, SurfaceMat.Atlas(iAtlasRES1).Scouts(i).Vertices ) = i;
    end    
%     for i=1:numel(SurfaceMat.Atlas(iAtlasRES2).Scouts)
%         RES2(SurfaceMat.Atlas(iAtlasRES2).Scouts(i).Vertices)=i;
%     end    

   
    
    % ===== GET ROC CURVES =====
    % Get ROC curve
    [TPR1, FPR1, thres1] = check_roc(GT, RES1);
%     [TPR2, FPR2, thres2] = check_roc(GT, RES2);
    
    % Get intersection point (*function from Matlab central (intersections)****)
%     [Cx,Cy] = intersections(FPR,TPR,linspace(1,0,Nt),linspace(0,1,Nt));
    
    % Get AUC (Area Under the Curve) !!!!
    
    % Get relative error !!!!
    
    % Plot ROC curve ------------------------------------------------------
    COLOR = {'b','g'}; %%%%% SEVERAL CURVES!!!!
    markerSize = 30;
    lineWidth = 3;
    
    figure()
    hold on;
    plot( [ 0 FPR1 ], [ 0 TPR1 ], '.-','color',COLOR{1},'MarkerSize',markerSize,'LineWidth', lineWidth);% 'MarkerFaceColor',
%     plot( [ 0 FPR2 ], [ 0 TPR2 ], '.-','color',COLOR{2},'MarkerSize',markerSize,'LineWidth', lineWidth);% 'MarkerFaceColor',
    plot( [ 0 1 ], [ 0 1 ], 'k:');
    xlabel('FPR');
    ylabel('TPR');
%     axis([0 1 0 1]); 
    axis square;
    legend(SurfaceMat.Atlas(iAtlasRES1).Name,'random');
    
    

    % ===== OUTPUT =====
    OutputFiles{1} = sInput.FileName;
    
    %%%% REMOVE!!!!!!!!
    roc_result = [ TPR1; FPR1 ];
    save roc_result roc_result

end 

