function varargout = process_random_scouts_generator( varargin )
% process_random_scouts_generator: Generates a fixed number of randomly 
% located brain scouts with certain number of vertices that also vary  
% randomly into the specified range (min max)
%
% USAGE: varargout = process_random_scouts_generator( varargin )  
 
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
%
% It calls: fill_new_atlas.m

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Generate random brain scouts';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Simulate'; 
    sProcess.Index       = 907; 
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'import'}; % No input
    sProcess.OutputTypes = {'import'}; % No output
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 0;
    sProcess.isSeparator = 0;

    % === SUBJECT NAME
    sProcess.options.subjectname.Comment = 'Subject name:';
    sProcess.options.subjectname.Type    = 'subjectname';
    sProcess.options.subjectname.Value   = 'Test';
    % === CONDITION NAME
    sProcess.options.condition.Comment = 'Condition name:';
    sProcess.options.condition.Type    = 'text';
    sProcess.options.condition.Value   = 'Simulation';
    % === NUMBER OF SCOUTS
    sProcess.options.nscouts.Comment = 'Number of scouts:';
    sProcess.options.nscouts.Type    = 'value';
    sProcess.options.nscouts.Value   = {5, '', 0};
    % === NUMBER OF VERTICES
    sProcess.options.nvertices.Comment = 'Number of vertices (min,max):';
    sProcess.options.nvertices.Type    = 'range';
    sProcess.options.nvertices.Value   = {[10 100], '', 0};
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput) %#ok<DEFNU>
    
    OutputFiles = {'import'};
    
    % === GET OPTIONS ===
    % Get subject name
    SubjectName = file_standardize(sProcess.options.subjectname.Value);
    if isempty(SubjectName)
        bst_report('Error', sProcess, sInputs, 'Subject name is empty.');
        return
    end
    % Get condition name
    Condition = file_standardize(sProcess.options.condition.Value);
    % Get number of scouts
    nscouts = sProcess.options.nscouts.Value{1};
    % Get number of vertices
    rangevert = sProcess.options.nvertices.Value{1:2};
    
    % === LOAD CORTEX SURFACE ===
    % Get subject
    sSubject = bst_get('Subject', SubjectName);
    if isempty(sSubject)
        bst_report('Error', sProcess, sInputs, 'Subject doesn''t exist.');
        return
    end
    if isempty(sSubject.Surface) || isempty(sSubject.iCortex)
        bst_report('Error', sProcess, sInputs, 'Subject has no default cortex surface.');
        return
    end
    % Get surface file
    CortexFile = sSubject.Surface(sSubject.iCortex).FileName;
    % Load surface
    sSurface = in_tess_bst(CortexFile);
    

    % ----
    LabelFunction = 'mean';
    % ----
    
    % === OUTPUT CONDITION ===
    % Default condition name
    if isempty(Condition)
        Condition = 'Simulation';
    end
    % Get condition asked by user
    [sStudy, iStudy] = bst_get('StudyWithCondition', bst_fullfile(SubjectName, Condition));
    % Condition does not exist: create it
    if isempty(sStudy)
        iStudy = db_add_condition(SubjectName, Condition, 1);
        sStudy = bst_get('Study', iStudy);
    end
    
    
    % ===== GENERATE SCOUTS =====
    
    % Vertices to classify
    N = size(sSurface.VertConn,1);
    iVertLeft = 1:N;
    
    % Labels
    Labels = zeros(1,N);
    
    for i=1:nscouts
        
        % Number of vertices of the scout
        if rangevert(2)>rangevert(1)
            nvertices = rangevert(1) + random('unid',rangevert(2)-rangevert(1));
        else % Fixed number of vertices (min=max)
            nvertices = rangevert(1);
        end
        
        % Start scout with the first vertex in the list
        iScout = iVertLeft( random('unid',length(iVertLeft)) );
        iNewVert = iScout;

        % Grow region until it's not growing anymore
        while ~isempty(iNewVert)
            
            iNewVert = tess_scout_swell(iScout, sSurface.VertConn);
            iScout1 = union(iScout, iNewVert);
            
            % Check it doesn't exceed the max number of vertices
            if length(iScout1) >=  rangevert(2)
                if length(iScout) <  rangevert(1)
                    iScout = iScout1(1:rangevert(1)); % It can give some isolated, but it's ok, we need to reach the minimum
                end
                break;
            end
            iScout = iScout1;
            
            % Stop when the desired number of vertices is reached 
            if length(iScout1) >=  nvertices
                break;
            end

        end

        Labels(iScout) = i;
        iVertLeft = setdiff(iVertLeft, iScout);

    end


    % ===== FILL NEW ATLAS =====
    sSurface = fill_new_atlas ( Labels, sSurface, 'RAND', LabelFunction );
    % Save the modifications back to the surface file
    bst_save(file_fullpath(CortexFile), sSurface, 'v6');
    
    
    % WE SHOULD UNLOAD THE MODIFIED SURFACE FILE
    


end 

