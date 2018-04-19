function sSurface = fill_new_atlas ( Labels, sSurface, Name, LabelFunction )

% FILL_NEW_ATLAS: Creates a Surface including a new atlas according to the labels
%
% USAGE:  sSurface = fill_new_atlas ( Labels, sSurface, Name, LabelFunction )
%
% INPUTS: 
%    - Labels        : Labels you want to assign 
%    - sSurface      : tesselation structure (Faces,Vertices,VertConn)
%    - Name          : Name for the atlas
%    - LabelFunction : function to group the signals of the vertices clustered
%
% OUTPUTS:
%    - sSurface : New Surface including the new atlas
%
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
    

% Labels
uLabels = unique(Labels); % clusters
uLabels(uLabels==0)=[];    % remove label==0
nLabels = length(uLabels); % nClusters


% ===== FILL NEW ATLAS =====

% Create a new atlas
sNewAtlas = db_template('Atlas');
sNewAtlas.Name = Name;


% String format depends on the number of clusters
if (nLabels > 99), strFormat = '%03d';
else               strFormat = '%02d';     end

% Ordered colors
if nLabels <= 64,   SurfaceColor = jet(length(uLabels));
else                SurfaceColor = rand(nLabels,3);     end

for iScout = 1:nLabels
    sNewAtlas.Scouts(iScout).Vertices = find(Labels == uLabels(iScout));
    sNewAtlas.Scouts(iScout).Seed     = sNewAtlas.Scouts(iScout).Vertices(1); 
    sNewAtlas.Scouts(iScout).Label    = sprintf(strFormat, iScout);
    sNewAtlas.Scouts(iScout).Color    = SurfaceColor(iScout,:);
    sNewAtlas.Scouts(iScout).Function = LabelFunction;
end

% Fix atlas structure
sNewAtlas = panel_scout('FixAtlasStruct', sNewAtlas);
% Add new atlas to an existing list
if ~isempty(sSurface.Atlas)
    % Create a unique atlas name
    sNewAtlas.Name = file_unique(sNewAtlas.Name, {sSurface.Atlas.Name});
    % Add atlas to surface structure
    sSurface.Atlas(end+1) = sNewAtlas;
else
    sSurface.Atlas = sNewAtlas;
end


end

