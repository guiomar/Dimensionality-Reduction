function [Labels,err] = tess_cluster_ik (sSurface, ImagingKernel, vertices, threshold)
% TESS_CLUSTER_IK: Cluster vertices according to the IK
%
% USAGE:  [Labels,err] = tess_cluster_ik (sSurface, ImagingKernel, vertices, threshold)
%
% INPUTS: 
%    - sSurface        : tesselation structure (Faces,Vertices,VertConn)
%    - ImagingKernel   : N-source x N-channel sensor to source projection matrix
%    - vertices        : vertices to cluster
%    - threshold       : minimum value of correlation to cluster
%
% OUTPUTS:
%    - Labels : N-source vector with each number corresponding to an
%    agregation of source (0 means not assigned).
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
% Authors: Guiomar Niso 2014


% === CHEKING INPUT =========
if isempty(sSurface),       err = 'Surface structure is empty';     return; end
if isempty(ImagingKernel),  err = 'ImagingKernel matrix is empty';  return; end

% === START ===
err = [];

% == IMAGING KERNEL ==
nVertices = length (vertices);
MAXs = zeros(nVertices,nVertices);

% Select only the number of vertices
ImagingKernel = ImagingKernel(vertices,:);

% Normalize IK
normIK = sqrt( sum( abs(ImagingKernel).^2,2 ) );
IK = ImagingKernel ./ repmat( normIK, 1, size(ImagingKernel,2) );


% == SPATIAL CONNECTIVITY ==
SPC = sSurface.VertConn (vertices, vertices);    
SPC( eye( size(SPC) )== 1 ) = 0; % Remove diagonal

% == CORRELATION ==
% Compute correlation
COR = corrcoef(IK');


% == MAX ==

%%% OPTION 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Every iteration finds new neighbours and adds the maxs

% % Maximum values of the neighbours
% [vmax, imax] = max(COR.*SPC);
% % Create the  corresponding matrix
% for i=1:nVertices
%     MAXs(imax(i),i)=1;
%     MAXs(i,imax(i))=1;
% end

%%% OPTION 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Every iteration finds new neighbours but only add maxs if also cor>LIM

% % Maximum values of the neighbours
% [vmax, imax] = max(COR.*SPC);
% % LIMIT the minimum value of correlation that should be reached to cluster
% imax( vmax<threshold ) = 0;
% % Create the  corresponding matrix
% for i=1:nVertices
%     if imax(i)>0
%         MAXs(imax(i),i)=1;
%         MAXs(i,imax(i))=1;
%     end
% end

%%% OPTION 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cluster all at once
% Only needs one iteration because it will cluster all that verify cor>LIM

A = COR.*SPC;
MAXs = A>threshold;


% Set the labels
Labels = tess_group_binary (MAXs);
            
end

