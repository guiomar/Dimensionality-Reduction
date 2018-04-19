function Labels = tess_group_binary (VertConn)

% TESS_GROUP_BINARY: Group the connected vertices of a binnary matrix
%
% USAGE:  Labels = tess_remove_small(VertConn)
% 
% INPUT:  VerConn = binary matrix: nVerticces x nVertices
%
% OUTPUT: Labels = name of the cluster each vertice belong to
% 
% Guiomar Niso, 2014

nVertices = size(VertConn,1);
Labels = zeros (1,nVertices);

lab = 1;

% Vertices to classify
iVertLeft = 1:nVertices;

while ~isempty(iVertLeft)
    
    % Start scout with the first vertex in the list
    iScout = iVertLeft(1);
    iNewVert = iScout;
    
    % Grow region until it's not growing anymore
    while ~isempty(iNewVert)
        iScout = union(iScout, iNewVert);
        iNewVert = tess_scout_swell(iScout, VertConn);
    end
    
    Labels(iScout) = lab;
    lab = lab+1;
    iVertLeft = setdiff(iVertLeft, iScout);
    
end


