function [LayerMapFinal,Metrics,err] = reduce_dimension_mra( SourceData, sSurface, sScouts, options )

% REDUCE_DIMENSION_MRA: Multiresolution analysis of connectivity
%
% USAGE:  [LayerMapFinal,Metrics,err] = bst_connectivity_mra( SourceData, sSurface, sScouts, options )
%
% INPUTS: 
%    - SourceData  : N-source x N-samples time-series for analysis 
%    - Surface     : tesselation structure (Faces,Vertices,VertConn)
%    - sScouts     : Structure of scouts where to apply
%    - Options     : Options for the MRA algorithm
%      * nIniScouts: Initial number of seeds
%      * SplitFactor: factor to slipt the scouts on each iteration
%      * LabelFunction: function to group the signals of the vertices clustered
%      * nComponents
%      * xyzFunction
%
% OUTPUTS:
%    - LayerMapFinal : N-source vector with each number corresponding to
%                      the iteration at which they left the process (indicative
%                      of their probability of high connectivity)
%    - Metrics       : N-iteration cells containing connectivity metric for
%                      each iteration
%    - err           : Error message, if any
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
% Authors: Guiomar Niso, Sebastien Dery 2014


Verbose = 0;

% Start the search with random vertices each iteration
ShuffleVertices = 1; 

% Options for the MRA algorithm
nIniScouts    = options.nInicial;
ScoutFunction = options.ScoutFunction; %'mean'
nComponents   = options.nComponents;
xyzFunction   = options.xyzFunction;
SplitFactor   = options.SplitFactor; 

% -------------------------------------------
options.method = 'corr';
isAbsolute = 1; % for correlation, for example

type = 'cor'; % 'cor', 'pca' %%%%%%%%%%%%%%
% -------------------------------------------

% Vertices in the selected scouts
iVert = [sScouts.Vertices];

nVert_sel = length(iVert);
nVertices = length(sSurface.Vertices);

% Init varout
LayerMap = zeros(nVert_sel,1);
LayerMapFinal = zeros(nVertices,1);
Metrics = [];
err = [];


% ======== START =====================

nSamples = size(SourceData,2);
nDipoles = round( nVert_sel / nIniScouts);

% Correlation matrix of the data 
%%% PARA EL MAX DE TODOS INSTEAD OF THE MEAN
if strcmp(type,'cor')
    RemoveMean = 1;
    [R, pValues] = bst_corrn(SourceData, SourceData, RemoveMean);
    corF = R .* (pValues<0.05);
end

% ==== Initial Tessellation =====

% VertConn: spatial adjacency (connection) between vertices
VertConn = sSurface.VertConn(iVert,iVert);
% Scouts: contains the label of the scout each vertex belongs to
% (clusterized into nIniScouts)
Scouts = tess_cluster(VertConn, nIniScouts, ShuffleVertices, Verbose); 

iLayer = 1;

while (nDipoles > 2)

    oldScouts = Scouts; % Switch variable to prevent mixing

    % Tesselate the scout vertices ------------------------------------
    
    nameScouts = unique( nonzeros(oldScouts) );
    Label = 0; % Label: labels of the scouts

    for i=1:length(nameScouts)

        % iV: vertices that belonged to Scout i
        iV = find( oldScouts==nameScouts(i) );
        % (cluster into SplitFactor scouts, each of the previous scouts)
        Scouts(iV) = tess_cluster( sSurface.VertConn(iVert(iV),iVert(iV)), SplitFactor, ShuffleVertices, Verbose ) + Label; 
        Label = Label + SplitFactor; 

    end

    % After spliting --------------------------------------------------

    nameScouts = unique( nonzeros(Scouts) );
    nScouts = length( nameScouts ); %%% nChannels

    % N-scouts x N-sample time-series
    Fs = zeros(nScouts, nSamples);

    for i=1:nScouts
        % Get the scout vertices
        iV = find( Scouts == nameScouts(i) );
        % Get the scout orientation
        ScoutOrient = sSurface.VertNormals(iVert(iV),:);
        % Get scout value
        Fs(i,:) = bst_scout_value( SourceData(iVert(iV),:), ScoutFunction, ScoutOrient, nComponents, xyzFunction );
    end
    
    switch type
        case 'pca'
            % BEFORE: PCA --> COR
            RemoveMean = 1;
            [R, pValues] = bst_corrn(Fs, Fs, RemoveMean);
            R = R .* (pValues<0.05);

        case 'cor'
            % MAX OF ALL COR
            R = zeros(nScouts,nScouts);
            for i = 1:nScouts-1
                for j = i+1:nScouts
                    R(i,j) = max(max( corF( iVert(Scouts == nameScouts(i)), iVert(Scouts == nameScouts(j)) ) ));
                end
            end
            R = R + R';
    end

    % Remove diagonal
    R( eye(length(R))==1 ) = 0;

    if isAbsolute,  R = abs(R); end


    % ===== Regions to keep =====

    % METHOD:

    % Higher than mean - Rows 
    maxChanConn = max( R,[],2 );
    meanMaxConn = mean( maxChanConn );

    % CRITERIA: keep channels with greater connetivity than the mean
    Keep = maxChanConn > meanMaxConn;

    % If directional data (e.g. 'granger'), check the other way around
    if (strcmpi(options.method,'granger'))

        % Higher than mean - Columns
        maxChanConn = max(R);
        meanMaxConn = mean( maxChanConn );

        Keep = Keep | ( maxChanConn > meanMaxConn )';
    end

    
    % PROCEED: 

    % Remove everything else
    % Scouts with low connectivity
    removeScout = find( ~Keep );        
    if isempty(removeScout),   break;  end
    % Vertices belonging to that scout
    removeVertices = ismember(Scouts,removeScout);
    % Keep track of when these dipoles left the iteration process
    LayerMap( removeVertices==1 ) = iLayer; 
    % Leave these vertices out
    Scouts( removeVertices==1 ) = 0;

    % Keep metric matrix
    % Metrics{iLayer} = R;

    % UPDATE:
    nDipoles = floor( sum(Scouts~=0) / ( nScouts-length(removeScout) ) );
    iLayer = iLayer + 1;

end 

% Assign victorious last
LayerMap( LayerMap==0 ) = iLayer;

% Bring LayerMap to the final one 
LayerMapFinal(iVert) = LayerMap;

end

